// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// TODO: allow to compile with a subset of solvers

const std = @import("std");
const umfpack = @import("umfpack.zig");
const petsc = @import("petsc_import.zig");
const RowCompressedMatrixSystem2d = @import("smooth.zig").RowCompressedMatrixSystem2d;

pub const Type = enum {
    umfpack,
    petsc,
};

pub const Solver = union(Type) {
    umfpack: UmfpackSolver,
    petsc: PetscSolver,

    pub fn init(backend: Type, system: RowCompressedMatrixSystem2d) Solver {
        switch (backend) {
            .umfpack => return .{ .umfpack = UmfpackSolver.init(system) },
            .petsc => return .{ .petsc = PetscSolver.init(system) },
        }
    }

    pub fn deinit(self: *Solver) void {
        switch (self.*) {
            .umfpack => {},
            .petsc => |*s| s.deinit(),
        }
    }

    pub fn solve(self: *Solver) !void {
        try switch (self.*) {
            .umfpack => |s| try s.solve(),
            .petsc => |*s| s.solve(),
        };
    }
};

const UmfpackSolver = struct {
    system: RowCompressedMatrixSystem2d,

    pub fn init(system: RowCompressedMatrixSystem2d) UmfpackSolver {
        return .{ .system = system };
    }

    pub fn solve(self: UmfpackSolver) !void {
        const dof = self.system.rhs_x.len;
        try umfpack.solve(@intCast(dof), @intCast(dof), self.system.lhs_p, self.system.lhs_i, self.system.lhs_values, self.system.rhs_x, self.system.x_new, self.system.rhs_y, self.system.y_new);
    }
};

const PetscSolver = struct {
    system: RowCompressedMatrixSystem2d,
    ksp: petsc.KSP,
    A: petsc.Mat,
    rhs_x: petsc.Vec,
    rhs_y: petsc.Vec,
    x_new: petsc.Vec,
    y_new: petsc.Vec,

    pub fn init(system: RowCompressedMatrixSystem2d) PetscSolver {
        _ = petsc.PetscInitialize(null, null, null, null);

        // // Current configuration: CG with algebraic multigrid
        // _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type cg -pc_type gamg -pc_gamg_type agg -pc_gamg_agg_nsmooths 1 -pc_gamg_threshold 0.01 -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        // Alternative configurations to try if the above doesn't work:
        // 1. GMRES with ILU preconditioner:
        // _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type gmres -pc_type ilu -pc_factor_levels 2 -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        // 2. CG with Jacobi preconditioner (simple but often works):
        // _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type cg -pc_type jacobi -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        // 3. GMRES with SOR preconditioner:
        // _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type gmres -pc_type sor -pc_sor_omega 1.5 -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        // 4. BiCGStab with ILU (good for non-symmetric problems):
        // _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type bcgs -pc_type ilu -pc_factor_levels 1 -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        _ = petsc.PetscOptionsInsertString(null, "-ksp_monitor -ksp_converged_reason -ksp_view -ksp_type bcgs -pc_type sor -pc_sor_omega 1.5 -ksp_rtol 1e-6 -ksp_atol 1e-8 -ksp_max_it 1000");

        const dof = system.rhs_x.len;
        const num_rows: petsc.PetscInt = @intCast(dof);
        const num_cols: petsc.PetscInt = @intCast(dof);
        const rows = system.lhs_p;
        const cols = system.lhs_i;
        const values = system.lhs_values;

        var s: PetscSolver = undefined;
        s.system = system;

        const PETSC_COMM_SELF = petsc.get_petsc_comm_self();

        _ = petsc.MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, num_rows, num_cols, rows.ptr, cols.ptr, values.ptr, &s.A);

        _ = petsc.KSPCreate(PETSC_COMM_SELF, &s.ksp);
        _ = petsc.KSPSetOperators(s.ksp, s.A, s.A);
        _ = petsc.KSPSetFromOptions(s.ksp);

        // Additional solver settings for better convergence
        // _ = petsc.KSPSetTolerances(ksp, 1e-6, 1e-8, 1e-10, 1000);
        _ = petsc.KSPSetInitialGuessNonzero(s.ksp, true); // Use previous solution as initial guess

        _ = petsc.KSPSetUp(s.ksp);

        _ = petsc.VecCreateSeqWithArray(PETSC_COMM_SELF, 1, num_rows, system.rhs_x.ptr, &s.rhs_x);
        _ = petsc.VecCreateSeqWithArray(PETSC_COMM_SELF, 1, num_rows, system.x_new.ptr, &s.x_new);

        _ = petsc.VecCreateSeqWithArray(PETSC_COMM_SELF, 1, num_rows, system.rhs_y.ptr, &s.rhs_y);
        _ = petsc.VecCreateSeqWithArray(PETSC_COMM_SELF, 1, num_rows, system.y_new.ptr, &s.y_new);

        return s;
    }

    pub fn deinit(self: *PetscSolver) void {
        _ = petsc.KSPDestroy(&self.ksp);
        _ = petsc.MatDestroy(&self.A);
        _ = petsc.VecDestroy(&self.rhs_x);
        _ = petsc.VecDestroy(&self.rhs_y);
        _ = petsc.VecDestroy(&self.x_new);
        _ = petsc.VecDestroy(&self.y_new);
        _ = petsc.PetscFinalize();
    }

    pub fn solve(self: *PetscSolver) !void {
        _ = petsc.KSPSetOperators(self.ksp, self.A, self.A);
        _ = petsc.KSPSolve(self.ksp, self.rhs_x, self.x_new);
        _ = petsc.KSPSolve(self.ksp, self.rhs_y, self.y_new);

        // Check convergence reason
        var reason: petsc.KSPConvergedReason = undefined;
        _ = petsc.KSPGetConvergedReason(self.ksp, &reason);
        if (reason < 0) {
            std.debug.print("Warning: KSP did not converge. Reason: {}\n", .{reason});
        }

        _ = petsc.MatDestroy(&self.A);
        const dof = self.system.rhs_x.len;
        const num_rows: petsc.PetscInt = @intCast(dof);
        const num_cols: petsc.PetscInt = @intCast(dof);
        const rows = self.system.lhs_p;
        const cols = self.system.lhs_i;
        const values = self.system.lhs_values;
        _ = petsc.MatCreateSeqAIJWithArrays(petsc.get_petsc_comm_self(), num_rows, num_cols, rows.ptr, cols.ptr, values.ptr, &self.A);
    }
};
