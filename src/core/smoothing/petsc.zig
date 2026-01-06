// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const petsc = @import("petsc_import.zig");

const RowCompressedMatrixSystem2d = @import("smooth.zig").RowCompressedMatrixSystem2d;

pub const PetscSolver = struct {
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

test "petsc simple" {
    if (petsc.PetscInitialize(null, null, null, null) != petsc.PETSC_SUCCESS) {
        return error.PetscInitializeFailed;
    }

    // Create Matrix

    var A: petsc.Mat = undefined;

    if (petsc.MatCreate(petsc.PETSC_COMM_WORLD, &A) != petsc.PETSC_SUCCESS) {
        return error.MatCreateFailed;
    }

    const m: petsc.PetscInt = 2;
    const n: petsc.PetscInt = 2;
    _ = petsc.MatSetSizes(A, petsc.PETSC_DECIDE, petsc.PETSC_DECIDE, m, n);

    // matrix format (?)

    _ = petsc.MatSetFromOptions(A);

    // TODO: prealloc!

    _ = petsc.MatSetUp(A);

    // Set Matrix A = diag(1, 2)
    const rows = [_]petsc.PetscInt{ 0, 1 };
    const cols = [_]petsc.PetscInt{ 0, 1 };
    const vals = [_]petsc.PetscScalar{ 1.0, 2.0 };
    _ = petsc.MatSetValues(A, 1, &rows[0], 1, &cols[0], &vals[0], petsc.INSERT_VALUES);
    _ = petsc.MatSetValues(A, 1, &rows[1], 1, &cols[1], &vals[1], petsc.INSERT_VALUES);

    _ = petsc.MatAssemblyBegin(A, petsc.MAT_FINAL_ASSEMBLY);
    _ = petsc.MatAssemblyEnd(A, petsc.MAT_FINAL_ASSEMBLY);

    // Create RHS

    var b: petsc.Vec = undefined;
    _ = petsc.VecCreate(petsc.PETSC_COMM_WORLD, &b);
    _ = petsc.VecSetSizes(b, petsc.PETSC_DECIDE, m);
    _ = petsc.VecSetFromOptions(b);

    // Set RHS b = [1, 4]^T
    const rhs = [_]petsc.PetscScalar{ 1, 4 };
    _ = petsc.VecSetValues(b, 2, rows[0..], rhs[0..], petsc.INSERT_VALUES); // TODO: use better version

    _ = petsc.VecAssemblyBegin(b);
    _ = petsc.VecAssemblyEnd(b);

    // Solver

    var ksp: petsc.KSP = undefined;
    _ = petsc.KSPCreate(petsc.PETSC_COMM_WORLD, &ksp);
    _ = petsc.KSPSetOperators(ksp, A, A);

    _ = petsc.KSPSolve(ksp, b, b);

    var res: [*c]const petsc.PetscScalar = undefined;
    _ = petsc.VecGetArrayRead(b, &res);

    try std.testing.expectApproxEqAbs(1.0, res[0], 1e-15);
    try std.testing.expectApproxEqAbs(2.0, res[1], 1e-15);

    _ = petsc.KSPDestroy(&ksp);
    _ = petsc.VecDestroy(&b);
    _ = petsc.MatDestroy(&A);

    _ = petsc.PetscFinalize();
}
