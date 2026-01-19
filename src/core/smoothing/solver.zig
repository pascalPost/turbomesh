// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const config = @import("config");

const Preconditioner = @import("preconditioner.zig").Preconditioner;
const RowCompressedMatrixSystem2d = @import("smooth.zig").RowCompressedMatrixSystem2d;

pub const Tag = enum {
    gmres,
    bicgstab,
    umfpack,
    petsc,
};

// TODO: move this into the input as this is used for parsing.
pub const Option = union(Tag) {
    gmres: struct {
        preconditioner: Preconditioner,
    },
    bicgstab: struct {
        preconditioner: Preconditioner,
    },
    umfpack: struct {},
    petsc: struct {},
};

pub const Solver = union(Tag) {
    const bicgstab_impl = @import("BiCGStab.zig");
    const gmres_impl = @import("GMRES.zig");
    const UmfpackSolver = if (config.use_umfpack) @import("umfpack.zig").UmfpackSolver else void;
    const PetscSolver = if (config.use_petsc) @import("petsc.zig").PetscSolver else void;

    gmres: gmres_impl.GMRESSolver,
    bicgstab: bicgstab_impl.BiCGStabSolver,
    umfpack: UmfpackSolver,
    petsc: PetscSolver,

    pub fn init(option: Option, system: RowCompressedMatrixSystem2d) !Solver {
        switch (option) {
            .gmres => |s| return .{ .gmres = gmres_impl.GMRESSolver.init(system, s.preconditioner) },
            .bicgstab => |s| return .{ .bicgstab = bicgstab_impl.BiCGStabSolver.init(system, s.preconditioner) },
            .umfpack => {
                if (config.use_umfpack) {
                    return .{ .umfpack = UmfpackSolver.init(system) };
                } else {
                    return error.ExternalSolverNotEnabled;
                }
            },
            .petsc => {
                if (config.use_petsc) {
                    return .{ .petsc = PetscSolver.init(system) };
                } else {
                    return error.ExternalSolverNotEnabled;
                }
            },
        }
    }

    pub fn deinit(self: *Solver) void {
        switch (self.*) {
            .bicgstab => |*s| s.deinit(),
            .gmres => |*s| s.deinit(),
            .umfpack => {},
            .petsc => |*s| {
                if (config.use_petsc) {
                    s.deinit();
                }
            },
        }
    }

    pub fn solve(self: *Solver) !void {
        switch (self.*) {
            .bicgstab => |*s| try s.solve(),
            .gmres => |*s| try s.solve(),
            .umfpack => |s| {
                if (config.use_umfpack) {
                    try s.solve();
                } else {
                    return error.ExternalSolverNotEnabled;
                }
            },
            .petsc => |*s| {
                if (config.use_petsc) {
                    try s.solve();
                } else {
                    return error.ExternalSolverNotEnabled;
                }
            },
        }
    }
};
