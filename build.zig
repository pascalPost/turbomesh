// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

const core = @import("build/core.zig");
const tests = @import("build/tests.zig");
const gui = @import("build/gui.zig");
const wasm = @import("build/wasm.zig");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // CLI options
    const use_cgns = b.option(bool, "use-cgns", "Enable cgns");
    const use_umfpack = b.option(bool, "use-umfpack", "Enable umfpack");
    const use_petsc = b.option(bool, "use-petsc", "Enable petsc");

    tests.addTests(b, target, optimize, .{
        .use_cgns = use_cgns orelse true,
        .use_umfpack = use_umfpack orelse true,
        .use_petsc = use_petsc orelse true,
    });

    gui.addDesktopGui(b, target, optimize, .{
        .use_cgns = use_cgns orelse true,
        .use_umfpack = use_umfpack orelse true,
        .use_petsc = use_petsc orelse true,
    });

    wasm.addWasm(b, target, optimize, .{
        .use_cgns = use_cgns orelse false,
        .use_umfpack = use_umfpack orelse false,
        .use_petsc = use_petsc orelse false,
    });
}
