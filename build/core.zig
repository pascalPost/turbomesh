// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

pub fn addCoreModule(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
) *std.Build.Module {
    const mod = b.addModule("core", .{
        .root_source_file = b.path("src/core/lib.zig"),
        .target = target,
        .optimize = optimize,
    });

    mod.linkSystemLibrary("cgns", .{});
    mod.linkSystemLibrary("umfpack", .{});
    mod.linkSystemLibrary("petsc", .{});
    mod.addCSourceFile(.{ .file = b.path("src/core/smoothing/petsc_shim.c") });
    mod.linkSystemLibrary("mpi", .{}); // needed by petsc

    mod.link_libc = true;

    // needed on macos
    mod.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    mod.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    const tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/core/tests.zig"),
            .target = target,
            .optimize = optimize,
        }),
        // needed for debugging for now, as the new zig based debug compiler does not yet provide debug symbols
        .use_llvm = true,
    });

    tests.linkSystemLibrary("cgns");
    tests.linkSystemLibrary("umfpack");
    tests.linkSystemLibrary("petsc");
    tests.addCSourceFile(.{ .file = b.path("src/core/smoothing/petsc_shim.c") });
    tests.linkSystemLibrary("mpi"); // needed by petsc
    tests.linkLibC();

    // needed on macos
    tests.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    tests.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    const run_tests = b.addRunArtifact(tests);
    run_tests.has_side_effects = true;

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);

    return mod;
}
