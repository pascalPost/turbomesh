const std = @import("std");

const core = @import("core.zig");

pub fn addTests(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    cfg: core.Config,
) void {
    const tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/core/tests.zig"),
            .target = target,
            .optimize = optimize,
        }),
        // needed for debugging for now, as the new zig based debug compiler does not yet provide debug symbols
        .use_llvm = true,
    });

    const opts = b.addOptions();
    opts.addOption(bool, "use_cgns", cfg.use_cgns);
    opts.addOption(bool, "use_umfpack", cfg.use_umfpack);
    opts.addOption(bool, "use_petsc", cfg.use_petsc);
    tests.root_module.addOptions("config", opts);

    if (cfg.use_cgns) {
        tests.linkSystemLibrary("cgns");
        tests.linkLibC();
    }
    if (cfg.use_umfpack) {
        tests.linkSystemLibrary("umfpack");
        tests.linkLibC();
    }
    if (cfg.use_petsc) {
        tests.linkSystemLibrary("petsc");
        tests.addCSourceFile(.{ .file = b.path("src/core/smoothing/petsc_shim.c") });
        tests.linkSystemLibrary("mpi"); // needed by petsc
        tests.linkLibC();
    }

    // needed on macos
    tests.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    tests.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    const run_tests = b.addRunArtifact(tests);
    run_tests.has_side_effects = true;

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);
}
