// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});

    const gl_bindings = @import("zigglgen").generateBindingsModule(b, .{
        .api = .gl,
        .version = .@"4.5",
        .profile = .core,
        .extensions = &.{},
    });

    const zglfw = b.dependency("zglfw", .{
        .target = target,
        .optimize = optimize,
        .shared = true,
    });

    const lib = b.addLibrary(.{
        .name = "turbomesh",
        .linkage = .dynamic,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/lib.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .use_llvm = true,
    });

    lib.root_module.addImport("gl", gl_bindings);
    lib.root_module.addImport("zglfw", zglfw.module("root"));

    if (target.result.os.tag != .emscripten) {
        lib.linkLibrary(zglfw.artifact("glfw"));
    }

    b.installArtifact(lib);

    const exe = b.addExecutable(.{
        .name = "turbomesh",
        .version = .{ .major = 0, .minor = 1, .patch = 0 },
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .use_llvm = true,
    });

    exe.root_module.addImport("zglfw", zglfw.module("root"));

    if (target.result.os.tag != .emscripten) {
        exe.linkLibrary(zglfw.artifact("glfw"));
    }

    exe.root_module.addImport("gl", gl_bindings);

    exe.linkSystemLibrary("cgns");
    exe.linkSystemLibrary("umfpack");
    exe.linkSystemLibrary("petsc");
    exe.addCSourceFile(.{ .file = b.path("src/smoothing/petsc_shim.c") });
    exe.linkSystemLibrary("mpi"); // needed by petsc
    exe.linkLibC();

    // needed on macos
    exe.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    exe.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    // This declares intent for the executable to be installed into the
    // standard location when the user invokes the "install" step (the default
    // step when running `zig build`).
    b.installArtifact(exe);

    // This *creates* a Run step in the build graph, to be executed when another
    // step is evaluated that depends on it. The next line below will establish
    // such a dependency.
    const run_cmd = b.addRunArtifact(exe);

    // By making the run step depend on the install step, it will be run from the
    // installation directory rather than directly from within the cache directory.
    // This is not necessary, however, if the application depends on other installed
    // files, this ensures they will be present and in the expected location.
    run_cmd.step.dependOn(b.getInstallStep());

    // This allows the user to pass arguments to the application in the build
    // command itself, like this: `zig build run -- arg1 arg2 etc`
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // This creates a build step. It will be visible in the `zig build --help` menu,
    // and can be selected like this: `zig build run`
    // This will evaluate the `run` step rather than the default, which is "install".
    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    const tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/tests.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .use_llvm = true,
    });

    tests.linkSystemLibrary("cgns");
    tests.linkSystemLibrary("umfpack");
    tests.linkSystemLibrary("petsc");
    tests.addCSourceFile(.{ .file = b.path("src/smoothing/petsc_shim.c") });
    tests.linkSystemLibrary("mpi"); // needed by petsc
    tests.linkLibC();

    // needed on macos
    tests.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    tests.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    const run_tests = b.addRunArtifact(tests);
    run_tests.has_side_effects = true;

    // Similar to creating the run step earlier, this exposes a `test` step to
    // the `zig build --help` menu, providing a way for the user to request
    // running the unit tests.
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);
}
