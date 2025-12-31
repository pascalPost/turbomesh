// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

pub fn addDesktopGui(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    core_mod: *std.Build.Module,
) void {
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
            .root_source_file = b.path("src/gui/lib.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .use_llvm = true,
    });

    lib.root_module.addImport("core", core_mod);
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
            .root_source_file = b.path("src/gui/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .use_llvm = true,
    });

    exe.root_module.addImport("core", core_mod);
    exe.root_module.addImport("zglfw", zglfw.module("root"));

    if (target.result.os.tag != .emscripten) {
        exe.linkLibrary(zglfw.artifact("glfw"));
    }

    exe.root_module.addImport("gl", gl_bindings);

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep()); // run from install directory
    if (b.args) |args| // pass arguments to the executable
    {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run the desktop GUI");
    run_step.dependOn(&run_cmd.step);
}
