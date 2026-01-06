// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

pub const Config = struct {
    use_cgns: bool,
    use_umfpack: bool,
    use_petsc: bool,
};

pub fn addCoreModule(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    cfg: Config,
) *std.Build.Module {
    const mod = b.addModule("core", .{
        .root_source_file = b.path("src/core/lib.zig"),
        .target = target,
        .optimize = optimize,
    });

    const opts = b.addOptions();
    opts.addOption(bool, "use_cgns", cfg.use_cgns);
    opts.addOption(bool, "use_umfpack", cfg.use_umfpack);
    opts.addOption(bool, "use_petsc", cfg.use_petsc);
    mod.addOptions("config", opts);

    if (cfg.use_cgns) {
        mod.linkSystemLibrary("cgns", .{});
        mod.link_libc = true;
    }
    if (cfg.use_umfpack) {
        mod.linkSystemLibrary("umfpack", .{});
        mod.link_libc = true;
    }
    if (cfg.use_petsc) {
        mod.linkSystemLibrary("petsc", .{});
        mod.addCSourceFile(.{ .file = b.path("src/core/smoothing/petsc_shim.c") });
        mod.linkSystemLibrary("mpi", .{}); // needed by petsc
        mod.link_libc = true;
    }

    // needed on macos
    mod.addLibraryPath(.{ .cwd_relative = "/usr/local/lib/" });
    mod.addSystemIncludePath(.{ .cwd_relative = "/usr/local/include/" });

    return mod;
}
