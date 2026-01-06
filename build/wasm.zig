const std = @import("std");
const core = @import("core.zig");

pub fn addWasm(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    cfg: core.Config,
) void {
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });
    const wasm = b.addExecutable(.{
        .name = "turbomesh",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/wasm/lib.zig"),
            .target = wasm_target,
            .optimize = optimize,
        }),
    });
    wasm.rdynamic = true;
    wasm.entry = .disabled;
    b.installArtifact(wasm);

    const core_mod = core.addCoreModule(b, target, optimize, cfg);
    wasm.root_module.addImport("core", core_mod);
}
