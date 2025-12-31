// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

const core = @import("build/core.zig");
const gui = @import("build/gui.zig");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const core_mod = core.addCoreModule(b, target, optimize);

    gui.addDesktopGui(b, target, optimize, core_mod);
}
