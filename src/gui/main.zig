// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const glfw = @import("zglfw");
pub const gl = @import("gl");
const reload = @import("reload.zig");
const cmd = @import("cmd.zig");
const State = @import("state.zig").State;
const core = @import("core");

const templates = core.templates;

var gl_proc_table: gl.ProcTable = undefined;

const Input = struct {
    template: templates.Template,
    smoothing: struct {
        iterations: usize = 0,
        solver: core.smoothing.solver.Option,
        wall_control_function: core.smoothing.wall_control_function.Algorithm = .{ .laplace = {} },
    },
    geometry: struct {
        scale: core.types.Float = 1.0,
        pitch: core.types.Float,
        profile: core.input.ProfileInput,
    },
    output: ?[:0]const u8 = null,
    gui: ?bool = null,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    const config_file = try cmd.parseArgs();
    defer config_file.close();

    var buffer: [1024]u8 = undefined;
    var config_file_reader = config_file.reader(&buffer);

    var json_reader = std.json.Reader.init(allocator, &config_file_reader.interface);
    defer json_reader.deinit();

    const parsed = try std.json.parseFromTokenSource(Input, allocator, &json_reader, .{});
    defer parsed.deinit();

    const input = parsed.value;

    // TODO: re-parse file on modification
    // TODO: add boundary layer thickness
    // TODO: add boundary stencil comm on walls for first mesh line smoothing
    // TODO: add batch option

    // TODO: add ini, toml, yaml config files to allow comments!

    // geometry
    const profile = try core.input.create_profile(allocator, input.geometry.profile, input.geometry.scale);
    defer profile.deinit();
    const geometry = core.machine.Geometry.init(input.geometry.scale * input.geometry.pitch, profile);

    // blocking
    var mesh = try input.template.run(allocator, geometry);
    defer mesh.deinit();

    // smoothing
    try core.smoothing.smooth.mesh(allocator, &mesh, input.smoothing.iterations, input.smoothing.solver, input.smoothing.wall_control_function);

    if (input.output) |filename| {
        try mesh.write(allocator, filename);
    }

    if (input.gui == null or !input.gui.?) std.process.exit(0);

    try glfw.init();
    defer glfw.terminate();

    const gl_major = 4;
    const gl_minor = 5;
    glfw.windowHint(.context_version_major, gl_major);
    glfw.windowHint(.context_version_minor, gl_minor);
    glfw.windowHint(.opengl_profile, .opengl_core_profile);
    glfw.windowHint(.opengl_forward_compat, true);
    glfw.windowHint(.client_api, .opengl_api);
    glfw.windowHint(.doublebuffer, true);

    var width: usize = 800;
    var height: usize = 600;
    const window = try glfw.createWindow(@intCast(width), @intCast(height), "turbomesh", null);
    defer glfw.destroyWindow(window);

    glfw.makeContextCurrent(window);

    glfw.swapInterval(1);

    gl.makeProcTableCurrent(&gl_proc_table);
    if (!gl_proc_table.init(glfw.getProcAddress)) {
        return error.LoadGlAddressesFailed;
    }
    defer gl.makeProcTableCurrent(null);

    {
        var width_: c_int = undefined;
        var height_: c_int = undefined;
        glfw.getWindowSize(window, &width_, &height_);

        var xscale: f32 = undefined;
        var yscale: f32 = undefined;
        glfw.getWindowContentScale(window, &xscale, &yscale);

        const w: c_int = @intFromFloat(@as(f32, @floatFromInt(width_)) * xscale);
        const h: c_int = @intFromFloat(@as(f32, @floatFromInt(height_)) * yscale);

        width = @intCast(w);
        height = @intCast(h);

        gl.Viewport(0, 0, w, h);
    }

    var state = State{
        .allocator = allocator,
        .window = window,
        .width = width,
        .height = height,
        .gl_proc_table_ptr = &gl_proc_table,
        .zoom = 1,
        .aspect_ratio = @as(f32, @floatFromInt(width)) / @as(f32, @floatFromInt(height)),
        .center = .{ 0, 0 },
        .mesh = mesh,
    };

    glfw.setWindowUserPointer(window, &state);

    var lib = try reload.Lib.open(&state);
    defer lib.close();

    while (!window.shouldClose()) {
        try lib.reload(&state);
        glfw.pollEvents();
        lib.update(&state);
        window.swapBuffers();
    }
}
