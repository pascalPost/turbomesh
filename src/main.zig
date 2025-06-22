const std = @import("std");
const builtin = @import("builtin");
const discrete = @import("discrete.zig");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");
const glfw = @import("zglfw");
pub const gl = @import("gl");
const o4h_template = @import("templates/O4H.zig");
const smooth = @import("smooth.zig");
const reload = @import("reload.zig");
const cmd = @import("cmd.zig");
const State = @import("state.zig").State;

var gl_proc_table: gl.ProcTable = undefined;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    const config_file = try cmd.parseArgs();
    defer config_file.close();

    var json_reader = std.json.reader(allocator, config_file.reader());
    defer json_reader.deinit();

    const parsed = try std.json.parseFromTokenSource(o4h_template.O4H, allocator, &json_reader, .{});
    defer parsed.deinit();

    // TODO: re-parse file on modification
    // TODO: add boundary layer thickness
    // TODO: add boundary stencil comm on walls for first mesh line smoothing
    // TODO: add batch option

    const template = parsed.value;

    // const template = o4h_template.O4H{
    //     .ps_csv_path = "./examples/T106/T106_ps.dat",
    //     .ss_csv_path = "./examples/T106/T106_ss.dat",
    //     .pitch = 0.08836, // m
    //     .blade_clustering = .{ .roberts = .{ .alpha = 0.5, .beta = 1.03 } },
    //     .num_cells = .{
    //         .o_grid = 20,
    //         .in_up_j = 30,
    //         .in_down_j = 10,
    //         .in_i = 10,
    //         .out_up_j = 40,
    //         .out_down_j = 10,
    //         .out_i = 10,
    //         .middle_i = 100,
    //         .down_j = 40,
    //         .bulge = 40,
    //         .upstream_i = 20,
    //         .downstream_i = 10,
    //     },
    // };

    // TODO: use this to write config file
    // const file = try std.fs.cwd().createFile("T106.json", .{});
    // defer file.close();
    //
    // try std.json.stringify(template, .{ .whitespace = .indent_4 }, file.writer());

    var mesh = try template.run(allocator);
    defer mesh.deinit();

    try smooth.mesh(allocator, &mesh, 20, .petsc);

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
