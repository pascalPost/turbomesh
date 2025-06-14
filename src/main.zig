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
const State = @import("state.zig").State;

var gl_proc_table: gl.ProcTable = undefined;

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    _ = allocator;

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
    //
    // var mesh = try template.run(allocator);
    // defer mesh.deinit();
    //
    // try smooth.mesh(allocator, &mesh, 10);
    //
    // const point_data = try createPointBuffer(allocator, mesh.blocks.items);
    // const point_buffer = point_data.points;
    // defer allocator.free(point_buffer);
    //
    // const data_width = point_data.range_x[1] - point_data.range_x[0];
    // const data_height = point_data.range_y[1] - point_data.range_y[0];
    //
    // const data_center = [2]f32{ point_data.range_x[0] + 0.5 * data_width, point_data.range_y[0] + 0.5 * data_height };
    // const data_scale: f32 = 2.0 / @max(data_width, data_height) * 0.8;
    //
    // const element_buffer = try createWireframeElementBuffer(allocator, mesh.blocks.items);
    // defer allocator.free(element_buffer);

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
        .window = window,
        .width = width,
        .height = height,
        .gl_proc_table_ptr = &gl_proc_table,
        // .program = program,
        // .vao = vao,
        // .offset_location = offset_location,
        // .center_location = center_location,
        // .scale_location = scale_location,
        // .color_location = color_location,
        .scale = 1,
        .aspect_ratio = @as(f32, @floatFromInt(width)) / @as(f32, @floatFromInt(height)),
        .center = .{ 0, 0 },
        // .point_buffer = point_buffer,
        // .element_buffer = element_buffer,
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
