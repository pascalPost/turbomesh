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

    const template = o4h_template.O4H{
        .ps_csv_path = "./examples/T106/T106_ps.dat",
        .ss_csv_path = "./examples/T106/T106_ss.dat",
        .pitch = 0.08836, // m
        .blade_clustering = .{ .roberts = .{ .alpha = 0.5, .beta = 1.03 } },
        .num_cells = .{
            .o_grid = 20,
            .in_up_j = 30,
            .in_down_j = 10,
            .in_i = 10,
            .out_up_j = 40,
            .out_down_j = 10,
            .out_i = 10,
            .middle_i = 100,
            .down_j = 40,
            .bulge = 40,
            .upstream_i = 20,
            .downstream_i = 10,
        },
    };

    var mesh = try template.run(allocator);
    defer mesh.deinit();

    try smooth.mesh(allocator, &mesh, 10);

    const point_data = try createPointBuffer(allocator, mesh.blocks.items);
    const point_buffer = point_data.points;
    defer allocator.free(point_buffer);

    const data_width = point_data.range_x[1] - point_data.range_x[0];
    const data_height = point_data.range_y[1] - point_data.range_y[0];

    const data_center = [2]f32{ point_data.range_x[0] + 0.5 * data_width, point_data.range_y[0] + 0.5 * data_height };
    const data_scale: f32 = 2.0 / @max(data_width, data_height) * 0.8;

    const element_buffer = try createWireframeElementBuffer(allocator, mesh.blocks.items);
    defer allocator.free(element_buffer);

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

    // dark mode detection for each platform
    // TODO: add dark mode detection
    // TODO: add mode selection
    const dark_mode = true;
    // if (builtin.os.tag == .linux) {
    //     std.process.Child.
    // }

    if (!dark_mode) {
        gl.ClearColor(1, 1, 1, 0); // white bg
    } else {
        gl.ClearColor(0, 0, 0, 0); // black bg
    }

    const vertex_shader_source =
        \\#version 330 core
        \\layout (location = 0) in vec2 aPos;
        \\uniform float uScale;
        \\uniform vec2 uCenter;
        \\uniform vec2 uOffset; 
        \\
        \\void main()
        \\{
        \\  vec2 pos = (aPos - uCenter) * uScale + uOffset;
        \\  gl_Position = vec4(pos, 0.0, 1.0);
        \\}
    ;

    const fragment_shader_source =
        \\#version 330 core
        \\out vec4 FragColor;
        \\uniform vec4 uColor;
        \\
        \\void main()
        \\{
        \\  FragColor = uColor;
        \\}
    ;

    const vertex_shader = gl.CreateShader(gl.VERTEX_SHADER);
    gl.ShaderSource(vertex_shader, 1, &.{vertex_shader_source}, null);
    gl.CompileShader(vertex_shader);

    {
        var info_log: [512:0]u8 = undefined;
        var success: c_int = undefined;
        gl.GetShaderiv(vertex_shader, gl.COMPILE_STATUS, &success);
        if (success == gl.FALSE) {
            gl.GetShaderInfoLog(vertex_shader, info_log.len, null, &info_log);
            std.debug.print("vertex shader compilation error:\n{str}\n", .{info_log});
        }
    }

    const fragment_shader = gl.CreateShader(gl.FRAGMENT_SHADER);
    gl.ShaderSource(fragment_shader, 1, &.{fragment_shader_source}, null);
    gl.CompileShader(fragment_shader);

    {
        var info_log: [512:0]u8 = undefined;
        var success: c_int = undefined;
        gl.GetShaderiv(fragment_shader, gl.COMPILE_STATUS, &success);
        if (success == gl.FALSE) {
            gl.GetShaderInfoLog(fragment_shader, info_log.len, null, &info_log);
            std.debug.print("fragment shader compilation error:\n{str}\n", .{info_log});
        }
    }

    const program = gl.CreateProgram();
    if (program == 0) {
        return error.CreateProgramFailed;
    }

    gl.AttachShader(program, vertex_shader);
    gl.AttachShader(program, fragment_shader);
    gl.LinkProgram(program);

    {
        var success: c_int = undefined;
        gl.GetProgramiv(program, gl.LINK_STATUS, &success);
        if (success == gl.FALSE) {
            var info_log: [512:0]u8 = undefined;
            gl.GetProgramInfoLog(program, info_log.len, null, &info_log);
            std.debug.print("program linkage error:\n{str}\n", .{info_log});
        }
    }

    const offset_location = try getUniformLocation(program, "uOffset");
    const center_location = try getUniformLocation(program, "uCenter");
    const scale_location = try getUniformLocation(program, "uScale");
    const color_location = try getUniformLocation(program, "uColor");

    gl.UseProgram(program);
    gl.DeleteShader(vertex_shader);
    gl.DeleteShader(fragment_shader);

    // TODO: allow to zoom with a scaling option

    var vao: gl.uint = undefined;
    gl.GenVertexArrays(1, (&vao)[0..1]);
    gl.BindVertexArray(vao);

    var vbo: gl.uint = undefined;
    gl.GenBuffers(1, (&vbo)[0..1]);

    gl.BindBuffer(gl.ARRAY_BUFFER, vbo);
    gl.BufferData(gl.ARRAY_BUFFER, @intCast(@sizeOf(f32) * point_buffer.len), point_buffer[0..].ptr, gl.STATIC_DRAW);

    gl.VertexAttribPointer(0, 2, gl.FLOAT, gl.FALSE, 2 * @sizeOf(gl.float), 0);
    gl.EnableVertexAttribArray(0);

    var ebo: gl.uint = undefined;
    gl.GenBuffers(1, (&ebo)[0..1]);

    gl.BindBuffer(gl.ELEMENT_ARRAY_BUFFER, ebo);
    gl.BufferData(gl.ELEMENT_ARRAY_BUFFER, @intCast(@sizeOf(gl.uint) * element_buffer.len), element_buffer.ptr, gl.STATIC_DRAW);

    gl.PointSize(10);

    var state = State{
        .window = window,
        .width = width,
        .height = height,
        .gl_proc_table_ptr = &gl_proc_table,
        .program = program,
        .vao = vao,
        .offset_location = offset_location,
        .center_location = center_location,
        .scale_location = scale_location,
        .color_location = color_location,
        .scale = data_scale,
        .center = data_center,
        .point_buffer = point_buffer,
        .element_buffer = element_buffer,
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

fn getUniformLocation(program: gl.uint, name: [:0]const u8) !gl.int {
    const location = gl.GetUniformLocation(program, name);
    if (location == -1) {
        return error.GetUniformLocationFailed;
    }
    return location;
}

fn createPointBuffer(allocator: std.mem.Allocator, blocks: []discrete.Block2d) !struct {
    points: []f32,
    range_x: [2]f32,
    range_y: [2]f32,
} {
    const num_points = blk: {
        var total: usize = 0;
        for (blocks) |block| {
            total += block.points.data.len;
        }
        break :blk total;
    };

    const buffer = try allocator.alloc(f32, 2 * num_points);

    var range_x = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };
    var range_y = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };

    var id: usize = 0;
    for (blocks) |block| {
        for (block.points.data) |point| {
            buffer[id] = @floatCast(point.data[0]);
            range_x[0] = @min(range_x[0], buffer[id]);
            range_x[1] = @max(range_x[1], buffer[id]);

            buffer[id + 1] = @floatCast(point.data[1]);
            range_y[0] = @min(range_y[0], buffer[id + 1]);
            range_y[1] = @max(range_y[1], buffer[id + 1]);

            id += 2;
        }
    }

    return .{
        .points = buffer,
        .range_x = range_x,
        .range_y = range_y,
    };
}

fn createWireframeElementBuffer(allocator: std.mem.Allocator, blocks: []const discrete.Block2d) ![]gl.uint {
    const num_lines = blk: {
        var total: usize = 0;
        for (blocks) |block| {
            const size = block.points.size;
            const block_lines = size[0] * (size[1] - 1) * 2 + size[1] * (size[0] - 1) * 2;
            total += block_lines;
        }
        break :blk total;
    };

    const buffer = try allocator.alloc(gl.uint, 8 * num_lines);

    var buffer_id: usize = 0;
    var point_offset: gl.uint = 0;
    for (blocks) |block| {
        const size = block.points.size;

        // loop over j
        {
            var point_id: gl.uint = 0;
            for (0..size[0]) |_| {
                for (1..size[1]) |_| {
                    buffer[buffer_id] = point_offset + point_id;
                    buffer[buffer_id + 1] = point_offset + point_id + 1;

                    point_id += 1;
                    buffer_id += 2;
                }

                point_id += 1;
            }
        }

        // loop over i
        {
            for (0..size[1]) |j| {
                var point_id: gl.uint = @intCast(j);
                for (1..size[0]) |_| {
                    buffer[buffer_id] = point_offset + point_id;
                    buffer[buffer_id + 1] = point_offset + point_id + @as(gl.uint, @intCast(size[1]));

                    point_id += @intCast(size[1]);
                    buffer_id += 2;
                }
            }
        }

        point_offset += @intCast(size[0] * size[1]);
    }

    return buffer;
}
