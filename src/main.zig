const std = @import("std");
const builtin = @import("builtin");
// const discrete = @import("discrete.zig");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");
const glfw = @import("zglfw");
pub const gl = @import("gl");

var gl_proc_table: gl.ProcTable = undefined;

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

var dragging = false;
var cursor_last: struct { f64, f64 } = .{ 0, 0 };
var offset: struct { f64, f64 } = .{ 0, 0 };
const sensitivity = 1.5;

var width: usize = 800;
var height: usize = 600;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    const block_points = [_]Mat2d{
        try Mat2d.init(allocator, .{ 21, 5 }),
    };

    defer {
        for (block_points) |block| {
            block.deinit(allocator);
        }
    }

    for (block_points) |block| {
        const size = block.size;

        var idx: usize = 0;
        var i: usize = 0;
        while (i < size[0]) : (i += 1) {
            var j: usize = 0;
            while (j < size[1]) : (j += 1) {
                block.data[idx] = Vec2d.init(@floatFromInt(i), @floatFromInt(j));
                idx += 1;
            }
        }
    }

    const point_data = try createPointBuffer(allocator, &block_points);
    const point_buffer = point_data.points;
    defer allocator.free(point_buffer);

    const data_width = point_data.range_x[1] - point_data.range_x[0];
    const data_height = point_data.range_y[1] - point_data.range_y[0];

    const data_center = [2]f32{ point_data.range_x[0] + 0.5 * data_width, point_data.range_y[0] + 0.5 * data_height };
    const data_scale: f32 = 2.0 / @max(data_width, data_height) * 0.8;

    const element_buffer = try createWireframeElementBuffer(allocator, &block_points);
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

    const window = try glfw.createWindow(@intCast(width), @intCast(height), "zig-gamedev: minimal_glfw_gl", null);
    defer glfw.destroyWindow(window);

    glfw.makeContextCurrent(window);

    glfw.swapInterval(1);

    _ = glfw.setMouseButtonCallback(window, struct {
        fn callback(win: *glfw.Window, button: glfw.MouseButton, action: glfw.Action, mods: glfw.Mods) callconv(.c) void {
            _ = win;
            _ = mods;
            dragging = button == .left and action == .press;
        }
    }.callback);

    _ = glfw.setCursorPosCallback(window, struct {
        fn callback(win: *glfw.Window, xpos: f64, ypos: f64) callconv(.c) void {
            _ = win;

            if (dragging) {
                offset[0] += (xpos - cursor_last[0]) * 2.0 / @as(f64, @floatFromInt(width)) * sensitivity;
                offset[1] -= (ypos - cursor_last[1]) * 2.0 / @as(f64, @floatFromInt(height)) * sensitivity;
            }

            cursor_last = .{ xpos, ypos };
        }
    }.callback);

    gl.makeProcTableCurrent(&gl_proc_table);
    if (!gl_proc_table.init(glfw.getProcAddress)) {
        return error.LoadGlAddressesFailed;
    }
    defer gl.makeProcTableCurrent(null);

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

    _ = glfw.setFramebufferSizeCallback(window, struct {
        fn func(win: *glfw.Window, width_: c_int, height_: c_int) callconv(.C) void {
            _ = win;
            width = @intCast(width_);
            height = @intCast(height_);
            gl.Viewport(0, 0, width_, height_);
        }
    }.func);

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

    while (!window.shouldClose()) {
        glfw.pollEvents();

        gl.Clear(gl.COLOR_BUFFER_BIT);

        gl.UseProgram(program);
        gl.Uniform1f(scale_location, data_scale);
        gl.Uniform2f(center_location, data_center[0], data_center[1]);
        gl.Uniform2f(offset_location, @floatCast(offset[0]), @floatCast(offset[1]));
        gl.BindVertexArray(vao);
        gl.Uniform4f(color_location, 1.0, 1.0, 1.0, 1.0);
        gl.DrawArrays(gl.POINTS, 0, @intCast(point_buffer.len / 2));
        gl.DrawElements(gl.LINES, @intCast(element_buffer.len), gl.UNSIGNED_INT, 0);

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

fn createPointBuffer(allocator: std.mem.Allocator, points: []const Mat2d) !struct {
    points: []f32,
    range_x: [2]f32,
    range_y: [2]f32,
} {
    const num_points = blk: {
        var total: usize = 0;
        for (points) |block| {
            total += block.data.len;
        }
        break :blk total;
    };

    const buffer = try allocator.alloc(f32, 2 * num_points);

    var range_x = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };
    var range_y = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };

    var id: usize = 0;
    for (points) |block| {
        for (block.data) |point| {
            buffer[id] = @floatCast(point.data[0]);
            range_x[0] = @min(range_x[0], buffer[id]);
            range_x[1] = @max(range_x[1], buffer[id]);

            buffer[id + 1] = @floatCast(point.data[1]);
            range_y[0] = @min(range_y[0], buffer[id]);
            range_y[1] = @max(range_y[1], buffer[id]);

            id += 2;
        }
    }

    return .{
        .points = buffer,
        .range_x = range_x,
        .range_y = range_y,
    };
}

fn createWireframeElementBuffer(allocator: std.mem.Allocator, points: []const Mat2d) ![]gl.uint {
    const num_lines = blk: {
        var total: usize = 0;
        for (points) |block| {
            const block_lines = block.size[0] * (block.size[1] - 1) * 2 + block.size[1] * (block.size[0] - 1) * 2;
            total += block_lines;
        }
        break :blk total;
    };

    const buffer = try allocator.alloc(gl.uint, 8 * num_lines);

    var buffer_id: usize = 0;
    for (points) |block| {
        const size = block.size;

        // loop over j
        {
            var point_id: gl.uint = 0;
            for (0..size[0]) |_| {
                for (1..size[1]) |_| {
                    buffer[buffer_id] = point_id;
                    buffer[buffer_id + 1] = point_id + 1;

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
                    buffer[buffer_id] = point_id;
                    buffer[buffer_id + 1] = point_id + @as(gl.uint, @intCast(size[1]));

                    point_id += @intCast(size[1]);
                    buffer_id += 2;
                }
            }
        }
    }

    return buffer;
}
