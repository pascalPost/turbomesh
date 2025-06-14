const std = @import("std");
const glfw = @import("zglfw");
const gl = @import("gl");
const State = @import("state.zig").State;

var program: gl.uint = undefined;
var vao: gl.uint = undefined;
var vbo: gl.uint = undefined;

var offset_location: gl.int = undefined;
var center_location: gl.int = undefined;
var scale_location: gl.int = undefined;
var aspect_location: gl.int = undefined;
var color_location: gl.int = undefined;

// TODO: move to state
var aspect: f32 = 1.0;

const vertices = [_]f32{
    -0.5, -0.5,
    0.5,  -0.5,
    0.0,  0.5,
};

export fn init(state: *State) callconv(.c) bool {
    // std.log.debug("lib init", .{});

    gl.makeProcTableCurrent(state.gl_proc_table_ptr);

    setCallbacks(state);

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
        \\uniform vec2 uCenter;
        \\uniform vec2 uOffset;
        \\uniform float uScale;
        \\uniform float uAspect;
        \\
        \\void main()
        \\{
        \\  vec2 scale = vec2(uScale, uScale / uAspect);
        \\  vec2 pos = (aPos - uCenter) * scale + uOffset;
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
            return false;
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
            return false;
        }
    }

    program = gl.CreateProgram();
    if (program == 0) {
        return false;
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
            return false;
        }
    }

    gl.UseProgram(program);
    gl.DeleteShader(vertex_shader);
    gl.DeleteShader(fragment_shader);

    offset_location = getUniformLocation(program, "uOffset") catch return false;
    center_location = getUniformLocation(program, "uCenter") catch return false;
    scale_location = getUniformLocation(program, "uScale") catch return false;
    aspect_location = getUniformLocation(program, "uAspect") catch return false;
    color_location = getUniformLocation(program, "uColor") catch return false;

    gl.GenVertexArrays(1, (&vao)[0..1]);
    gl.BindVertexArray(vao);

    gl.GenBuffers(1, (&vbo)[0..1]);

    gl.BindBuffer(gl.ARRAY_BUFFER, vbo);
    gl.BufferData(gl.ARRAY_BUFFER, @intCast(@sizeOf(f32) * vertices[0..].len), vertices[0..].ptr, gl.STATIC_DRAW);

    gl.VertexAttribPointer(0, 2, gl.FLOAT, gl.FALSE, 2 * @sizeOf(gl.float), 0);
    gl.EnableVertexAttribArray(0);

    // var ebo: gl.uint = undefined;
    // gl.GenBuffers(1, (&ebo)[0..1]);
    //
    // gl.BindBuffer(gl.ELEMENT_ARRAY_BUFFER, ebo);
    // gl.BufferData(gl.ELEMENT_ARRAY_BUFFER, @intCast(@sizeOf(gl.uint) * element_buffer.len), element_buffer.ptr, gl.STATIC_DRAW);

    gl.PointSize(10);

    return true;
}

export fn deinit() callconv(.c) void {
    // std.log.debug("lib deinit", .{});

    gl.BindVertexArray(0);
    gl.BindBuffer(gl.ARRAY_BUFFER, 0);
    gl.UseProgram(0);

    gl.DeleteProgram(program);
    gl.DeleteBuffers(1, (&vbo)[0..1]);
    gl.DeleteVertexArrays(1, (&vao)[0..1]);
}

export fn update(state: *State) callconv(.c) void {
    gl.makeProcTableCurrent(state.gl_proc_table_ptr);

    gl.Clear(gl.COLOR_BUFFER_BIT);

    gl.UseProgram(program);
    gl.BindVertexArray(vao);

    gl.Uniform2f(center_location, state.center[0], state.center[1]);
    gl.Uniform2f(offset_location, @floatCast(state.offset[0]), @floatCast(state.offset[1]));
    gl.Uniform1f(scale_location, state.scale);
    gl.Uniform1f(aspect_location, aspect);

    // gl.Uniform1f(state.scale_location, state.scale);
    // gl.Uniform2f(state.center_location, state.center[0], state.center[1]);
    // gl.Uniform2f(state.offset_location, @floatCast(state.offset[0]), @floatCast(state.offset[1]));
    gl.Uniform4f(color_location, 1.0, 1.0, 1.0, 1.0);

    gl.DrawArrays(gl.POINTS, 0, @intCast(vertices.len / 2));
    // gl.DrawElements(gl.LINES, @intCast(state.element_buffer.len), gl.UNSIGNED_INT, 0);
}

fn getUniformLocation(prog: gl.uint, name: [:0]const u8) !gl.int {
    const location = gl.GetUniformLocation(prog, name);
    if (location == -1) {
        return error.GetUniformLocationFailed;
    }
    return location;
}

// fn createPointBuffer(allocator: std.mem.Allocator, blocks: []discrete.Block2d) !struct {
//     points: []f32,
//     range_x: [2]f32,
//     range_y: [2]f32,
// } {
//     const num_points = blk: {
//         var total: usize = 0;
//         for (blocks) |block| {
//             total += block.points.data.len;
//         }
//         break :blk total;
//     };
//
//     const buffer = try allocator.alloc(f32, 2 * num_points);
//
//     var range_x = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };
//     var range_y = [2]f32{ std.math.floatMax(f32), std.math.floatMin(f32) };
//
//     var id: usize = 0;
//     for (blocks) |block| {
//         for (block.points.data) |point| {
//             buffer[id] = @floatCast(point.data[0]);
//             range_x[0] = @min(range_x[0], buffer[id]);
//             range_x[1] = @max(range_x[1], buffer[id]);
//
//             buffer[id + 1] = @floatCast(point.data[1]);
//             range_y[0] = @min(range_y[0], buffer[id + 1]);
//             range_y[1] = @max(range_y[1], buffer[id + 1]);
//
//             id += 2;
//         }
//     }
//
//     return .{
//         .points = buffer,
//         .range_x = range_x,
//         .range_y = range_y,
//     };
// }
//
// fn createWireframeElementBuffer(allocator: std.mem.Allocator, blocks: []const discrete.Block2d) ![]gl.uint {
//     const num_lines = blk: {
//         var total: usize = 0;
//         for (blocks) |block| {
//             const size = block.points.size;
//             const block_lines = size[0] * (size[1] - 1) * 2 + size[1] * (size[0] - 1) * 2;
//             total += block_lines;
//         }
//         break :blk total;
//     };
//
//     const buffer = try allocator.alloc(gl.uint, 8 * num_lines);
//
//     var buffer_id: usize = 0;
//     var point_offset: gl.uint = 0;
//     for (blocks) |block| {
//         const size = block.points.size;
//
//         // loop over j
//         {
//             var point_id: gl.uint = 0;
//             for (0..size[0]) |_| {
//                 for (1..size[1]) |_| {
//                     buffer[buffer_id] = point_offset + point_id;
//                     buffer[buffer_id + 1] = point_offset + point_id + 1;
//
//                     point_id += 1;
//                     buffer_id += 2;
//                 }
//
//                 point_id += 1;
//             }
//         }
//
//         // loop over i
//         {
//             for (0..size[1]) |j| {
//                 var point_id: gl.uint = @intCast(j);
//                 for (1..size[0]) |_| {
//                     buffer[buffer_id] = point_offset + point_id;
//                     buffer[buffer_id + 1] = point_offset + point_id + @as(gl.uint, @intCast(size[1]));
//
//                     point_id += @intCast(size[1]);
//                     buffer_id += 2;
//                 }
//             }
//         }
//
//         point_offset += @intCast(size[0] * size[1]);
//     }
//
//     return buffer;
// }

fn setCallbacks(state: *State) void {
    _ = glfw.setFramebufferSizeCallback(state.window, struct {
        fn func(win: *glfw.Window, width_: c_int, height_: c_int) callconv(.C) void {
            if (glfw.getWindowUserPointer(win, State)) |s| {
                s.width = @intCast(width_);
                s.height = @intCast(height_);
                s.aspect_ratio = @as(f32, @floatFromInt(width_)) / @as(f32, @floatFromInt(height_));

                // std.debug.print("width: {}, height: {}\n", .{ s.width, s.height });

                gl.Viewport(0, 0, width_, height_);
            }
        }
    }.func);

    _ = glfw.setMouseButtonCallback(state.window, struct {
        fn callback(win: *glfw.Window, button: glfw.MouseButton, action: glfw.Action, mods: glfw.Mods) callconv(.c) void {
            _ = mods;

            if (glfw.getWindowUserPointer(win, State)) |s| {
                s.dragging = button == .left and action == .press;
            }
        }
    }.callback);

    _ = glfw.setCursorPosCallback(state.window, struct {
        fn callback(win: *glfw.Window, xpos: f64, ypos: f64) callconv(.c) void {
            if (glfw.getWindowUserPointer(win, State)) |s| {
                if (s.dragging) {
                    s.offset[0] += (xpos - s.cursor_last[0]) * 2.0 / @as(f64, @floatFromInt(s.width)) * s.scroll_sensitivity;
                    s.offset[1] -= (ypos - s.cursor_last[1]) * 2.0 / @as(f64, @floatFromInt(s.height)) * s.scroll_sensitivity;
                }

                s.cursor_last = .{ xpos, ypos };
            }
        }
    }.callback);

    _ = glfw.setScrollCallback(state.window, struct {
        fn callback(win: *glfw.Window, xoffset: f64, yoffset: f64) callconv(.c) void {
            _ = xoffset;
            if (glfw.getWindowUserPointer(win, State)) |s| {
                s.scale += @floatCast(1.0 * yoffset);
            }
        }
    }.callback);
}
