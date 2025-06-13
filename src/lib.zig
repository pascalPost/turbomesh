const std = @import("std");
const glfw = @import("zglfw");
const gl = @import("gl");
const State = @import("state.zig").State;

export fn init(state: *State) callconv(.c) void {
    std.log.debug("lib init\n", .{});

    _ = glfw.setFramebufferSizeCallback(state.window, struct {
        fn func(win: *glfw.Window, width_: c_int, height_: c_int) callconv(.C) void {
            if (glfw.getWindowUserPointer(win, State)) |s| {
                s.width = @intCast(width_);
                s.height = @intCast(height_);

                std.debug.print("width: {}, height: {}\n", .{ s.width, s.height });

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

export fn update(state: *State) callconv(.c) void {
    gl.makeProcTableCurrent(state.gl_proc_table_ptr);

    gl.Clear(gl.COLOR_BUFFER_BIT);

    gl.UseProgram(state.program);
    gl.Uniform1f(state.scale_location, state.scale);
    gl.Uniform2f(state.center_location, state.center[0], state.center[1]);
    gl.Uniform2f(state.offset_location, @floatCast(state.offset[0]), @floatCast(state.offset[1]));
    gl.BindVertexArray(state.vao);
    gl.Uniform4f(state.color_location, 1.0, 1.0, 1.0, 1.0);
    // gl.DrawArrays(gl.POINTS, 0, @intCast(state.point_buffer.len / 2));
    gl.DrawElements(gl.LINES, @intCast(state.element_buffer.len), gl.UNSIGNED_INT, 0);
}
