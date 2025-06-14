const std = @import("std");
const discrete = @import("discrete.zig");
const glfw = @import("zglfw");
const gl = @import("gl");

pub const State = struct {
    allocator: std.mem.Allocator,

    window: *glfw.Window,
    width: usize,
    height: usize,

    gl_proc_table_ptr: *const gl.ProcTable,

    // program: gl.uint,
    // vao: gl.uint,
    //
    // offset_location: gl.int,
    // center_location: gl.int,
    // scale_location: gl.int,
    // color_location: gl.int,

    scroll_sensitivity: f32 = 1.5,

    dragging: bool = false,
    cursor_last: struct { f64, f64 } = .{ 0, 0 },

    scale: f32,
    aspect_ratio: f32,
    offset: struct { f64, f64 } = .{ 0, 0 },
    center: [2]f32,

    // point_buffer: []f32,
    // element_buffer: []gl.uint,

    // data

    mesh: ?discrete.Mesh = null,
};
