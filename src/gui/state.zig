// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const core = @import("core");
const glfw = @import("zglfw");
const gl = @import("gl");

const discrete = core.discrete;

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

    mouse_sensitivity: f32 = 1.5,
    dragging: bool = false,
    cursor_last: struct { f64, f64 } = .{ 0, 0 },

    aspect_ratio: f32,
    offset: struct { f64, f64 } = .{ 0, 0 },
    center: [2]f32,

    zoom_speed: f32 = 0.1,
    zoom: f32,

    // data

    mesh: ?discrete.Mesh = null,
};
