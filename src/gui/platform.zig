const std = @import("std");
const wayland = @import("wayland.zig");
pub const gl = @import("gl");

var gl_proc_table: gl.ProcTable = undefined;

const Type = enum {
    wayland,
};

const Platform = union(Type) {
    wayland: void,
};

const platform = Platform{ .wayland = {} };

// const MouseState = struct {
//
// };

pub fn init(env: std.process.EnvMap, size: struct { usize, usize }) !void {
    switch (platform) {
        Type.wayland => try wayland.init(env, size, &gl_proc_table),
    }
    gl.makeProcTableCurrent(&gl_proc_table);
}

pub fn deinit() void {
    gl.makeProcTableCurrent(null);
    switch (platform) {
        Type.wayland => wayland.deinit(),
    }
}

pub fn swapBuffersAndReturnStopSignal() !bool {
    switch (platform) {
        Type.wayland => {
            try wayland.swapBuffers();
            return wayland.signal_stop;
        },
    }
}
