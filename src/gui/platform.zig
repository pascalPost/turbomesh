const std = @import("std");
const wayland = @import("wayland.zig");

const Type = enum {
    wayland,
};

const Platform = union(Type) {
    wayland: void,
};

const platform = Platform{ .wayland = {} };

pub fn init(env: std.process.EnvMap, size: struct { usize, usize }) !void {
    switch (platform) {
        Type.wayland => try wayland.init(env, size),
    }
}

pub fn deinit() void {
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
