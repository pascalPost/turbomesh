const std = @import("std");
const o4h = @import("O4H.zig");
const discrete = @import("../discrete.zig");

const Type = enum {
    O4H,
};

pub const Template = union(Type) {
    O4H: o4h.O4H,

    pub fn run(self: Template, allocator: std.mem.Allocator) !discrete.Mesh {
        try switch (self) {
            .O4H => |t| return t.run(allocator),
        };
    }
};
