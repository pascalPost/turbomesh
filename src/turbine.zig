const std = @import("std");
const types = @import("types.zig");

const Float = types.Float;
const Vec2d = types.Vec2d;

const TurbineTemplate = struct {
    ps_csv_path: []const u8,
    ss_csv_path: []const u8,

    fn run() void {}
};

// test "execute turbine template" {
//     const allocator = std.testing.allocator;
// }
