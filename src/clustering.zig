const std = @import("std");
const types = @import("types.zig");
const geometry = @import("geometry.zig");

const Index = types.Index;
const Float = types.Float;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Mat2d = types.Mat2d;

fn uniform(points: Index) void {
    var i: Index = 0;
    while (i < points) : (i += 1) {
        const u: Float = i / (points - 1);
    }
}

const Edge2d = struct {
    points: []Vec2d,

    fn init(allocator: std.mem.Allocator, size: Index2d) !Edge2d {
        const points = try Mat2d.init(allocatior, size);
        return Edge2d{ points };
    }
};

test "uniform clustering between 0 and 1" {
    const allocator = std.testing.allocator;

    // given a curve
    const line = geometry.Line2d{ .start = .{ 0.0, 0.0 }, .end = .{ 2.0, 0.0 } };

    // and an edge (that will be the discrete reprensentation of the curve)
    const edge = try Edge2d.init(allocator, .{})
}
