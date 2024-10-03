const std = @import("std");

pub const Index = u32;
pub const Float = f32;
const nan = std.math.nan(Float);

pub const Index2d = struct { Index, Index };
pub const Vec2d = struct { Float, Float };

pub fn eql(a: Vec2d, b: Vec2d) bool {
    return a[0] == b[0] and a[1] == b[1];
}

/// A 2D matrix type with storage in row-major order
pub const Mat2d = struct {
    size: Index2d,
    data: []Vec2d,

    // TODO implement an iterator (?)

    pub fn init(allocator: std.mem.Allocator, size: Index2d) !Mat2d {
        const data = try allocator.alloc(Vec2d, size[0] * size[1]);
        @memset(data, .{ nan, nan });
        return Mat2d{ .size = size, .data = data };
    }

    pub fn deinit(self: Mat2d, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }

    pub fn index(self: Mat2d, i: Index2d) Index {
        return i[1] + self.size[1] * i[0];
    }
};

test "Mat2d" {
    const allocator = std.testing.allocator;
    const size = .{ 20, 10 };
    const mat = try Mat2d.init(allocator, size);
    defer mat.deinit(allocator);

    try std.testing.expect(mat.data.len == size[0] * size[1]);

    for (mat.data) |vec| {
        try std.testing.expect(std.math.isNan(vec[0]) and std.math.isNan(vec[1]));
    }
}
