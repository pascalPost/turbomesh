const std = @import("std");

pub const Index = usize;
pub const Float = f64;
const nan = std.math.nan(Float);

/// 2D point index type.
pub const Index2d = struct { Index, Index };

/// Mesh global 2D block point index type.
pub const MeshIndex2d = struct { block: usize, point: Index2d };

pub const Vec2d = struct {
    data: [2]Float,

    pub fn init(v0: Float, v1: Float) Vec2d {
        return .{ .data = .{ v0, v1 } };
    }

    pub fn add(self: *Vec2d, v: Vec2d) void {
        self.data[0] += v.data[0];
        self.data[1] += v.data[1];
    }
};

pub fn castConst(data: []const Vec2d) []const [2]Float {
    const bytes = std.mem.sliceAsBytes(data[0..]);
    return std.mem.bytesAsSlice([2]Float, bytes);
}

pub fn cast(data: []Vec2d) [][2]Float {
    const bytes = std.mem.sliceAsBytes(data[0..]);
    return std.mem.bytesAsSlice([2]Float, bytes);
}

pub fn eql(a: Vec2d, b: Vec2d) bool {
    return std.mem.eql(Float, a.data[0..], b.data[0..]);
}

pub fn eqlApprox(a: Vec2d, b: Vec2d, tol: Float) bool {
    return @abs(a.data[0] - b.data[0]) <= tol and @abs(a.data[1] - b.data[1]) <= tol;
}

pub fn add(a: Vec2d, b: Vec2d) Vec2d {
    return .{ .data = .{ a.data[0] + b.data[0], a.data[1] + b.data[1] } };
}

pub fn addAll(vec: anytype) Vec2d {
    var res = Vec2d.init(0, 0);
    inline for (vec) |v| res = add(res, v);
    return res;
}

pub fn sub(a: Vec2d, b: Vec2d) Vec2d {
    return .{ .data = .{ a.data[0] - b.data[0], a.data[1] - b.data[1] } };
}

pub fn scale(s: Float, v: Vec2d) Vec2d {
    return .{ .data = .{ s * v.data[0], s * v.data[1] } };
}

pub fn negate(v: Vec2d) Vec2d {
    return .{ .data = .{ -v.data[0], -v.data[1] } };
}

pub fn abs(v: Vec2d) Float {
    return std.math.sqrt(v.data[0] * v.data[0] + v.data[1] * v.data[1]);
}

pub fn neg(v: Vec2d) Vec2d {
    return .{ .data = .{ -v.data[0], -v.data[1] } };
}

/// A 2D matrix type with storage in column-major order
pub const Mat2d = struct {
    size: Index2d,
    data: []Vec2d,

    // TODO implement an iterator (?)

    pub fn init(allocator: std.mem.Allocator, size: Index2d) !Mat2d {
        var data = try allocator.alloc(Vec2d, size[0] * size[1]);
        for (data[0..]) |*d| d.data = [2]Float{ nan, nan };
        return Mat2d{ .size = size, .data = data };
    }

    pub fn deinit(self: Mat2d, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }

    pub fn index(self: Mat2d, i: Index2d) Index {
        return i[1] + self.size[1] * i[0];
    }

    pub fn getIndex(self: Mat2d, i: Index2d) Vec2d {
        return self.data[self.index(i)];
    }
};

test "Mat2d" {
    const allocator = std.testing.allocator;
    const size = .{ 20, 10 };
    const mat = try Mat2d.init(allocator, size);
    defer mat.deinit(allocator);

    try std.testing.expect(mat.data.len == size[0] * size[1]);

    for (mat.data) |vec| {
        try std.testing.expect(std.math.isNan(vec.data[0]) and std.math.isNan(vec.data[1]));
    }
}
