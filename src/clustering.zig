const std = @import("std");
const types = @import("types.zig");

const Float = types.Float;

const Uniform = struct {
    fn compute(self: Uniform, data: []Float) void {
        _ = self;
        const n = data.len;
        for (0..n) |i| {
            data[i] = @as(Float, @floatFromInt(i)) / @as(Float, @floatFromInt(n - 1));
        }
    }
};

/// Roberts cluster function, see
/// https://github.com/luohancfd/CFCFD-NG/blob/dev/lib/nm/source/fobject.cxx
/// alpha = 0.5 cluster at both ends
/// alpha = 0.0 cluster toward t=1.0
/// stretching factor 1.0 < beta < +inf, closer to 1.0 gives stronger clustering.
pub const Roberts = struct {
    alpha: Float,
    beta: Float,

    pub fn compute(self: Roberts, data: []Float) void {
        const n = data.len;
        const alpha = self.alpha;
        const beta = self.beta;

        std.debug.assert(n > 1);

        for (0..n) |i| {
            const u = @as(Float, @floatFromInt(i)) / @as(Float, @floatFromInt(n - 1));
            const tmp = std.math.pow(Float, (beta + 1.0) / (beta - 1.0), (u - alpha) / (1.0 - alpha));
            const tbar = (beta + 2.0 * alpha) * tmp - beta + 2.0 * alpha;
            data[i] = tbar / ((2.0 * alpha + 1.0) * (1.0 + tmp));
        }
    }
};

pub const Tag = enum {
    uniform,
    roberts,
};

pub const Function = union(Tag) {
    // NOTE: a pointer based implementation optimizes space.
    uniform: Uniform,
    roberts: Roberts,
};

pub fn create(allocator: std.mem.Allocator, clustering: Function, n: types.Index) ![]types.Float {
    var u = try allocator.alloc(types.Float, n);
    switch (clustering) {
        inline else => |f| f.compute(u[0..]),
    }
    return u;
}
