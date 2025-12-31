// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

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

/// SingleHyperbolicTangentClustering implements the hyperbolic tangent clustering function
/// matching the specified spacing of the first cell approximately.
///
/// Details can be find in:
///
/// Vinokur, Marcel. “On One-Dimensional Stretching Functions for Finite-Difference Calculations.”
/// Journal of Computational Physics 50, no. 2 (May 1983): 215–34.
/// https://doi.org/10.1016/0021-9991(83)90065-7.
///
/// Thompson, Joe F. “A General Three-Dimensional Elliptic Grid Generation System on a Composite Block Structure.”
/// Computer Methods in Applied Mechanics and Engineering 64, no. 1–3 (October 1987): 377–411.
/// https://doi.org/10.1016/0045-7825(87)90047-8.
pub const SingleHyperbolicClustering = struct {
    delta_s: Float,

    pub fn compute(self: SingleHyperbolicClustering, data: []Float) void {
        // TODO: enhance implementation: it seems that only the derivative at the wall is specified not
        // the spacing itself, see https://www.cfd-online.com/Wiki/Structured_mesh_generation
        // or in https://www.osti.gov/etdeweb/servlets/purl/632740 page 9:
        // ds/deta is approximately equal to the normalized height of the first cell delta_1 / delta_s
        // For this, the more general implementation of the Vinokur clustering might be needed, where
        // the spacing can be specified at any point

        // compute B
        const n_1: Float = @floatFromInt(data.len - 1);
        const b = n_1 * self.delta_s;

        const y = 1.0 / b;

        // eq. 63 to 67 in Vinokur 1983
        const delta = if (y < 2.7829681) blk: {
            const y_bar = y - 1.0;
            break :blk @sqrt(6.0 * y_bar) * (1.0 + y_bar * (-0.15 + y_bar * (0.057321429 + y_bar * (-0.024907295 + y_bar * (0.0077424461 - 0.0010794123 * y_bar)))));
        } else blk: {
            const w = 1.0 / y - 0.028527431;
            const v = @log(y);
            break :blk v + (1.0 + 1.0 / v) * @log(2.0 * v) - 0.02041793 + w * (0.24902722 + w * (1.9496443 + w * (-2.6294547 + 8.56795911 * w)));
        };

        for (data, 0..) |*xi, i| {
            xi.* = @as(Float, @floatFromInt(i)) / (n_1);
        }

        for (data[1..]) |*xi| {
            const s = 1.0 + std.math.tanh(0.5 * delta * (xi.* - 1.0)) / std.math.tanh(0.5 * delta);
            xi.* = s;
        }

        std.debug.assert(data[0] == 0.0);
        std.debug.assert(data[data.len - 1] == 1.0);
    }
};

pub const Tag = enum {
    uniform,
    roberts,
    single_hyperbolic_clustering,
};

pub const Function = union(Tag) {
    // NOTE: a pointer based implementation optimizes space.
    uniform: Uniform,
    roberts: Roberts,
    single_hyperbolic_clustering: SingleHyperbolicClustering,
};

pub fn create(allocator: std.mem.Allocator, clustering: Function, n: types.Index) ![]types.Float {
    var u = try allocator.alloc(types.Float, n);
    switch (clustering) {
        inline else => |f| f.compute(u[0..]),
    }
    return u;
}
