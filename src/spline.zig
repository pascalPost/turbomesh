// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

/// Natural cubic spline that parameterizes the curve by (approximate) arc length
/// over the range [0, 1]. The input points are re-parameterized by cumulative
/// chord length and a lookup table is built to map normalized arc length back
/// into the spline parameter domain.
pub fn FittingSpline(comptime dim: usize) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        params: []f64, // normalized chord-length parameters for control points
        points: []const [dim]f64,
        second_derivs: [dim][]f64,
        sample_params: []f64, // parameter samples used for arc-length lookup
        sample_arc: []f64, // cumulative arc length over samples, normalized to [0, 1]
        total_length: f64,

        const sample_count = 200;

        pub fn init(allocator: std.mem.Allocator, points: []const [dim]f64, degree: usize) !Self {
            if (degree != 3) return error.UnsupportedDegree;
            if (points.len < 2) return error.NotEnoughPoints;

            var params = try allocator.alloc(f64, points.len);
            errdefer allocator.free(params);

            const points_copy = try allocator.alloc([dim]f64, points.len);
            errdefer allocator.free(points_copy);
            @memcpy(points_copy, points);

            const total_chord = computeChordParams(points, params[0..]);

            // build natural spline second derivatives for each dimension
            var second_derivs: [dim][]f64 = undefined;
            inline for (0..dim) |d| {
                second_derivs[d] = try computeSecondDerivs(allocator, params, points_copy, d);
            }
            errdefer {
                inline for (0..dim) |d| allocator.free(second_derivs[d]);
            }

            const sample_params = try allocator.alloc(f64, sample_count + 1);
            errdefer allocator.free(sample_params);
            const sample_arc = try allocator.alloc(f64, sample_count + 1);
            errdefer allocator.free(sample_arc);

            var spline = Self{
                .allocator = allocator,
                .params = params,
                .points = points_copy,
                .second_derivs = second_derivs,
                .sample_params = sample_params,
                .sample_arc = sample_arc,
                .total_length = total_chord,
            };

            try spline.buildArcLengthTable();

            return spline;
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.params);
            self.allocator.free(self.points);
            inline for (0..dim) |d| self.allocator.free(self.second_derivs[d]);
            self.allocator.free(self.sample_params);
            self.allocator.free(self.sample_arc);
        }

        pub fn interpolate(self: *const Self, u: []const f64, values: [][dim]f64) !void {
            if (u.len != values.len) return error.Mismatch;

            for (u, values) |u_val, *out| {
                const param = self.paramAtArcFraction(u_val);
                out.* = self.eval(param);
            }
        }

        pub fn integrate(self: *const Self) f64 {
            return self.total_length;
        }

        fn buildArcLengthTable(self: *Self) !void {
            // sample parameters evenly in the spline domain
            for (self.sample_params, 0..) |*p, i| {
                p.* = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(self.sample_params.len - 1));
            }

            var length: f64 = 0.0;
            self.sample_arc[0] = 0.0;
            var prev = self.eval(self.sample_params[0]);
            for (self.sample_params[1..], self.sample_arc[1..]) |p, *arc| {
                const cur = self.eval(p);
                length += distance(prev, cur);
                arc.* = length;
                prev = cur;
            }

            self.total_length = length;
            if (length == 0.0) {
                @memset(self.sample_arc, 0.0);
                return;
            }

            for (self.sample_arc) |*v| v.* /= length;
        }

        fn paramAtArcFraction(self: *const Self, u: f64) f64 {
            if (self.total_length == 0.0 or self.sample_arc.len == 0) return 0.0;

            const target = std.math.clamp(u, 0.0, 1.0);

            // upper-bound binary search on normalized arc table
            var lo: usize = 0;
            var hi: usize = self.sample_arc.len - 1;
            while (lo < hi) {
                const mid = (lo + hi) / 2;
                if (self.sample_arc[mid] < target) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }

            if (lo == 0) return self.sample_params[0];
            if (lo >= self.sample_arc.len) return self.sample_params[self.sample_params.len - 1];

            const a0 = self.sample_arc[lo - 1];
            const a1 = self.sample_arc[lo];
            const p0 = self.sample_params[lo - 1];
            const p1 = self.sample_params[lo];

            const t = if (a1 > a0) (target - a0) / (a1 - a0) else 0.0;
            return p0 + t * (p1 - p0);
        }

        fn computeChordParams(points: []const [dim]f64, params: []f64) f64 {
            params[0] = 0.0;
            var total: f64 = 0.0;
            for (points[1..], 1..) |p, i| {
                total += distance(points[i - 1], p);
                params[i] = total;
            }
            if (total == 0.0) {
                const denom = @as(f64, @floatFromInt(params.len - 1));
                for (params, 0..) |*p, i| p.* = @as(f64, @floatFromInt(i)) / denom;
                return 0.0;
            }
            for (params) |*p| p.* /= total;
            return total;
        }

        fn computeSecondDerivs(allocator: std.mem.Allocator, params: []const f64, points: []const [dim]f64, component: usize) ![]f64 {
            const n = params.len;
            var z = try allocator.alloc(f64, n);
            errdefer allocator.free(z);

            if (n == 2) {
                @memset(z, 0.0);
                return z;
            }

            var tmp = try allocator.alloc(f64, n);
            defer allocator.free(tmp);

            z[0] = 0.0;
            tmp[0] = 0.0;

            var i: usize = 1;
            while (i < n - 1) : (i += 1) {
                const h_im1 = params[i] - params[i - 1];
                const h_i = params[i + 1] - params[i];
                if (h_im1 == 0.0 or h_i == 0.0) return error.CoincidentParameters;

                const dy_im1 = points[i][component] - points[i - 1][component];
                const dy_i = points[i + 1][component] - points[i][component];

                const alpha = (dy_i / h_i) - (dy_im1 / h_im1);
                const denom = 2.0 * (params[i + 1] - params[i - 1]) - h_im1 * tmp[i - 1];
                const mu = h_i / denom;

                tmp[i] = mu;
                z[i] = (6.0 * alpha - h_im1 * z[i - 1]) / denom;
            }

            z[n - 1] = 0.0;

            var k = n - 2;
            while (true) {
                z[k] = z[k] - tmp[k] * z[k + 1];
                if (k == 0) break;
                k -= 1;
            }

            return z;
        }

        fn eval(self: *const Self, param: f64) [dim]f64 {
            const u = std.math.clamp(param, 0.0, 1.0);
            var idx: usize = 0;
            while (idx + 1 < self.params.len and self.params[idx + 1] < u) idx += 1;
            if (idx >= self.params.len - 1) idx = self.params.len - 2;

            const h = self.params[idx + 1] - self.params[idx];
            std.debug.assert(h > 0.0);
            const a = (self.params[idx + 1] - u) / h;
            const b = (u - self.params[idx]) / h;

            var result: [dim]f64 = undefined;
            inline for (0..dim) |d| {
                const y0 = self.points[idx][d];
                const y1 = self.points[idx + 1][d];
                const z0 = self.second_derivs[d][idx];
                const z1 = self.second_derivs[d][idx + 1];
                result[d] = a * y0 + b * y1 + ((a * a * a - a) * z0 + (b * b * b - b) * z1) * (h * h) / 6.0;
            }
            return result;
        }

        fn distance(a: [dim]f64, b: [dim]f64) f64 {
            var sum: f64 = 0.0;
            inline for (0..dim) |i| {
                const d = b[i] - a[i];
                sum += d * d;
            }
            return std.math.sqrt(sum);
        }
    };
}

test "spline interpolating a straight line" {
    const allocator = std.testing.allocator;
    const dim = 2;
    const degree = 3;
    const pts = [_][dim]f64{
        .{ 0.0, 0.0 },
        .{ 0.5, 0.5 },
        .{ 1.0, 1.0 },
        .{ 2.0, 2.0 },
        .{ 3.0, 3.0 },
        .{ 4.0, 4.0 },
    };

    var spline = try FittingSpline(dim).init(allocator, pts[0..], degree);
    defer spline.deinit();

    var values: [6][dim]f64 = undefined;
    const u = [_]f64{ 0.0, 0.125, 0.25, 0.5, 0.75, 1.0 };
    try spline.interpolate(&u, &values);

    const expected = pts;

    for (values, expected) |v, e| {
        try std.testing.expectApproxEqAbs(e[0], v[0], 1e-9);
        try std.testing.expectApproxEqAbs(e[1], v[1], 1e-9);
    }

    try std.testing.expectApproxEqAbs(std.math.sqrt(2.0) * 4.0, spline.integrate(), 1e-9);
}

test "monotonic arc-length mapping on curve" {
    const allocator = std.testing.allocator;
    const dim = 2;
    const degree = 3;
    const pts = [_][dim]f64{
        .{ 0.0, 0.0 },
        .{ 1.0, 0.5 },
        .{ 2.0, 1.5 },
        .{ 2.5, 3.0 },
    };

    var spline = try FittingSpline(dim).init(allocator, pts[0..], degree);
    defer spline.deinit();

    var values: [3][dim]f64 = undefined;
    const u = [_]f64{ 0.0, 0.5, 1.0 };
    try spline.interpolate(&u, &values);

    try std.testing.expect(values[0][0] <= values[1][0]);
    try std.testing.expect(values[1][0] <= values[2][0]);
    try std.testing.expectApproxEqAbs(pts[0][0], values[0][0], 1e-9);
    try std.testing.expectApproxEqAbs(pts[0][1], values[0][1], 1e-9);
    try std.testing.expectApproxEqAbs(pts[pts.len - 1][0], values[2][0], 1e-9);
    try std.testing.expectApproxEqAbs(pts[pts.len - 1][1], values[2][1], 1e-9);
}

test "length of two-point spline" {
    const allocator = std.testing.allocator;
    const dim = 2;
    const degree = 3;
    const pts = [_][dim]f64{
        .{ 0.0, 0.0 },
        .{ 0.0, 3.0 },
    };

    var spline = try FittingSpline(dim).init(allocator, pts[0..], degree);
    defer spline.deinit();

    try std.testing.expectApproxEqAbs(3.0, spline.integrate(), 1e-9);
}
