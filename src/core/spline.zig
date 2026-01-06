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

        pub fn deinit(self: *const Self) void {
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

test "T106 blade coordinate integration" {
    const allocator = std.testing.allocator;
    const dim = 2;
    const degree = 3;

    // T106 blade coordinates over chord (Table I-2, page 213, "The Effects of Wakes on Separating Boundary Layers in Low Pressure Turbines", Rory Douglas Stieger,
    // https://www-g.eng.cam.ac.uk/whittle/T106/T106A_cascade/T106A_cascade_publications/RStieger_PhD.pdf)
    var pts = [_][dim]f64{
        .{ 0.852022, -0.499894 },
        .{ 0.854287, -0.499152 },
        .{ 0.856024, -0.497879 },
        .{ 0.857388, -0.496202 },
        .{ 0.858425, -0.494252 },
        .{ 0.858921, -0.492134 },
        .{ 0.858849, -0.489945 },
        .{ 0.858336, -0.487832 },
        .{ 0.856981, -0.484433 },
        .{ 0.854884, -0.479210 },
        .{ 0.852085, -0.472232 },
        .{ 0.848607, -0.463560 },
        .{ 0.844469, -0.453274 },
        .{ 0.839713, -0.441446 },
        .{ 0.834337, -0.428134 },
        .{ 0.828372, -0.413425 },
        .{ 0.821836, -0.397388 },
        .{ 0.814726, -0.380103 },
        .{ 0.807069, -0.361644 },
        .{ 0.798872, -0.342097 },
        .{ 0.790120, -0.321545 },
        .{ 0.780839, -0.300063 },
        .{ 0.771018, -0.277748 },
        .{ 0.760651, -0.254690 },
        .{ 0.749746, -0.230987 },
        .{ 0.738285, -0.206733 },
        .{ 0.726255, -0.182047 },
        .{ 0.713682, -0.156981 },
        .{ 0.700516, -0.131702 },
        .{ 0.686776, -0.106292 },
        .{ 0.672425, -0.080891 },
        .{ 0.657461, -0.055603 },
        .{ 0.641855, -0.030573 },
        .{ 0.625580, -0.005925 },
        .{ 0.617178, 0.006200 },
        .{ 0.608581, 0.018188 },
        .{ 0.599805, 0.029981 },
        .{ 0.590809, 0.041609 },
        .{ 0.581625, 0.052984 },
        .{ 0.572194, 0.064157 },
        .{ 0.562565, 0.075027 },
        .{ 0.552660, 0.085647 },
        .{ 0.542536, 0.095906 },
        .{ 0.532108, 0.105857 },
        .{ 0.521465, 0.115384 },
        .{ 0.510501, 0.124541 },
        .{ 0.499342, 0.133204 },
        .{ 0.487867, 0.141440 },
        .{ 0.476249, 0.149133 },
        .{ 0.464341, 0.156368 },
        .{ 0.452322, 0.163044 },
        .{ 0.440047, 0.169240 },
        .{ 0.427639, 0.174893 },
        .{ 0.415012, 0.180039 },
        .{ 0.402272, 0.184630 },
        .{ 0.389355, 0.188696 },
        .{ 0.376273, 0.192225 },
        .{ 0.363052, 0.195189 },
        .{ 0.349624, 0.197578 },
        .{ 0.336102, 0.199358 },
        .{ 0.322457, 0.200535 },
        .{ 0.308772, 0.201087 },
        .{ 0.294991, 0.201010 },
        .{ 0.281229, 0.200303 },
        .{ 0.267419, 0.198967 },
        .{ 0.253684, 0.197008 },
        .{ 0.239943, 0.194412 },
        .{ 0.226333, 0.191202 },
        .{ 0.212761, 0.187358 },
        .{ 0.199372, 0.182919 },
        .{ 0.186068, 0.177858 },
        .{ 0.172994, 0.172231 },
        .{ 0.160048, 0.166003 },
        .{ 0.147370, 0.159246 },
        .{ 0.134837, 0.151909 },
        .{ 0.122597, 0.144092 },
        .{ 0.110472, 0.135697 },
        .{ 0.098656, 0.126874 },
        .{ 0.086915, 0.117442 },
        .{ 0.075496, 0.107622 },
        .{ 0.064090, 0.097111 },
        .{ 0.053025, 0.086239 },
        .{ 0.041933, 0.074534 },
        .{ 0.031245, 0.062457 },
        .{ 0.020435, 0.049183 },
        .{ 0.010278, 0.035401 },
        .{ 0.007569, 0.031203 },
        .{ 0.005124, 0.026842 },
        .{ 0.003235, 0.022872 },
        .{ 0.001625, 0.018778 },
        .{ 0.000631, 0.015282 },
        .{ 0.000000, 0.011694 },
        .{ 0.000088, 0.008846 },
        .{ 0.000351, 0.006020 },
        .{ 0.001225, 0.003898 },
        .{ 0.002663, 0.002099 },
        .{ 0.004054, 0.001048 },
        .{ 0.005660, 0.000389 },
        .{ 0.007649, 0.000031 },
        .{ 0.009663, 0.000000 },
        .{ 0.012095, 0.000289 },
        .{ 0.014480, 0.000826 },
        .{ 0.017259, 0.001730 },
        .{ 0.019942, 0.002874 },
        .{ 0.025860, 0.006081 },
        .{ 0.039158, 0.014055 },
        .{ 0.052488, 0.021970 },
        .{ 0.065205, 0.029003 },
        .{ 0.078096, 0.035710 },
        .{ 0.090332, 0.041618 },
        .{ 0.102755, 0.047127 },
        .{ 0.114539, 0.051907 },
        .{ 0.126483, 0.056272 },
        .{ 0.137856, 0.059993 },
        .{ 0.149354, 0.063309 },
        .{ 0.160399, 0.066097 },
        .{ 0.171531, 0.068513 },
        .{ 0.182359, 0.070508 },
        .{ 0.193243, 0.072168 },
        .{ 0.203964, 0.073488 },
        .{ 0.214716, 0.074522 },
        .{ 0.225362, 0.075276 },
        .{ 0.236021, 0.075789 },
        .{ 0.257060, 0.076098 },
        .{ 0.277659, 0.075493 },
        .{ 0.287634, 0.074850 },
        .{ 0.297594, 0.074004 },
        .{ 0.307184, 0.072977 },
        .{ 0.316734, 0.071626 },
        .{ 0.335108, 0.068196 },
        .{ 0.353287, 0.063847 },
        .{ 0.373463, 0.058063 },
        .{ 0.393279, 0.051654 },
        .{ 0.412745, 0.044491 },
        .{ 0.431923, 0.036586 },
        .{ 0.441443, 0.032270 },
        .{ 0.450826, 0.027661 },
        .{ 0.460367, 0.022685 },
        .{ 0.469785, 0.017480 },
        .{ 0.479443, 0.011805 },
        .{ 0.488953, 0.005886 },
        .{ 0.498720, -0.000528 },
        .{ 0.508326, -0.007180 },
        .{ 0.518154, -0.014344 },
        .{ 0.527832, -0.021709 },
        .{ 0.537943, -0.029707 },
        .{ 0.547860, -0.037944 },
        .{ 0.557378, -0.046301 },
        .{ 0.566712, -0.054864 },
        .{ 0.576342, -0.063994 },
        .{ 0.585810, -0.073291 },
        .{ 0.595251, -0.082951 },
        .{ 0.604502, -0.092792 },
        .{ 0.613686, -0.102927 },
        .{ 0.622688, -0.113225 },
        .{ 0.640270, -0.134439 },
        .{ 0.657151, -0.156283 },
        .{ 0.673265, -0.178579 },
        .{ 0.688595, -0.201154 },
        .{ 0.703073, -0.223816 },
        .{ 0.716738, -0.246414 },
        .{ 0.729583, -0.268769 },
        .{ 0.741636, -0.290734 },
        .{ 0.752938, -0.312161 },
        .{ 0.763519, -0.332903 },
        .{ 0.773423, -0.352847 },
        .{ 0.782696, -0.371894 },
        .{ 0.791327, -0.389909 },
        .{ 0.799366, -0.406804 },
        .{ 0.806821, -0.422490 },
        .{ 0.813669, -0.436886 },
        .{ 0.819910, -0.449898 },
        .{ 0.825496, -0.461463 },
        .{ 0.830405, -0.471509 },
        .{ 0.834580, -0.479967 },
        .{ 0.837962, -0.486768 },
        .{ 0.840543, -0.491860 },
        .{ 0.842197, -0.495145 },
        .{ 0.843535, -0.497014 },
        .{ 0.845263, -0.498560 },
        .{ 0.846593, -0.499318 },
        .{ 0.848086, -0.499770 },
        .{ 0.849626, -0.499988 },
        .{ 0.852022, -0.499894 },
    };
    const chord = 198.0 * 1e-3; // [m]
    const suction_surface_length = 264.7 * 1e-3; // [m]
    const pressure_surface_length = 230.0 * 1e-3; // [m]

    for (&pts) |*p| {
        p[0] *= chord;
        p[1] *= chord;
    }

    var spline = try FittingSpline(dim).init(allocator, pts[0..], degree);
    defer spline.deinit();

    const length = spline.integrate();

    try std.testing.expectApproxEqAbs(suction_surface_length + pressure_surface_length, length, 1e-2);
}
