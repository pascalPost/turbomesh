const std = @import("std");
const types = @import("types.zig");

const Float = types.Float;
const Index = types.Index;
const Vec2d = types.Vec2d;
const eql = types.eql;
const eqlApprox = types.eqlApprox;
const Mat2d = types.Mat2d;
const add = types.add;
const addAll = types.addAll;
const sub = types.sub;
const scale = types.scale;

/// linear TFI as described in chapter 3.5.1 of Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
/// This function takes the clusterings (mappings into intermediate space)
/// s1 (i_min), s2 (i_max), t1 (j_min), t2 (j_min) and the mappings into
/// physical space x_* and computes the coordinates in the internal of the
/// given 2D field
/// as described in chapter 3.6.5 Boundary-Blended Control Functions
pub fn linearBoundaryBlendedControlFunction(data: *Mat2d, x_i_min: []const Vec2d, x_i_max: []const Vec2d, x_j_min: []const Vec2d, x_j_max: []const Vec2d, s1: []const Float, s2: []const Float, t1: []const Float, t2: []const Float) void {
    // TODO add better error handling!

    const n = x_i_min.len;
    std.debug.assert(x_i_max.len == n);
    std.debug.assert(s1.len == n);
    std.debug.assert(s2.len == n);

    const m = x_j_min.len;
    std.debug.assert(x_j_max.len == m);
    std.debug.assert(t1.len == m);
    std.debug.assert(t2.len == m);

    std.debug.assert(s1[0] == 0);
    std.debug.assert(s1[s1.len - 1] == 1.0);

    std.debug.assert(s2[0] == 0);
    std.debug.assert(s2[s2.len - 1] == 1.0);

    std.debug.assert(t1[0] == 0);
    std.debug.assert(t1[t1.len - 1] == 1.0);

    std.debug.assert(t2[0] == 0);
    std.debug.assert(t2[t2.len - 1] == 1.0);

    var idx: Index = 0;
    var i: Index = 0;

    const tol = 1e-10;

    const x_0_0 = x_i_min[0];
    std.debug.assert(eqlApprox(x_0_0, x_j_min[0], tol));

    const x_n_0 = x_i_min[n - 1];
    std.debug.assert(eqlApprox(x_n_0, x_j_max[0], tol));

    const x_0_m = x_j_min[m - 1];
    std.debug.assert(eqlApprox(x_0_m, x_i_max[0], tol));

    const x_n_m = x_i_max[n - 1];
    std.debug.assert(eqlApprox(x_n_m, x_j_max[m - 1], tol));

    while (i < data.size[0]) : (i += 1) {
        const s1_i = s1[i];
        const s2_i = s2[i];

        std.debug.assert(s1_i >= 0.0 and s1_i <= 1.0);
        std.debug.assert(s2_i >= 0.0 and s2_i <= 1.0);

        const x_i_0 = x_i_min[i];
        const x_i_m = x_i_max[i];

        var j: Index = 0;
        while (j < data.size[1]) : (j += 1) {
            const t1_j = t1[j];
            const t2_j = t2[j];

            std.debug.assert(t1_j >= 0.0 and t1_j <= 1.0);
            std.debug.assert(t2_j >= 0.0 and t2_j <= 1.0);

            const x_0_j = x_j_min[j];
            const x_n_j = x_j_max[j];

            const u: Float = ((1.0 - t1_j) * s1_i + t1_j * s2_i) / (1.0 - (s2_i - s1_i) * (t2_j - t1_j));
            const v: Float = ((1.0 - s1_i) * t1_j + s1_i * t2_j) / (1.0 - (t2_j - t1_j) * (s2_i - s1_i));

            const u_ij: Vec2d = add(scale(1.0 - u, x_0_j), scale(u, x_n_j));
            const v_ij: Vec2d = add(scale(1.0 - v, x_i_0), scale(v, x_i_m));
            const uv_ij: Vec2d = addAll(&[_]Vec2d{
                scale(u * v, x_n_m),
                scale(u * (1.0 - v), x_n_0),
                scale((1.0 - u) * v, x_0_m),
                scale((1.0 - u) * (1.0 - v), x_0_0),
            });

            const x_ij: Vec2d = sub(add(u_ij, v_ij), uv_ij);

            std.debug.assert(!std.math.isNan(x_ij.data[0]) and !std.math.isNan(x_ij.data[1]));

            data.data[idx] = x_ij;

            idx += 1;
        }
    }

    std.debug.assert(idx == data.data.len);
}

fn arclengthControlFunction(control_fn: []Float, edge: []const Vec2d) void {
    std.debug.assert(control_fn.len == edge.len);

    control_fn[0] = 0.0;

    var s: Float = 0.0;
    var i: Index = 1;
    while (i < edge.len) : (i += 1) {
        const x_im1 = edge[i - 1];
        const x_i = edge[i];
        const ds = @sqrt((x_i[0] - x_im1[0]) * (x_i[0] - x_im1[0]) + (x_i[1] - x_im1[1]) * (x_i[1] - x_im1[1]));
        s += ds;
        control_fn[i] = s;
    }

    for (control_fn[1..]) |*v| {
        v.* /= s;
    }
}

// test "tfi" {
//     const edge_i_min = [_]Vec2d{ .{ 0.0, 0.0 }, .{ 0.5, 0.0 }, .{ 1.0, 0.0 } };
//     const edge_i_max = [_]Vec2d{ .{ 0.0, 2.0 }, .{ 0.5, 2.0 }, .{ 1.0, 2.0 } };
//     const edge_j_min = [_]Vec2d{ .{ 0.0, 0.0 }, .{ 0.0, 1.0 }, .{ 0.0, 2.0 } };
//     const edge_j_max = [_]Vec2d{ .{ 1.0, 0.0 }, .{ 1.0, 1.0 }, .{ 1.0, 2.0 } };

//     // compute clustering
//     var cf_i_min: [edge_i_min.len]Float = undefined;
//     var cf_i_max: [edge_i_max.len]Float = undefined;
//     var cf_j_min: [edge_j_min.len]Float = undefined;
//     var cf_j_max: [edge_j_max.len]Float = undefined;

//     arclengthControlFunction(&cf_i_min, &edge_i_min);
//     arclengthControlFunction(&cf_i_max, &edge_i_max);
//     arclengthControlFunction(&cf_j_min, &edge_j_min);
//     arclengthControlFunction(&cf_j_max, &edge_j_max);

//     try std.testing.expect(std.mem.eql(Float, &cf_i_min, &[_]Float{ 0.0, 0.5, 1.0 }));
//     try std.testing.expect(std.mem.eql(Float, &cf_i_max, &[_]Float{ 0.0, 0.5, 1.0 }));
//     try std.testing.expect(std.mem.eql(Float, &cf_j_min, &[_]Float{ 0.0, 0.5, 1.0 }));
//     try std.testing.expect(std.mem.eql(Float, &cf_j_max, &[_]Float{ 0.0, 0.5, 1.0 }));

//     // run tfi
//     const allocator = std.testing.allocator;
//     var data = try Mat2d.init(allocator, .{ edge_i_min.len, edge_j_min.len });
//     defer data.deinit(allocator);

//     linearBoundaryBlendedControlFunction(&data, &edge_i_min, &edge_i_max, &edge_j_min, &edge_j_max, &cf_i_min, &cf_i_max, &cf_j_min, &cf_j_max);

//     // check result
// }
