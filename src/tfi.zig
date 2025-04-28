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
const negate = types.negate;

pub fn linear2d(
    data: []Vec2d,
    edge_i_min: []const Vec2d,
    edge_i_max: []const Vec2d,
    edge_j_min: []const Vec2d,
    edge_j_max: []const Vec2d,
) !void {

    // TODO: make this either comptime t or allow to do it for n consecutive float

    const size = types.Index2d{ edge_i_min.len, edge_j_min.len };
    if (edge_i_max.len != size[0] or edge_j_max.len != size[1]) return error.InconsistentSize;

    // NOTE: corner values are taken from the i edges;
    // consistency with the values on the j edges is not checked.
    const corner_0_0 = edge_i_min[0];
    const corner_1_0 = edge_i_min[size[0] - 1];
    const corner_0_1 = edge_i_max[0];
    const corner_1_1 = edge_i_max[size[0] - 1];

    for (0..size[0]) |i| {
        const v_xi_0 = edge_i_min[i];
        const v_xi_1 = edge_i_max[i];

        const xi: Float = @as(Float, @floatFromInt(i)) / @as(Float, @floatFromInt(size[0] - 1));

        for (0..size[1]) |j| {
            const v_0_eta = edge_j_min[j];
            const v_1_eta = edge_j_max[j];

            const eta: Float = @as(Float, @floatFromInt(j)) / @as(Float, @floatFromInt(size[1] - 1));

            const u_ij = add(scale(1.0 - xi, v_0_eta), scale(xi, v_1_eta));
            const v_ij = add(scale(1.0 - eta, v_xi_0), scale(eta, v_xi_1));
            const uv_ij = addAll(.{
                scale(xi * eta, corner_1_1),
                scale(xi * (1.0 - eta), corner_1_0),
                scale((1.0 - xi) * eta, corner_0_1),
                scale((1.0 - xi) * (1.0 - eta), corner_0_0),
            });

            data[i * size[1] + j] = addAll(.{
                u_ij,
                v_ij,
                negate(uv_ij),
            }); // (i,j)
        }
    }
}

// /// Applies the linear tfi based on the values set on the edges of the data block
// /// with the given size. A column major memory layout is assumed.
// /// The linear TFI is described in chapter 3.5.1 of Thompson et al., Eds.,
// /// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
// pub fn linear2dInPlace(data: []types.Vec2d, size: types.Index2d) void {
//     // NOTE: we assume column major layout.
//
//     // corner values
//     const v_0_0 = data[0]; // (0,0)
//     const v_0_1 = data[size[1] - 1]; // (0,J-1)
//     const v_1_0 = data[(size[0] - 1) * size[1]]; // (I-1,0)
//     const v_1_1 = data[(size[0] - 1) * size[1] + size[1] - 1]; // (I-1,J-1)
//
//     const im1: Float = @floatFromInt(size[0] - 1);
//     const jm1: Float = @floatFromInt(size[1] - 1);
//
//     for (0..size[0]) |i| {
//         const v_xi_0 = data[i * size[1]]; // (i,0)
//         const v_xi_1 = data[i * size[1] + size[1] - 1]; // (i,J-1)
//
//         for (0..size[1]) |j| {
//             const v_0_eta = data[j]; // (0,j)
//             const v_1_eta = data[(size[0] - 1) * size[0] + j]; // (I-1,j)
//
//             const xi: Float = @as(Float, @floatFromInt(i)) / im1;
//             const eta: Float = @as(Float, @floatFromInt(j)) / jm1;
//
//             const u_ij = (1.0 - xi) * v_0_eta + xi * v_1_eta;
//             const v_ij = (1.0 - eta) * v_xi_0 + eta * v_xi_1;
//             const uv_ij = xi * eta * v_1_1 + xi * (1.0 - eta) * v_1_0 + (1.0 - xi) * eta * v_0_1 + (1.0 - xi) * (1.0 - eta) * v_0_0;
//
//             data[i * size[1] + j] = u_ij + v_ij - uv_ij; // (i,j)
//         }
//     }
// }

/// linear TFI as described in chapter 3.5.1 of Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
/// This function takes the clusterings (mappings into intermediate space)
/// s1 (i_min), s2 (i_max), t1 (j_min), t2 (j_min) and the mappings into
/// physical space x_* and computes the coordinates in the internal of the
/// given 2D field
/// as described in chapter 3.6.5 Boundary-Blended Control Functions
pub fn linear2dBoundaryBlendedControlFunction(
    data: *Mat2d,
    x_i_min: []const Vec2d,
    x_i_max: []const Vec2d,
    x_j_min: []const Vec2d,
    x_j_max: []const Vec2d,
    s1: []const Float,
    s2: []const Float,
    t1: []const Float,
    t2: []const Float,
) void {
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
