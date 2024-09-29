const std = @import("std");
const types = @import("types.zig");

const Vec2d = types.Vec2d;
const Mat2d = types.Mat2d;

// /// helper function for tfi_linear_2d.
// /// returns a Vec containing the given clustering in computational space
// /// (u in [0,1]).
// fn copy_and_rescale(u: &[Scalar]) -> Vec<Scalar> {
// let mut u = u.to_vec();
//
// let scale = *u.last().unwrap();
// u.iter_mut().for_each(|u| *u /= scale);
// u
// }

fn check_edge_clustering() void {

}

/// linear TFI as described in chapter 3.5.1 of Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
/// This function takes the clusterings (mappings into intermediate space)
/// s1 (i_min), s2 (i_max), t1 (j_min), t2 (j_min) and the mappings into
/// physical space x_* and computes the coordinates in the internal of the
/// given 2D field
/// as described in chapter 3.6.5 Boundary-Blended Control Functions
fn tfi_linear_boundary_blended_control_function(data: Mat2d, s1: []Vec2d, s2: []Vec2d, t1: []Vec2d, t2: []Vec2d) void {
    var i = 0;
    var j = 0;

    while (i < data.size[0]) : (i += 1) {
        const s1_i = s1[i];
        const s2_i = s2[i];

        std.debug.assert(s1_i  >= 0.0 and s1_i <= 1.0);
        std.debug.assert(s2_i  >= 0.0 and s2_i <= 1.0);

        while (j < data.size[1]) : (j += 1) {
            const t1_j = t1[j];
            const t2_j = t2[j];

            std.debug.assert(t1_j  >= 0.0 and t1_j <= 1.0);
            std.debug.assert(t2_j  >= 0.0 and t2_j <= 1.0);

            const u = ((1.0 - t1_j) * s1_i + t1_j * s2_i) / (1.0 - (s2_i - s1_i) * (t2_j - t1_j));
            const v = ((1.0 - s1_i) * t1_j + s1_i * t2_j) / (1.0 - (t2_j - t1_j) * (s2_i - s1_i));
        }
    }
}
