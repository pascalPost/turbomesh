// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Scalar, Vec2d};
use float_cmp::approx_eq;
use ndarray::Array2;

// /// Boundary-Blended Control Functions eq. (3.21) from Thompson et al., Eds.,
// /// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
// pub fn boundary_blended_control_function(uv: &mut Array2d<Vec2d>) {
//     // loop over internal block coordinates
//     let [i_len, j_len] = uv.shape;

//     for i in 1..i_len - 1 {
//         let s1 = uv[[i, 0]].0;
//         let s2 = uv[[i, j_len - 1]].0;

//         for j in 1..j_len - 1 {
//             let t1 = uv[[0, j]].1;
//             let t2 = uv[[i_len - 1, j]].1;

//             let Vec2d(u, v) = &mut uv[[i, j]];
//             *u = ((1.0 - t1) * s1 + t1 * s2) / (1.0 - (s2 - s1) * (t2 - t1));
//             *v = ((1.0 - s1) * t1 + s1 * t2) / (1.0 - (t2 - t1) * (s2 - s1));
//         }
//     }
// }

// /// arclength control function, see section 3.6.4, Thompson et
// /// al., Eds., Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
// fn arclength_control_function(uv: &mut Array2d<Vec2d>, xy: &Array2d<Vec2d>) {
//     // loop over boundary nodes
//     let [i_len, j_len] = xy.shape;

//     assert_eq!(
//         [i_len, j_len],
//         uv.shape,
//         "shape mismatch in mapping detected."
//     );

//     {
//         for i in 1..i_len {
//             // bottom
//             uv[[i, 0]].0 = uv[[i - 1, 0]].0 + (xy[[i, 0]] - xy[[i - 1, 0]]).abs();

//             // top
//             uv[[i, j_len - 1]].0 =
//                 uv[[i - 1, j_len - 1]].0 + (xy[[i, j_len - 1]] - xy[[i - 1, j_len - 1]]).abs();
//         }

//         let u_fac_bottom = uv[[i_len - 1, 0]].0;
//         let u_fac_top = uv[[i_len - 1, j_len - 1]].0;

//         for i in 1..i_len {
//             uv[[i, 0]].0 /= u_fac_bottom;

//             uv[[i, j_len - 1]].0 /= u_fac_top;
//             uv[[i, j_len - 1]].1 = 1.0;
//         }
//     }

//     {
//         for j in 1..j_len {
//             // left
//             uv[[0, j]].1 = uv[[0, j - 1]].1 + (xy[[0, j]] - xy[[0, j - 1]]).abs();

//             // right
//             uv[[i_len - 1, j]].1 =
//                 uv[[i_len - 1, j - 1]].1 + (xy[[i_len - 1, j]] - xy[[i_len - 1, j - 1]]).abs();
//         }

//         let v_fac_left = uv[[0, j_len - 1]].1;
//         let v_fac_right = uv[[i_len - 1, j_len - 1]].1;

//         for j in 1..j_len {
//             uv[[0, j]].1 /= v_fac_left;

//             uv[[i_len - 1, j]].0 = 1.0;
//             uv[[i_len - 1, j]].1 /= v_fac_right;
//         }
//     }
// }

// /// creation of the intermediate control domain, see section 3.6, Thompson et
// /// al., Eds., Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
// fn intermediate_control_domain(xy: &Array2d<Vec2d>) -> Array2d<Vec2d> {
//     let mut uv = Array2d::new(xy.shape);

//     arclength_control_function(&mut uv, xy);
//     boundary_blended_control_function(&mut uv);

//     uv
// }

fn check_edge_clustering(u: &[Scalar], name: &str) {
    for (i, u) in u.iter().enumerate() {
        assert!(!u.is_nan(), "The clustering for edge {name} contains a NAN at i = {i}. Check that this value was properly set.");
    }

    assert_eq!(
        u.first(),
        Some(&0.0),
        "The clustering for edge {name} does not satisfy u[0] == 0.0."
    );
    assert_eq!(
        u.last(),
        Some(&1.0),
        "The clustering for edge {name} does not satisfy u[-1] == 1.0."
    );
    assert!(
        u.windows(2).all(|u| u[0] < u[1]),
        "The clustering for edge {name} does not satisfy u[i] < u[i+1]."
    );
}

fn check_edge_mapping(v: &[Vec2d], name: &str) {
    v.iter().enumerate().for_each(|(i, Vec2d(x, y))| {
        assert!(
            !x.is_nan() || !y.is_nan(),
            "The edge {name} contains a NAN at point index i = {i} with x = ({x}, {y})"
        );
    });
}

fn check_and_get_corners(
    coords_i_min: &[Vec2d],
    coords_i_max: &[Vec2d],
    coords_j_min: &[Vec2d],
    coords_j_max: &[Vec2d],
) -> (Vec2d, Vec2d, Vec2d, Vec2d) {
    // const ULPS: i64 = 100;
    const EPSILON: f64 = 1e-10;

    // x[0,0]
    let x_0_0 = *coords_i_min.first().unwrap();
    assert!(
        approx_eq!(
            &Vec2d,
            &x_0_0,
            coords_j_min.first().unwrap(),
            epsilon = EPSILON // ulps = ULPS
        ),
        "The corners of edges i_min[0] = {} and 
            j_min[0] = {} bound to constituting x[0, 0] do not match.",
        x_0_0,
        coords_j_min.first().unwrap()
    );

    // x[-1, -1]
    let x_1_1 = *coords_i_max.last().unwrap();
    assert!(
        approx_eq!(
            &Vec2d,
            &x_1_1,
            coords_j_max.last().unwrap(),
            epsilon = EPSILON // ulps = ULPS
        ),
        "The corners of edges i_max[-1] = {} and
            j_max[-1] = {} bound to constitute x[-1, -1] do not match.",
        x_1_1,
        coords_j_max.last().unwrap()
    );

    // x[-1, 0]
    let x_1_0 = *coords_i_min.last().unwrap();
    assert!(
        approx_eq!(
            &Vec2d,
            &x_1_0,
            coords_j_max.first().unwrap(),
            epsilon = EPSILON // ulps = ULPS
        ),
        "The corners of edges i_min[-1] = {} and 
            j_max[0] = {} bound to constitute x[-1, 0] do not match.",
        x_1_0,
        coords_j_max.first().unwrap()
    );

    // x[0, -1]
    let x_0_1 = *coords_i_max.first().unwrap();
    assert!(
        approx_eq!(
            &Vec2d,
            &x_0_1,
            coords_j_min.last().unwrap(),
            epsilon = EPSILON // ulps = ULPS
        ),
        "The corners of edges i_max[0] = {} and 
            j_min[-1] = {} bound to constitute x[0, -1] do not match.",
        x_0_1,
        coords_j_min.last().unwrap()
    );

    (x_0_0, x_1_0, x_0_1, x_1_1)
}

/// linear TFI as described in chapter 3.5.1 of Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
/// This function takes the clusterings (mappings into intermediate space)
/// s1 (i_min), s2 (i_max), t1 (j_min), t2 (j_min) and the mappings into
/// physical space x_* and computes the coordinates in the internal of the
/// given 2D field
pub fn tfi_linear_2d(
    s1: &[Scalar],
    s2: &[Scalar],
    t1: &[Scalar],
    t2: &[Scalar],
    x_i_min: &[Vec2d],
    x_i_max: &[Vec2d],
    x_j_min: &[Vec2d],
    x_j_max: &[Vec2d],
    x: &mut Array2<Vec2d>,
) {
    let s1 = copy_and_rescale(s1);
    let s2 = copy_and_rescale(s2);
    let t1 = copy_and_rescale(t1);
    let t2 = copy_and_rescale(t2);

    check_edge_clustering(&s1, "i_min");
    check_edge_clustering(&s2, "i_max");
    check_edge_clustering(&t1, "j_min");
    check_edge_clustering(&t2, "j_max");

    check_edge_mapping(&x_i_min, "i_min");
    check_edge_mapping(&x_i_max, "i_max");
    check_edge_mapping(&x_j_min, "j_min");
    check_edge_mapping(&x_j_max, "j_max");

    let (x_0_0, x_1_0, x_0_1, x_1_1) =
        check_and_get_corners(&x_i_min, &x_i_max, &x_j_min, &x_j_max);

    // TODO switch to iterator syntax to get independent from mem layout and
    // better suited for parallel exec

    let (i_len, j_len) = x.dim();

    // TODO refactor looping for higher efficiency
    for j in 0..j_len {
        let t1 = t1[j];
        let t2 = t2[j];

        let x_0_eta = x_j_min[j];
        let x_1_eta = x_j_max[j];

        for i in 0..i_len {
            let s1 = s1[i];
            let s2 = s2[i];

            let x_xi_0 = x_i_min[i];
            let x_xi_1 = x_i_max[i];

            let xi = ((1.0 - t1) * s1 + t1 * s2) / (1.0 - (s2 - s1) * (t2 - t1));
            let eta = ((1.0 - s1) * t1 + s1 * t2) / (1.0 - (t2 - t1) * (s2 - s1));

            let u_ij = (1.0 - xi) * x_0_eta + xi * x_1_eta;
            let v_ij = (1.0 - eta) * x_xi_0 + eta * x_xi_1;
            let uv_ij = xi * eta * x_1_1
                + xi * (1.0 - eta) * x_1_0
                + (1.0 - xi) * eta * x_0_1
                + (1.0 - xi) * (1.0 - eta) * x_0_0;

            x[[i, j]] = u_ij + v_ij - uv_ij
        }
    }

    /// helper function for tfi_linear_2d.
    /// returns a Vec containing the given clustering in computational space
    /// (u in [0,1]).
    fn copy_and_rescale(u: &[Scalar]) -> Vec<Scalar> {
        let mut u = u.to_vec();

        let scale = *u.last().unwrap();
        u.iter_mut().for_each(|u| *u /= scale);
        u
    }
}

// TODO if there is something like constexpr this function can be merged with the version including an intermediate grid
pub fn tfi_linear_2d_simple(x: &mut Array2<Vec2d>) {
    let (m, n) = x.dim();

    let x_0_0 = x[[0, 0]];
    let x_1_0 = x[[m - 1, 0]];
    let x_0_1 = x[[0, n - 1]];
    let x_1_1 = x[[m - 1, n - 1]];

    // TODO switch to iterator syntax to get independent from mem layout and
    // better suited for parallel exec

    let (i_len, j_len) = x.dim();

    // TODO refactor looping for higher efficiency
    for j in 0..j_len {
        let x_0_eta = x[[0, j]];
        let x_1_eta = x[[m - 1, j]];

        for i in 0..i_len {
            let x_xi_0 = x[[i, 0]];
            let x_xi_1 = x[[i, n - 1]];

            let xi = i as Scalar / (i_len - 1) as Scalar;
            let eta = j as Scalar / (j_len - 1) as Scalar;

            let u_ij = (1.0 - xi) * x_0_eta + xi * x_1_eta;
            let v_ij = (1.0 - eta) * x_xi_0 + eta * x_xi_1;
            let uv_ij = xi * eta * x_1_1
                + xi * (1.0 - eta) * x_1_0
                + (1.0 - xi) * eta * x_0_1
                + (1.0 - xi) * (1.0 - eta) * x_0_0;

            x[[i, j]] = u_ij + v_ij - uv_ij
        }
    }
}
