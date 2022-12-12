// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Array2d, Block2d, Vec2d};

/// Boundary-Blended Control Functions eq. (3.21) from Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
fn boundary_blended_control_function(uv: &mut Array2d<Vec2d>) {
    // loop over internal block coordinates
    let i_len = uv.shape.0;
    let j_len = uv.shape.1;

    for i in 1..i_len - 1 {
        let s1 = uv[[i, 0]].0;
        let s2 = uv[[i, j_len - 1]].0;

        for j in 1..j_len - 1 {
            let t1 = uv[[0, j]].1;
            let t2 = uv[[i_len - 1, j]].1;

            let Vec2d(u, v) = &mut uv[[i, j]];
            *u = ((1.0 - t1) * s1 + t1 * s2) / (1.0 - (s2 - s1) * (t2 - t1));
            *v = ((1.0 - s1) * t1 + s1 * t2) / (1.0 - (t2 - t1) * (s2 - s1));
        }
    }
}

/// arclength control function, see section 3.6.4, Thompson et
/// al., Eds., Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
fn arclength_control_function(uv: &mut Array2d<Vec2d>, xy: &Array2d<Vec2d>) {
    // loop over boundary nodes
    let i_len = xy.shape.0;
    let j_len = xy.shape.1;

    assert_eq!(i_len, uv.shape.0);
    assert_eq!(j_len, uv.shape.1);

    {
        for i in 1..i_len {
            // bottom
            uv[[i, 0]].0 = uv[[i - 1, 0]].0 + (xy[[i, 0]] - xy[[i - 1, 0]]).abs();

            // top
            uv[[i, j_len - 1]].0 =
                uv[[i - 1, j_len - 1]].0 + (xy[[i, j_len - 1]] - xy[[i - 1, j_len - 1]]).abs();
        }

        let u_fac_bottom = uv[[i_len - 1, 0]].0;
        let u_fac_top = uv[[i_len - 1, j_len - 1]].0;

        for i in 1..i_len {
            uv[[i, 0]].0 /= u_fac_bottom;

            uv[[i, j_len - 1]].0 /= u_fac_top;
            uv[[i, j_len - 1]].1 = 1.0;
        }
    }

    {
        for j in 1..j_len {
            // left
            uv[[0, j]].1 = uv[[0, j - 1]].1 + (xy[[0, j]] - xy[[0, j - 1]]).abs();

            // right
            uv[[i_len - 1, j]].1 =
                uv[[i_len - 1, j - 1]].1 + (xy[[i_len - 1, j]] - xy[[i_len - 1, j - 1]]).abs();
        }

        let v_fac_left = uv[[0, j_len - 1]].1;
        let v_fac_right = uv[[i_len - 1, j_len - 1]].1;

        for j in 1..j_len {
            uv[[0, j]].1 /= v_fac_left;

            uv[[i_len - 1, j]].0 = 1.0;
            uv[[i_len - 1, j]].1 /= v_fac_right;
        }
    }
}

/// creation of the intermediate control domain, see section 3.6, Thompson et
/// al., Eds., Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
fn intermediate_control_domain(xy: &Array2d<Vec2d>) -> Array2d<Vec2d> {
    let mut uv = Array2d::new(xy.shape);

    arclength_control_function(&mut uv, xy);
    boundary_blended_control_function(&mut uv);

    uv
}

/// linear TFI as described in chapter 3.5.1 of Thompson et al., Eds.,
/// Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.
///
/// The coordinates of the boundaries of the block must be set before
/// calling this function. On this basis, the coordinates inside the block
/// are computed.
pub fn tfi_linear_2d(block: &mut Block2d) {
    let uv = intermediate_control_domain(&block.coords);
    let xy = &mut block.coords;

    // loop over internal block coordinates
    let i_len = xy.shape.0;
    let j_len = xy.shape.1;

    for i in 1..i_len - 1 {
        for j in 1..j_len - 1 {
            let Vec2d(xi, eta) = uv[[i, j]];

            // TODO can be optimized by reducing memory accesses

            let u_ij = (1.0 - xi) * xy[[0, j]] + xi * xy[[i_len - 1, j]];
            let v_ij = (1.0 - eta) * xy[[i, 0]] + eta * xy[[i, j_len - 1]];
            let uv_ij = xi * eta * xy[[i_len - 1, j_len - 1]]
                + xi * (1.0 - eta) * xy[[i_len - 1, 0]]
                + (1.0 - xi) * eta * xy[[0, j_len - 1]]
                + (1.0 - xi) * (1.0 - eta) * xy[[0, 0]];

            xy[[i, j]] = u_ij + v_ij - uv_ij
        }
    }
}
