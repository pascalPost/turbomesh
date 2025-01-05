const std = @import("std");
const types = @import("types.zig");

pub fn block(points: *types.Mat2d) void {
    const dof = (points.size[0] - 2) * (points.size[1] - 2);

    // // allocate solution vectors
    // let mut x_new = Vector::new(dof);
    // let mut y_new = Vector::new(dof);

    // // allocate right hand side containing the current internal coordinates
    // let mut rhs_x = Vector::new(dof);
    // let mut rhs_y = Vector::new(dof);

    // // allocate left hand side as square matrix
    // let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;

    // for n in 0..iterations {
    //     print!("  iteration\t{n}");

    //     lhs.reset();
    //     rhs_x.fill(0.0);
    //     rhs_y.fill(0.0);

    //     // fill lhs & rhs

    //     // TODO refactor looping for higher efficiency
    //     for j in 1..shape.1 - 1 {
    //         for i in 1..shape.0 - 1 {
    //             // laplace conditions (zero intialization)
    //             let s = 0.0;
    //             let t = 0.0;

    //             // let Vec2d(x_i_j, y_i_j) = coords[[i, j]];
    //             let Vec2d(x_ip1_j, y_ip1_j) = coords[[i + 1, j]];
    //             let Vec2d(x_im1_j, y_im1_j) = coords[[i - 1, j]];
    //             let Vec2d(x_i_jp1, y_i_jp1) = coords[[i, j + 1]];
    //             let Vec2d(x_i_jm1, y_i_jm1) = coords[[i, j - 1]];

    //             let Vec2d(x_ip1_jp1, y_ip1_jp1) = coords[[i + 1, j + 1]];
    //             let Vec2d(x_im1_jm1, y_im1_jm1) = coords[[i - 1, j - 1]];
    //             let Vec2d(x_ip1_jm1, y_ip1_jm1) = coords[[i + 1, j - 1]];
    //             let Vec2d(x_im1_jp1, y_im1_jp1) = coords[[i - 1, j + 1]];

    //             let x_xi = 0.5 * (x_ip1_j - x_im1_j);
    //             let x_eta = 0.5 * (x_i_jp1 - x_i_jm1);
    //             let y_xi = 0.5 * (y_ip1_j - y_im1_j);
    //             let y_eta = 0.5 * (y_i_jp1 - y_i_jm1);

    //             let p = x_eta * x_eta + y_eta * y_eta;
    //             let q = x_xi * x_eta + y_xi * y_eta;
    //             let r = x_xi * x_xi + y_xi * y_xi;

    //             let a_i_j = -2.0 * p - 2.0 * r;
    //             let a_ip1_j = p + 0.5 * s;
    //             let a_im1_j = p - 0.5 * s;
    //             let a_i_jp1 = r + 0.5 * t;
    //             let a_i_jm1 = r - 0.5 * t;
    //             let a_ip1_jp1 = -0.5 * q;
    //             let a_ip1_jm1 = 0.5 * q;
    //             let a_im1_jp1 = 0.5 * q;
    //             let a_im1_jm1 = -0.5 * q;

    //             let idx = matrix_index(i, j, i_size);

    //             lhs.put(idx, idx, a_i_j)?;

    //             // boundary values
    //             if i == 1 {
    //                 rhs_x[idx] -= a_im1_j * x_im1_j;
    //                 rhs_y[idx] -= a_im1_j * y_im1_j;
    //             } else {
    //                 lhs.put(idx, matrix_index(i - 1, j, i_size), a_im1_j)?;
    //             }

    //             if i == i_size - 2 {
    //                 rhs_x[idx] -= a_ip1_j * x_ip1_j;
    //                 rhs_y[idx] -= a_ip1_j * y_ip1_j;
    //             } else {
    //                 lhs.put(idx, matrix_index(i + 1, j, i_size), a_ip1_j)?;
    //             }

    //             if j == 1 {
    //                 rhs_x[idx] -= a_i_jm1 * x_i_jm1;
    //                 rhs_y[idx] -= a_i_jm1 * y_i_jm1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i, j - 1, i_size), a_i_jm1)?;
    //             }

    //             if j == j_size - 2 {
    //                 rhs_x[idx] -= a_i_jp1 * x_i_jp1;
    //                 rhs_y[idx] -= a_i_jp1 * y_i_jp1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i, j + 1, i_size), a_i_jp1)?;
    //             }

    //             if i == i_size - 2 || j == j_size - 2 {
    //                 rhs_x[idx] -= a_ip1_jp1 * x_ip1_jp1;
    //                 rhs_y[idx] -= a_ip1_jp1 * y_ip1_jp1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i + 1, j + 1, i_size), a_ip1_jp1)?;
    //             }

    //             if i == i_size - 2 || j == 1 {
    //                 rhs_x[idx] -= a_ip1_jm1 * x_ip1_jm1;
    //                 rhs_y[idx] -= a_ip1_jm1 * y_ip1_jm1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i + 1, j - 1, i_size), a_ip1_jm1)?;
    //             }

    //             if i == 1 || j == j_size - 2 {
    //                 rhs_x[idx] -= a_im1_jp1 * x_im1_jp1;
    //                 rhs_y[idx] -= a_im1_jp1 * y_im1_jp1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i - 1, j + 1, i_size), a_im1_jp1)?;
    //             }

    //             if i == 1 || j == 1 {
    //                 rhs_x[idx] -= a_im1_jm1 * x_im1_jm1;
    //                 rhs_y[idx] -= a_im1_jm1 * y_im1_jm1;
    //             } else {
    //                 lhs.put(idx, matrix_index(i - 1, j - 1, i_size), a_im1_jm1)?;
    //             }
    //         }
    //     }

    //     // solve
    //     let mut solver = Solver::new(config)?;
    //     solver.initialize(&lhs)?;
    //     solver.factorize()?;
    //     solver.solve(&mut x_new, &rhs_x)?;
    //     solver.solve(&mut y_new, &rhs_y)?;

    //     let mut x_norm_sqr = 0.0;
    //     let mut y_norm_sqr = 0.0;

    //     for j in 1..shape.1 - 1 {
    //         for i in 1..shape.0 - 1 {
    //             let idx = matrix_index(i, j, i_size);
    //             let Vec2d(x, y) = coords[[i, j]];
    //             let dx = x - x_new[idx];
    //             let dy = y - y_new[idx];

    //             x_norm_sqr += dx * dx;
    //             y_norm_sqr += dy * dy;
    //         }
    //     }

    //     let norm = (x_norm_sqr + y_norm_sqr).sqrt();

    //     println!("\tresidual {norm}");

    //     // copy into coordinate field
    //     // TODO consider different memory layout for higher efficiency
    //     for j in 1..shape.1 - 1 {
    //         for i in 1..shape.0 - 1 {
    //             let idx = matrix_index(i, j, i_size);
    //             coords[[i, j]] = Vec2d(x_new[idx], y_new[idx]);
    //         }
    //     }
    // }

    // Ok(())
}
