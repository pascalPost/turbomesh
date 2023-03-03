// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    types::{BlockBoundary, BlockConnection, EdgeIndex},
    Block2d, Mesh, Scalar, Vec2d,
};
use ndarray::Array2;
use russell_lab::{vector_norm, NormVec, Vector};
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::error;

fn matrix_index(i: usize, j: usize, i_size: usize) -> usize {
    (i - 1) + (j - 1) * (i_size - 2)
}

// fn derivatives(
//     x: &Array2d<Scalar>,
//     y: &Array2d<Scalar>,
//     x_xi: &mut Array2d<Scalar>,
//     x_eta: &mut Array2d<Scalar>,
//     y_xi: &mut Array2d<Scalar>,
//     y_eta: &mut Array2d<Scalar>,
// ) {
//     // TODO compute derivatives only on internal points and reduce memory
//     // allocation

//     // TODO enhance memory access

//     for j in 1..x.shape[1] - 1 {
//         for i in 1..x.shape[0] - 1 {
//             x_xi[[i, j]] = 0.5 * (x[[i + 1, j]] - x[[i - 1, j]]);
//             x_eta[[i, j]] = 0.5 * (x[[i, j + 1]] - x[[i, j - 1]]);
//             y_xi[[i, j]] = 0.5 * (y[[i + 1, j]] - y[[i - 1, j]]);
//             y_eta[[i, j]] = 0.5 * (y[[i, j + 1]] - y[[i, j - 1]]);
//         }
//     }
// }

pub fn smooth_block(block: &mut Block2d, iterations: usize) -> Result<(), Box<dyn error::Error>> {
    let coords = &mut block.coords;

    let config = ConfigSolver::new();

    // let [I, J] = X.shape;
    let shape = coords.dim();
    let (i_size, j_size) = shape;

    // // copy x and y coordinate into seperate fields to allow for more effective
    // // computation. To reduce memory footprint, might be done sequentially.
    // let mut x: Vec<Scalar> = block.coords.iter().map(|x| x.0).collect();
    // let mut y: Vec<Scalar> = block.coords.iter().map(|x| x.1).collect();

    // // allocate fields containing derivatives (might be unnecessary?)
    // let x_xi = Array2d::<Scalar>::new(shape);
    // let x_eta = Array2d::<Scalar>::new(shape);
    // let y_xi = Array2d::<Scalar>::new(shape);
    // let y_eta = Array2d::<Scalar>::new(shape);

    // DOF
    let dof = (shape.0 - 2) * (shape.1 - 2);

    // allocate solution vectors
    let mut x_new = Vector::new(dof);
    let mut y_new = Vector::new(dof);

    // allocate right hand side containing the current internal coordinates
    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    // allocate left hand side as square matrix
    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;

    for n in 0..iterations {
        print!("  iteration\t{n}");

        lhs.reset();
        rhs_x.fill(0.0);
        rhs_y.fill(0.0);

        // fill lhs & rhs

        // TODO refactor looping for higher efficiency
        for j in 1..shape.1 - 1 {
            for i in 1..shape.0 - 1 {
                // laplace conditions (zero intialization)
                let s = 0.0;
                let t = 0.0;

                // let Vec2d(x_i_j, y_i_j) = coords[[i, j]];
                let Vec2d(x_ip1_j, y_ip1_j) = coords[[i + 1, j]];
                let Vec2d(x_im1_j, y_im1_j) = coords[[i - 1, j]];
                let Vec2d(x_i_jp1, y_i_jp1) = coords[[i, j + 1]];
                let Vec2d(x_i_jm1, y_i_jm1) = coords[[i, j - 1]];

                let Vec2d(x_ip1_jp1, y_ip1_jp1) = coords[[i + 1, j + 1]];
                let Vec2d(x_im1_jm1, y_im1_jm1) = coords[[i - 1, j - 1]];
                let Vec2d(x_ip1_jm1, y_ip1_jm1) = coords[[i + 1, j - 1]];
                let Vec2d(x_im1_jp1, y_im1_jp1) = coords[[i - 1, j + 1]];

                let x_xi = 0.5 * (x_ip1_j - x_im1_j);
                let x_eta = 0.5 * (x_i_jp1 - x_i_jm1);
                let y_xi = 0.5 * (y_ip1_j - y_im1_j);
                let y_eta = 0.5 * (y_i_jp1 - y_i_jm1);

                let p = x_eta * x_eta + y_eta * y_eta;
                let q = x_xi * x_eta + y_xi * y_eta;
                let r = x_xi * x_xi + y_xi * y_xi;

                let a_i_j = -2.0 * p - 2.0 * r;
                let a_ip1_j = p + 0.5 * s;
                let a_im1_j = p - 0.5 * s;
                let a_i_jp1 = r + 0.5 * t;
                let a_i_jm1 = r - 0.5 * t;
                let a_ip1_jp1 = -0.5 * q;
                let a_ip1_jm1 = 0.5 * q;
                let a_im1_jp1 = 0.5 * q;
                let a_im1_jm1 = -0.5 * q;

                let idx = matrix_index(i, j, i_size);

                lhs.put(idx, idx, a_i_j)?;

                // boundary values
                if i == 1 {
                    rhs_x[idx] -= a_im1_j * x_im1_j;
                    rhs_y[idx] -= a_im1_j * y_im1_j;
                } else {
                    lhs.put(idx, matrix_index(i - 1, j, i_size), a_im1_j)?;
                }

                if i == i_size - 2 {
                    rhs_x[idx] -= a_ip1_j * x_ip1_j;
                    rhs_y[idx] -= a_ip1_j * y_ip1_j;
                } else {
                    lhs.put(idx, matrix_index(i + 1, j, i_size), a_ip1_j)?;
                }

                if j == 1 {
                    rhs_x[idx] -= a_i_jm1 * x_i_jm1;
                    rhs_y[idx] -= a_i_jm1 * y_i_jm1;
                } else {
                    lhs.put(idx, matrix_index(i, j - 1, i_size), a_i_jm1)?;
                }

                if j == j_size - 2 {
                    rhs_x[idx] -= a_i_jp1 * x_i_jp1;
                    rhs_y[idx] -= a_i_jp1 * y_i_jp1;
                } else {
                    lhs.put(idx, matrix_index(i, j + 1, i_size), a_i_jp1)?;
                }

                if i == i_size - 2 || j == j_size - 2 {
                    rhs_x[idx] -= a_ip1_jp1 * x_ip1_jp1;
                    rhs_y[idx] -= a_ip1_jp1 * y_ip1_jp1;
                } else {
                    lhs.put(idx, matrix_index(i + 1, j + 1, i_size), a_ip1_jp1)?;
                }

                if i == i_size - 2 || j == 1 {
                    rhs_x[idx] -= a_ip1_jm1 * x_ip1_jm1;
                    rhs_y[idx] -= a_ip1_jm1 * y_ip1_jm1;
                } else {
                    lhs.put(idx, matrix_index(i + 1, j - 1, i_size), a_ip1_jm1)?;
                }

                if i == 1 || j == j_size - 2 {
                    rhs_x[idx] -= a_im1_jp1 * x_im1_jp1;
                    rhs_y[idx] -= a_im1_jp1 * y_im1_jp1;
                } else {
                    lhs.put(idx, matrix_index(i - 1, j + 1, i_size), a_im1_jp1)?;
                }

                if i == 1 || j == 1 {
                    rhs_x[idx] -= a_im1_jm1 * x_im1_jm1;
                    rhs_y[idx] -= a_im1_jm1 * y_im1_jm1;
                } else {
                    lhs.put(idx, matrix_index(i - 1, j - 1, i_size), a_im1_jm1)?;
                }
            }
        }

        // solve
        let mut solver = Solver::new(config)?;
        solver.initialize(&lhs)?;
        solver.factorize()?;
        solver.solve(&mut x_new, &rhs_x)?;
        solver.solve(&mut y_new, &rhs_y)?;

        let x_norm = vector_norm(&x_new, NormVec::Euc);
        let y_norm = vector_norm(&y_new, NormVec::Euc);

        let norm = (x_norm * x_norm + y_norm * y_norm).sqrt();

        println!("\tresidual {norm}");

        // copy into coordinate field
        // TODO consider different memory layout for higher efficiency
        for j in 1..shape.1 - 1 {
            for i in 1..shape.0 - 1 {
                let idx = matrix_index(i, j, i_size);
                coords[[i, j]] = Vec2d(x_new[idx], y_new[idx]);
            }
        }
    }

    Ok(())
}

pub fn compute_derivatives(mesh: &Mesh) {
    // allocate derivative fields
    let n = mesh.blocks.len();

    let mut x_xi_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
    let mut x_eta_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
    let mut y_xi_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
    let mut y_eta_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);

    mesh.blocks.iter().for_each(|block| {
        let shape = block.coords.dim();
        let x = &block.coords;

        let mut x_xi = Array2::<Scalar>::zeros(shape);
        let mut x_eta = Array2::<Scalar>::zeros(shape);
        let mut y_xi = Array2::<Scalar>::zeros(shape);
        let mut y_eta = Array2::<Scalar>::zeros(shape);

        // loop internal points
        for j in 1..shape.1 - 1 {
            for i in 1..shape.0 - 1 {
                let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
                let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
                let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
                let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

                x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
                x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
                y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
                y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
            }
        }

        // add own contributions on corners
        {
            let j = 0;
            let i = 0;

            let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
            let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

            x_xi[[i, j]] = 0.5 * x_ip1_j;
            x_eta[[i, j]] = 0.5 * x_i_jp1;
            y_xi[[i, j]] = 0.5 * y_ip1_j;
            y_eta[[i, j]] = 0.5 * y_i_jp1;
        }

        {
            let j = shape.1 - 1;
            let i = shape.0 - 1;

            let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
            let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

            x_xi[[i, j]] = -0.5 * x_im1_j;
            x_eta[[i, j]] = -0.5 * x_i_jm1;
            y_xi[[i, j]] = -0.5 * y_im1_j;
            y_eta[[i, j]] = -0.5 * y_i_jm1;
        }

        {
            let j = shape.1 - 1;
            let i = 0;

            let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
            let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

            x_xi[[i, j]] = 0.5 * x_ip1_j;
            x_eta[[i, j]] = -0.5 * x_i_jm1;
            y_xi[[i, j]] = 0.5 * y_ip1_j;
            y_eta[[i, j]] = -0.5 * y_i_jm1;
        }

        {
            let j = 0;
            let i = shape.0 - 1;

            let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
            let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

            x_xi[[i, j]] = -0.5 * x_im1_j;
            x_eta[[i, j]] = 0.5 * x_i_jp1;
            y_xi[[i, j]] = -0.5 * y_im1_j;
            y_eta[[i, j]] = 0.5 * y_i_jp1;
        }

        // add own contributions on edges excluding corners

        // edge i min & i max
        for j in 1..shape.1 - 1 {
            {
                // edge i min
                let i = 0;

                let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
                let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
                let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

                x_xi[[i, j]] = 0.5 * x_ip1_j;
                x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
                y_xi[[i, j]] = 0.5 * y_ip1_j;
                y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
            }

            {
                // edge i max
                let i = shape.0 - 1;

                let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
                let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
                let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

                x_xi[[i, j]] = -0.5 * x_im1_j;
                x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
                y_xi[[i, j]] = -0.5 * y_im1_j;
                y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
            }
        }

        // edge j min & j max
        for i in 1..shape.0 - 1 {
            {
                // edge j min
                let j = 0;

                let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
                let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
                let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

                x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
                x_eta[[i, j]] = 0.5 * x_i_jp1;
                y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
                y_eta[[i, j]] = 0.5 * y_i_jp1;
            }

            {
                // edge j max
                let j = shape.1 - 1;

                let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
                let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
                let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

                x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
                x_eta[[i, j]] = -0.5 * x_i_jm1;
                y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
                y_eta[[i, j]] = -0.5 * y_i_jm1;
            }
        }

        x_xi_mesh.push(x_xi);
        x_eta_mesh.push(x_eta);
        y_xi_mesh.push(y_xi);
        y_eta_mesh.push(y_eta);
    });

    // add contribution from non-owned connected points
    mesh.boundaries.iter().for_each(|boundary| {
        match boundary {
            BlockBoundary::Connection(ranges) => {
                let mut data: Vec<Scalar> = vec![];
                let BlockConnection(own_range, other_range) = ranges;

                {
                    let x_xi = &x_xi_mesh[other_range.block];

                    assert!(other_range.end > other_range.start);

                    data.reserve(other_range.end - other_range.start);

                    let (range_i, range_j) = match other_range.edge {
                        EdgeIndex::IMin => (0..1, other_range.start..other_range.end + 1),
                        EdgeIndex::IMax => (
                            x_xi.dim().0 - 1..x_xi.dim().0,
                            other_range.start..other_range.end + 1,
                        ),
                        EdgeIndex::JMin => (other_range.start..other_range.end + 1, 0..1),
                        EdgeIndex::JMax => (
                            other_range.start..other_range.end + 1,
                            x_xi.dim().1 - 1..x_xi.dim().1,
                        ),
                    };

                    for j in range_j {
                        for i in range_i.clone() {
                            data.push(x_xi[[i, j]]);
                        }
                    }
                }

                data.reverse();

                {
                    let x_xi = &mut x_xi_mesh[own_range.block];

                    assert!(own_range.end > own_range.start);

                    assert!(own_range.end - own_range.start + 1 == data.len());

                    let (range_i, range_j) = match own_range.edge {
                        EdgeIndex::IMin => (0..1, own_range.start..own_range.end + 1),
                        EdgeIndex::IMax => (
                            x_xi.dim().0 - 1..x_xi.dim().0,
                            own_range.start..own_range.end + 1,
                        ),
                        EdgeIndex::JMin => (own_range.start..own_range.end + 1, 0..1),
                        EdgeIndex::JMax => (
                            own_range.start..own_range.end + 1,
                            x_xi.dim().1 - 1..x_xi.dim().1,
                        ),
                    };

                    for j in range_j {
                        for i in range_i.clone() {
                            x_xi[[i, j]] += data.pop().unwrap();
                        }
                    }
                }
            }
            // BlockBoundary::PeriodicConnection { connection, translation },
            _ => (),
        }
    });

    // copy data to non-owned connected points
}
