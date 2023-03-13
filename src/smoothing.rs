// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    types::{edge_view_mut, BlockBoundary, BlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use float_cmp::approx_eq;
use ndarray::{Array1, Array2};
use russell_lab::{vector_norm, NormVec, Vector};
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::error::Error;

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

pub fn smooth_block(block: &mut Block2d, iterations: usize) -> Result<(), Box<dyn Error>> {
    let coords = &mut block.coords;

    let config = ConfigSolver::new();

    let shape = coords.dim();
    let (i_size, j_size) = shape;

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

// fn compute_derivatives(mesh: &Mesh) {
//     // allocate derivative fields
//     let n = mesh.blocks.len();

//     let mut x_xi_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
//     let mut x_eta_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
//     let mut y_xi_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);
//     let mut y_eta_mesh: Vec<Array2<Scalar>> = Vec::with_capacity(n);

//     mesh.blocks.iter().for_each(|block| {
//         let shape = block.coords.dim();
//         let x = &block.coords;

//         let mut x_xi = Array2::<Scalar>::zeros(shape);
//         let mut x_eta = Array2::<Scalar>::zeros(shape);
//         let mut y_xi = Array2::<Scalar>::zeros(shape);
//         let mut y_eta = Array2::<Scalar>::zeros(shape);

//         // loop internal points
//         for j in 1..shape.1 - 1 {
//             for i in 1..shape.0 - 1 {
//                 let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//                 let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//                 let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
//                 let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//                 x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
//                 x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
//                 y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
//                 y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
//             }
//         }

//         // add own contributions on corners
//         {
//             let j = 0;
//             let i = 0;

//             let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//             let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

//             x_xi[[i, j]] = 0.5 * x_ip1_j;
//             x_eta[[i, j]] = 0.5 * x_i_jp1;
//             y_xi[[i, j]] = 0.5 * y_ip1_j;
//             y_eta[[i, j]] = 0.5 * y_i_jp1;
//         }

//         {
//             let j = shape.1 - 1;
//             let i = shape.0 - 1;

//             let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//             let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//             x_xi[[i, j]] = -0.5 * x_im1_j;
//             x_eta[[i, j]] = -0.5 * x_i_jm1;
//             y_xi[[i, j]] = -0.5 * y_im1_j;
//             y_eta[[i, j]] = -0.5 * y_i_jm1;
//         }

//         {
//             let j = shape.1 - 1;
//             let i = 0;

//             let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//             let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//             x_xi[[i, j]] = 0.5 * x_ip1_j;
//             x_eta[[i, j]] = -0.5 * x_i_jm1;
//             y_xi[[i, j]] = 0.5 * y_ip1_j;
//             y_eta[[i, j]] = -0.5 * y_i_jm1;
//         }

//         {
//             let j = 0;
//             let i = shape.0 - 1;

//             let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//             let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

//             x_xi[[i, j]] = -0.5 * x_im1_j;
//             x_eta[[i, j]] = 0.5 * x_i_jp1;
//             y_xi[[i, j]] = -0.5 * y_im1_j;
//             y_eta[[i, j]] = 0.5 * y_i_jp1;
//         }

//         // add own contributions on edges excluding corners

//         // edge i min & i max
//         for j in 1..shape.1 - 1 {
//             {
//                 // edge i min
//                 let i = 0;

//                 let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//                 let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
//                 let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//                 x_xi[[i, j]] = 0.5 * x_ip1_j;
//                 x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
//                 y_xi[[i, j]] = 0.5 * y_ip1_j;
//                 y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
//             }

//             {
//                 // edge i max
//                 let i = shape.0 - 1;

//                 let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//                 let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];
//                 let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//                 x_xi[[i, j]] = -0.5 * x_im1_j;
//                 x_eta[[i, j]] = 0.5 * (x_i_jp1 - x_i_jm1);
//                 y_xi[[i, j]] = -0.5 * y_im1_j;
//                 y_eta[[i, j]] = 0.5 * (y_i_jp1 - y_i_jm1);
//             }
//         }

//         // edge j min & j max
//         for i in 1..shape.0 - 1 {
//             {
//                 // edge j min
//                 let j = 0;

//                 let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//                 let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//                 let Vec2d(x_i_jp1, y_i_jp1) = x[[i, j + 1]];

//                 x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
//                 x_eta[[i, j]] = 0.5 * x_i_jp1;
//                 y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
//                 y_eta[[i, j]] = 0.5 * y_i_jp1;
//             }

//             {
//                 // edge j max
//                 let j = shape.1 - 1;

//                 let Vec2d(x_ip1_j, y_ip1_j) = x[[i + 1, j]];
//                 let Vec2d(x_im1_j, y_im1_j) = x[[i - 1, j]];
//                 let Vec2d(x_i_jm1, y_i_jm1) = x[[i, j - 1]];

//                 x_xi[[i, j]] = 0.5 * (x_ip1_j - x_im1_j);
//                 x_eta[[i, j]] = -0.5 * x_i_jm1;
//                 y_xi[[i, j]] = 0.5 * (y_ip1_j - y_im1_j);
//                 y_eta[[i, j]] = -0.5 * y_i_jm1;
//             }
//         }

//         x_xi_mesh.push(x_xi);
//         x_eta_mesh.push(x_eta);
//         y_xi_mesh.push(y_xi);
//         y_eta_mesh.push(y_eta);
//     });

//     // add contribution from non-owned connected points
//     mesh.edges.iter().for_each(|boundary| {
//         match boundary {
//             BlockBoundary::Connection(ranges) => {
//                 add_non_owned(&mut x_xi_mesh, ranges);
//                 add_non_owned(&mut x_eta_mesh, ranges);
//                 add_non_owned(&mut y_xi_mesh, ranges);
//                 add_non_owned(&mut y_eta_mesh, ranges);
//             }
//             // BlockBoundary::PeriodicConnection { connection, translation },
//             _ => (),
//         }
//     });

//     // copy data to non-owned connected points
//     mesh.edges.iter().for_each(|boundary| {
//         match boundary {
//             BlockBoundary::Connection(ranges) => {
//                 copy_to_non_owned(&mut x_xi_mesh, ranges);
//                 copy_to_non_owned(&mut x_eta_mesh, ranges);
//                 copy_to_non_owned(&mut y_xi_mesh, ranges);
//                 copy_to_non_owned(&mut y_eta_mesh, ranges);
//             }
//             // BlockBoundary::PeriodicConnection { connection, translation },
//             _ => (),
//         }
//     });
// }

// fn add_non_owned(fields: &mut Vec<Array2<Scalar>>, ranges: &BlockConnection) {
//     let BlockConnection(own_range, other_range) = ranges;

//     let edge_view_other = crate::types::edge_view(&fields, other_range);

//     let data: Vec<Scalar> = edge_view_other.iter().cloned().collect();

//     let mut edge_view_own = crate::types::edge_view_mut(fields, own_range);

//     edge_view_own
//         .iter_mut()
//         .zip(data.iter())
//         .for_each(|(target, source)| *target += *source);
// }

// fn copy_to_non_owned(fields: &mut Vec<Array2<Scalar>>, ranges: &BlockConnection) {
//     let BlockConnection(own_range, other_range) = ranges;

//     let edge_view_own = crate::types::edge_view(fields, own_range);
//     let data: Vec<Scalar> = edge_view_own.iter().cloned().collect();

//     let mut edge_view_other = crate::types::edge_view_mut(fields, other_range);

//     edge_view_other
//         .iter_mut()
//         .zip(data.iter())
//         .for_each(|(target, source)| *target = *source);
// }

#[derive(Debug, Clone)]
enum PointProps {
    Solve { stencil: PointStencilData },
    Fix,
    Connect { donor_index: usize }, // SolvePeriodic,
}

// impl PointProps {
//     fn matrix_index(&self) -> usize {
//         match self {
//             PointProps::Solve { stencil } => stencil.at(PointStencilIndex::Center).index,
//             PointProps::Fix { matrix_index } => *matrix_index,
//             PointProps::Connect { matrix_index, .. } => *matrix_index,
//         }
//     }
// }

struct MatrixEntry {
    index: usize,
    prop: PointProps,
}

impl MatrixEntry {
    fn new(index: usize, prop: PointProps) -> Self {
        Self { index, prop }
    }
}

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    // assign matrix indices
    let n = mesh.blocks.len();

    // sets the matrix properties for every point. Initialize with fixed points
    // on the boundary and solved for (moved) inner points
    let mut points_to_matrix: Vec<Array2<MatrixEntry>> = Vec::with_capacity(n);

    let start_indices: Vec<usize> = mesh
        .blocks
        .iter()
        .map(|block| {
            let (size_i, size_j) = block.coords.dim();
            (size_i - 1) * (size_j - 1)
        })
        .scan(0 as usize, |state, internal_size| {
            let idx = *state;
            *state += internal_size;
            Some(idx)
        })
        .collect();

    mesh.blocks
        .iter()
        .enumerate()
        .zip(start_indices.iter())
        .for_each(|((block_idx, block), start_idx)| {
            let coords = &block.coords;

            // matrix indices of all block points
            let point_indices = Array1::<usize>::from_shape_fn(coords.len(), |i| i + start_idx)
                .into_shape(coords.dim())
                .unwrap();

            // initialize point to matrix mapping to fixed points
            let mut block_points_to_matrix =
                Array2::<MatrixEntry>::from_shape_fn(coords.raw_dim(), |(i, j)| {
                    MatrixEntry::new(point_indices[[i, j]], PointProps::Fix)
                });

            // set internal points to be solved
            block_points_to_matrix
                .slice_mut(ndarray::s![1..coords.dim().0 - 1, 1..coords.dim().1 - 1])
                .indexed_iter_mut()
                .for_each(|((i, j), p)| {
                    // transform array view index to block index
                    let i = i + 1;
                    let j = j + 1;
                    let stencil = PointStencilData::new([
                        PointIndex::new(block_idx, i, j),
                        PointIndex::new(block_idx, i + 1, j),
                        PointIndex::new(block_idx, i - 1, j),
                        PointIndex::new(block_idx, i, j + 1),
                        PointIndex::new(block_idx, i, j - 1),
                        PointIndex::new(block_idx, i + 1, j + 1),
                        PointIndex::new(block_idx, i + 1, j - 1),
                        PointIndex::new(block_idx, i - 1, j + 1),
                        PointIndex::new(block_idx, i - 1, j - 1),
                    ]);
                    p.prop = PointProps::Solve { stencil };
                });

            points_to_matrix.push(block_points_to_matrix);
        });

    // add the edge connections (set donor range to be solved and receiver range
    // to connect)

    mesh.edges.iter().for_each(|edge| {
        if let BlockBoundary::Connection(connection) = edge {
            // set the donor indices to be solved
            connection.donor.iter().for_each(|(i, j)| {
                // assemble point stencil containing parts in other block to be
                // changed
                let block_idx = connection.donor.block;

                let mut stencil = PointStencilData::new([
                    PointIndex::new(block_idx, i, j),         // 0
                    PointIndex::new(block_idx, i + 1, j),     // 1
                    PointIndex::new(block_idx, i - 1, j),     // 2
                    PointIndex::new(block_idx, i, j + 1),     // 3
                    PointIndex::new(block_idx, i, j - 1),     // 4
                    PointIndex::new(block_idx, i + 1, j + 1), // 5
                    PointIndex::new(block_idx, i + 1, j - 1), // 6
                    PointIndex::new(block_idx, i - 1, j + 1), // 7
                    PointIndex::new(block_idx, i - 1, j - 1), // 8
                ]);

                let rec_block_idx = connection.receiver.block;

                if connection.donor.start[[0, 0]] == connection.donor.end[[0, 0]] {
                    if connection.donor.start[[0, 0]] == 0 {
                        // j_min edge
                        let rec_idx =
                            connection.get_index_in_receiver_block(&[i as isize - 1, j as isize]);
                        stencil.stencil[2] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize - 1, j as isize + 1]);
                        stencil.stencil[7] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize - 1, j as isize - 1]);
                        stencil.stencil[8] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);
                    } else {
                        // j_max edge
                        let rec_idx =
                            connection.get_index_in_receiver_block(&[i as isize + 1, j as isize]);
                        stencil.stencil[1] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize + 1, j as isize + 1]);
                        stencil.stencil[5] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize + 1, j as isize - 1]);
                        stencil.stencil[6] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);
                    }
                } else {
                    if connection.donor.start[[1, 0]] == 0 {
                        // i_min edge
                        let rec_idx =
                            connection.get_index_in_receiver_block(&[i as isize, j as isize - 1]);
                        stencil.stencil[4] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize + 1, j as isize - 1]);
                        stencil.stencil[6] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize - 1, j as isize - 1]);
                        stencil.stencil[8] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);
                    } else {
                        // i_max edge
                        let rec_idx =
                            connection.get_index_in_receiver_block(&[i as isize, j as isize + 1]);
                        stencil.stencil[3] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize + 1, j as isize + 1]);
                        stencil.stencil[5] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);

                        let rec_idx = connection
                            .get_index_in_receiver_block(&[i as isize - 1, j as isize + 1]);
                        stencil.stencil[7] =
                            PointIndex::new(rec_block_idx, rec_idx.0 as usize, rec_idx.1 as usize);
                    }
                }

                points_to_matrix[connection.donor.block][[i, j]].prop =
                    PointProps::Solve { stencil };
            });

            // // remove end points of the connection if needed
            // {
            //     let (i, j) = connection.receiver.first();
            //     // let mut stencil = &mut points_to_matrix[connection.donor.block][[i, j]].prop;
            //     if let PointProps::Solve { stencil } =
            //         &points_to_matrix[connection.donor.block][[i, j]].prop
            //     {
            //         if connection.donor.start[[0, 0]] == connection.donor.end[[0, 0]] {
            //             if connection.donor.start[[0, 0]] == 0 {
            //             } else {
            //             }
            //         } else {
            //             if connection.donor.start[[1, 0]] == 0 {
            //             } else {
            //             }
            //         }
            //     };
            // }
            // {
            //     let (i, j) = connection.receiver.last();
            // }

            // set the receiver indices to be connected to the donor indices
            connection
                .receiver
                .iter()
                .zip(connection.donor.iter())
                .for_each(|(rec_idx, donor_idx)| {
                    let donor_index =
                        points_to_matrix[connection.donor.block][[donor_idx.0, donor_idx.1]].index;
                    points_to_matrix[connection.receiver.block][[rec_idx.0, rec_idx.1]].prop =
                        PointProps::Connect { donor_index };
                });
        };
    });

    // reinforce fixed data to precede over connections on the outer points of
    // the edge, which can be part of two segments due to overlap

    mesh.edges.iter().for_each(|edge| {
        let mut set_fixed_closure = |range| {
            let mut points = edge_view_mut(&mut points_to_matrix, range);

            points
                .iter_mut()
                .for_each(|point| point.prop = PointProps::Fix {});
        };

        match edge {
            BlockBoundary::Connection(_) => (),
            BlockBoundary::PeriodicConnection { connection, .. } => {
                set_fixed_closure(&connection.0);
                set_fixed_closure(&connection.1);
            }
            BlockBoundary::Inlet(range) => set_fixed_closure(range),
            BlockBoundary::Outlet(range) => set_fixed_closure(range),
            BlockBoundary::Wall(range) => set_fixed_closure(range),
        };
    });

    let mut matrix_indices: Vec<Array2<Option<usize>>> = Vec::with_capacity(n);
    let mut matrix_connectivity: Vec<PointStencilData> = vec![];
    let mut matrix_fixed: Vec<PointIndex> = vec![];

    // TODO reserve matrix_fixed

    // compute matrix indices

    // assign internal points

    let start_indices: Vec<usize> = mesh
        .blocks
        .iter()
        .map(|block| {
            let (size_i, size_j) = block.coords.dim();
            (size_i - 2) * (size_j - 2)
        })
        .scan(0 as usize, |state, internal_size| {
            let idx = *state;
            *state += internal_size;
            Some(idx)
        })
        .collect();

    mesh.blocks
        .iter()
        .zip(start_indices.iter())
        .for_each(|(block, start_idx)| {
            let coords = &block.coords;
            let mut block_indices = Array2::<Option<usize>>::from_elem(coords.dim(), None);

            // add internal points to be solved
            block_indices
                .slice_mut(ndarray::s![1..coords.dim().0 - 1, 1..coords.dim().1 - 1])
                .iter_mut()
                .enumerate()
                .for_each(|(counter, index)| *index = Some(start_idx + counter));

            matrix_indices.push(block_indices);
        });

    // assign edge points

    // mesh.edges.iter().for_each(|bound| match bound {
    //     BlockBoundary::Connection() => (),
    //     BlockBoundary::PeriodicConnection {
    //         connection,
    //         translation,
    //     } => (),
    //     BlockBoundary::Inlet() => (),
    //     BlockBoundary::Outlet() => (),
    //     BlockBoundary::Wall() => (),
    // });

    // mesh.blocks.iter().for_each(|block|{
    //     block.
    // });

    // println!("{:?}", matrix_indices);

    // compute_derivatives(mesh);

    // compute DOF

    let dof = mesh.points();

    // let dof = matrix_indices
    //     .iter()
    //     .map(|block| block.iter().filter(|i| i.is_some()).count())
    //     .sum();

    // // compute connectivity

    // matrix_connectivity.reserve(mesh.points());

    // matrix_indices.iter().enumerate().for_each(|(b, block)| {
    //     block.indexed_iter().for_each(|((i, j), index)| {
    //         if let Some(idx) = index {
    //             matrix_connectivity.push(PointStencilData::new([
    //                 PointIndex::new(*idx, b, i, j),
    //                 PointIndex::new(*idx, b, i + 1, j),
    //                 PointIndex::new(*idx, b, i - 1, j),
    //                 PointIndex::new(*idx, b, i, j + 1),
    //                 PointIndex::new(*idx, b, i, j - 1),
    //                 PointIndex::new(*idx, b, i + 1, j + 1),
    //                 PointIndex::new(*idx, b, i + 1, j - 1),
    //                 PointIndex::new(*idx, b, i - 1, j + 1),
    //                 PointIndex::new(*idx, b, i - 1, j - 1),
    //             ]));
    //         }
    //     })
    // });

    // // assemble LHS and RHS

    // let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;

    // let mut rhs_x = Vector::new(dof);
    // let mut rhs_y = Vector::new(dof);

    // points_to_matrix
    //     .iter()
    //     .enumerate()
    //     .for_each(|(block_index, block_points)| {
    //         block_points.indexed_iter().for_each(|((i, j), point)| {
    //             match point {
    //                 PointProps::Solve {
    //                     matrix_index,
    //                     stencil,
    //                 } => {
    //                     //     let i_j = point.at(PointConnectivityIndex::i_j);
    //                     //     let ip1_j = point.at(PointConnectivityIndex::ip1_j);
    //                     //     let im1_j = point.at(PointConnectivityIndex::im1_j);
    //                     //     let i_jp1 = point.at(PointConnectivityIndex::i_jp1);
    //                     //     let i_jm1 = point.at(PointConnectivityIndex::i_jm1);
    //                     //     let ip1_jp1 = point.at(PointConnectivityIndex::ip1_jp1);
    //                     //     let ip1_jm1 = point.at(PointConnectivityIndex::ip1_jm1);
    //                     //     let im1_jp1 = point.at(PointConnectivityIndex::im1_jp1);
    //                     //     let im1_jm1 = point.at(PointConnectivityIndex::im1_jm1);

    //                     //     lhs.put(idx, idx, a_i_j)?;
    //                 }
    //                 PointProps::Fix { matrix_index } => {
    //                     let Vec2d(x, y) = mesh.blocks[block_index].coords[[i, j]];
    //                     lhs.put(*matrix_index, *matrix_index, 1.0).unwrap();
    //                     rhs_x[*matrix_index] = x;
    //                     rhs_x[*matrix_index] = y;
    //                 }
    //                 PointProps::Connect {
    //                     matrix_index,
    //                     donor_index,
    //                 } => {
    //                     lhs.put(*matrix_index, *matrix_index, 1.0);
    //                     lhs.put(*donor_index, *donor_index, 1.0);
    //                 }
    //             }
    //         });
    //     });

    // // matrix_connectivity.iter().for_each(|point| {
    // //     let i_j = point.at(PointConnectivityIndex::i_j);
    // //     let ip1_j = point.at(PointConnectivityIndex::ip1_j);
    // //     let im1_j = point.at(PointConnectivityIndex::im1_j);
    // //     let i_jp1 = point.at(PointConnectivityIndex::i_jp1);
    // //     let i_jm1 = point.at(PointConnectivityIndex::i_jm1);
    // //     let ip1_jp1 = point.at(PointConnectivityIndex::ip1_jp1);
    // //     let ip1_jm1 = point.at(PointConnectivityIndex::ip1_jm1);
    // //     let im1_jp1 = point.at(PointConnectivityIndex::im1_jp1);
    // //     let im1_jm1 = point.at(PointConnectivityIndex::im1_jm1);

    // //     lhs.put(idx, idx, a_i_j)?;
    // // });

    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct PointIndex {
    pub block: usize,
    pub i: usize,
    pub j: usize,
}

impl PointIndex {
    fn new(block: usize, i: usize, j: usize) -> Self {
        Self { block, i, j }
    }
}

#[derive(Debug, Clone, Copy)]
enum PointStencilIndex {
    Center,
    East,
    West,
    North,
    South,
    NorthEast,
    SouthEast,
    SouthWest,
    NorthWest,
}

impl PointStencilIndex {
    pub fn index(&self) -> usize {
        *self as usize
    }

    pub fn ip1_j() -> Self {
        PointStencilIndex::East
    }

    pub fn im1_j() -> Self {
        PointStencilIndex::West
    }

    pub fn i_jp1() -> Self {
        PointStencilIndex::North
    }

    pub fn i_jm1() -> Self {
        PointStencilIndex::South
    }

    pub fn ip1_jp1() -> Self {
        PointStencilIndex::NorthEast
    }

    pub fn ip1_jm1() -> Self {
        PointStencilIndex::SouthEast
    }

    pub fn im1_jp1() -> Self {
        PointStencilIndex::SouthWest
    }

    pub fn im1_jm1() -> Self {
        PointStencilIndex::NorthWest
    }
}

#[derive(Debug, Clone)]
struct PointStencilData {
    stencil: [PointIndex; 9],
}

impl PointStencilData {
    pub fn new(stencil: [PointIndex; 9]) -> Self {
        Self { stencil }
    }

    pub fn at(&self, position: PointStencilIndex) -> &PointIndex {
        &self.stencil[position.index()]
    }
}
