// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    types::{edge_view_mut, BlockBoundary, BlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use float_cmp::approx_eq;
use ndarray::{s, Array1, Array2};
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
    Solve,
    Fix,
    Connect { donor_index: usize }, // SolvePeriodic,
}

#[derive(Debug, Clone)]
struct MatrixEntry {
    index: usize,
    prop: PointProps,
}

impl MatrixEntry {
    fn new(index: usize, prop: PointProps) -> Self {
        Self { index, prop }
    }
}

fn coordinates_with_ghost_points(mesh: &Mesh) -> Vec<Array2<Vec2d>> {
    let mut coords = Vec::<Array2<Vec2d>>::with_capacity(mesh.blocks.len());

    // allocation and internal data
    mesh.blocks.iter().for_each(|block| {
        let dim = block.points();
        // let mut block_coords = Array2::<Vec2d>::zeros([dim[0] + 2, dim[1] + 2]);
        let mut block_coords =
            Array2::<Vec2d>::from_elem([dim[0] + 2, dim[1] + 2], Vec2d(Scalar::NAN, Scalar::NAN));

        // copy over to internal
        block
            .coords
            .iter()
            .zip(
                block_coords
                    .slice_mut(s![1..dim[0] + 1, 1..dim[1] + 1])
                    .iter_mut(),
            )
            .for_each(|(x_ref, x)| {
                *x = *x_ref;
            });

        coords.push(block_coords);
    });

    // fill ghost points
    mesh.edges.iter().for_each(|edge| {
        if let BlockBoundary::Connection(connection) = edge {
            let ghost_point_modifyer = connection
                .donor
                .get_ghost_layer_modifyer(&mesh.blocks[connection.donor.block].points());

            connection.donor.iter().for_each(|(i, j)| {
                let i = i as isize + ghost_point_modifyer.0;
                let j = j as isize + ghost_point_modifyer.1;

                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                let x =
                    mesh.blocks[connection.receiver.block].coords[[i_rec as usize, j_rec as usize]];

                coords[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]] = x;
            });
        }
    });

    coords
}

fn matrix_entries(mesh: &Mesh) -> Vec<Array2<MatrixEntry>> {
    let mut matrix_entries: Vec<Array2<MatrixEntry>> = Vec::with_capacity(mesh.blocks.len());

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

    // internal indices
    mesh.blocks
        .iter()
        .zip(start_indices.iter())
        .for_each(|(block, start_idx)| {
            let dim = block.points();

            let mut block_entries = Array2::<MatrixEntry>::from_elem(
                [dim[0] + 2, dim[1] + 2],
                MatrixEntry::new(usize::MAX, PointProps::Fix),
            );

            block_entries
                .slice_mut(s![1..dim[0] + 1, 1..dim[1] + 1])
                .iter_mut()
                .enumerate()
                .for_each(|(counter, e)| {
                    *e = MatrixEntry::new(start_idx + counter, PointProps::Solve);
                });

            matrix_entries.push(block_entries);
        });

    // fill ghost point indices
    mesh.edges.iter().for_each(|edge| {
        if let BlockBoundary::Connection(connection) = edge {
            println!("{:?} {:?}", connection.donor.start, connection.donor.end);
            println!(
                "{:?} {:?}",
                connection.receiver.start, connection.receiver.end
            );

            let ghost_point_modifyer = connection
                .donor
                .get_ghost_layer_modifyer(&mesh.blocks[connection.donor.block].points());

            connection.donor.iter().for_each(|(i, j)| {
                let i = i as isize + ghost_point_modifyer.0;
                let j = j as isize + ghost_point_modifyer.1;

                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                let index_rec = matrix_entries[connection.receiver.block]
                    [[(i_rec + 1) as usize, (j_rec + 1) as usize]]
                .index;

                matrix_entries[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]]
                    .index = index_rec;
            });

            connection.donor.iter().for_each(|(i, j)| {
                let (i_rec, j_rec) =
                    connection.get_index_in_receiver_block(&[i as isize, j as isize]);

                let donor_index = matrix_entries[connection.donor.block]
                    [[(i + 1) as usize, (j + 1) as usize]]
                .index;

                matrix_entries[connection.receiver.block]
                    [[(i_rec + 1) as usize, (j_rec + 1) as usize]]
                .prop = PointProps::Connect { donor_index };
            });
        }
    });

    // set all other edges to be fixed
    mesh.edges.iter().for_each(|edge| {
        let mut set_fixed_closure = |range| {
            let mut points = edge_view_mut(&mut matrix_entries, range);

            points
                .iter_mut()
                .for_each(|point| point.prop = PointProps::Fix);
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

    matrix_entries
}

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    let coords = coordinates_with_ghost_points(&mesh);
    let entries = matrix_entries(&mesh);

    let dof = mesh.points();

    // assemble LHS and RHS

    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;

    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    entries
        .iter()
        .enumerate()
        .for_each(|(block_index, block_points)| {
            let dim = mesh.blocks[block_index].points();

            println!("{} {:?}", block_index, mesh.blocks[block_index].points());

            block_points
                .slice(s![1..dim[0] - 1, 1..dim[1] - 1])
                .indexed_iter()
                .for_each(|((i, j), point)| {
                    // convert to index w/ ghost points
                    let i = i + 1;
                    let j = j + 1;

                    let index = point.index;

                    match point.prop {
                        PointProps::Solve => {
                            // laplace conditions (zero intialization)
                            let s = 0.0;
                            let t = 0.0;

                            // let Vec2d(x_i_j, y_i_j) = coords[[i, j]];
                            let Vec2d(x_ip1_j, y_ip1_j) = coords[block_index][[i + 1, j]];
                            let Vec2d(x_im1_j, y_im1_j) = coords[block_index][[i - 1, j]];
                            let Vec2d(x_i_jp1, y_i_jp1) = coords[block_index][[i, j + 1]];
                            let Vec2d(x_i_jm1, y_i_jm1) = coords[block_index][[i, j - 1]];

                            // let Vec2d(x_ip1_jp1, y_ip1_jp1) = coords[block_index][[i + 1, j + 1]];
                            // let Vec2d(x_im1_jm1, y_im1_jm1) = coords[block_index][[i - 1, j - 1]];
                            // let Vec2d(x_ip1_jm1, y_ip1_jm1) = coords[block_index][[i + 1, j - 1]];
                            // let Vec2d(x_im1_jp1, y_im1_jp1) = coords[block_index][[i - 1, j + 1]];

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

                            println!("{:?}", entries[block_index]);
                            println!(
                                "{} {} : {}",
                                i - 1,
                                j,
                                entries[block_index][[i - 1, j]].index
                            );

                            lhs.put(index, index, a_i_j).unwrap();
                            lhs.put(index, entries[block_index][[i - 1, j]].index, a_im1_j)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i + 1, j]].index, a_ip1_j)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i, j - 1]].index, a_i_jm1)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i, j + 1]].index, a_i_jp1)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i + 1, j + 1]].index, a_ip1_jp1)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i + 1, j - 1]].index, a_ip1_jm1)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i - 1, j + 1]].index, a_im1_jp1)
                                .unwrap();
                            lhs.put(index, entries[block_index][[i - 1, j - 1]].index, a_im1_jm1)
                                .unwrap();
                        }
                        PointProps::Fix => {
                            let Vec2d(x, y) = mesh.blocks[block_index].coords[[i, j]];
                            lhs.put(index, index, 1.0).unwrap();
                            rhs_x[index] = x;
                            rhs_y[index] = y;
                        }
                        PointProps::Connect { donor_index } => {
                            lhs.put(index, index, 1.0).unwrap();
                            lhs.put(donor_index, donor_index, 1.0).unwrap();
                        }
                    }
                });
        });

    // // solve
    // let mut solver = Solver::new(config)?;
    // solver.initialize(&lhs)?;
    // solver.factorize()?;
    // solver.solve(&mut x_new, &rhs_x)?;
    // solver.solve(&mut y_new, &rhs_y)?;

    // let x_norm = vector_norm(&x_new, NormVec::Euc);
    // let y_norm = vector_norm(&y_new, NormVec::Euc);

    // let norm = (x_norm * x_norm + y_norm * y_norm).sqrt();

    // println!("\tresidual {norm}");

    Ok(())
}
