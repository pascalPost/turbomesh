// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    types::{BlockBoundary, BlockBoundaryRange, BlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use float_cmp::approx_eq;
use ndarray::{s, Array2};
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::error::Error;

pub enum SmoothingMethod {
    Global,
    BlockInternal,
}

fn matrix_index(i: usize, j: usize, i_size: usize) -> usize {
    (i - 1) + (j - 1) * (i_size - 2)
}

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

        let mut x_norm_sqr = 0.0;
        let mut y_norm_sqr = 0.0;

        for j in 1..shape.1 - 1 {
            for i in 1..shape.0 - 1 {
                let idx = matrix_index(i, j, i_size);
                let Vec2d(x, y) = coords[[i, j]];
                let dx = x - x_new[idx];
                let dy = y - y_new[idx];

                x_norm_sqr += dx * dx;
                y_norm_sqr += dy * dy;
            }
        }

        let norm = (x_norm_sqr + y_norm_sqr).sqrt();

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

#[derive(PartialEq, Debug, Clone)]
enum PointProps {
    Solve,
    Fix,
    Connect {
        donor_index: usize,
    },
    ConnectPeriodic {
        donor_index: usize,
        translation: Vec2d,
    },
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
            let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

            connection.donor.iter().for_each(|(i, j)| {
                let i = i as isize + ghost_point_modifyer.0;
                let j = j as isize + ghost_point_modifyer.1;

                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                let x =
                    mesh.blocks[connection.receiver.block].coords[[i_rec as usize, j_rec as usize]];

                coords[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]] = x;
            });
        }
        // TODO reduce code duplication
        else if let BlockBoundary::PeriodicConnection(periodic_connection) = edge {
            let connection = &periodic_connection.connection;
            let translation = &periodic_connection.translation;

            let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

            connection.donor.iter().for_each(|(i, j)| {
                let i = i as isize + ghost_point_modifyer.0;
                let j = j as isize + ghost_point_modifyer.1;

                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                let x =
                    mesh.blocks[connection.receiver.block].coords[[i_rec as usize, j_rec as usize]];

                coords[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]] =
                    x - *translation;
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
        .map(|block| block.coords.len())
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

            // set index for all block points
            block_entries
                .slice_mut(s![1..dim[0] + 1, 1..dim[1] + 1])
                .iter_mut()
                .enumerate()
                .for_each(|(counter, e)| {
                    e.index = start_idx + counter;
                });

            // set all internal block points to solve
            block_entries
                .slice_mut(s![2..dim[0], 2..dim[1]])
                .iter_mut()
                .for_each(|e| {
                    e.prop = PointProps::Solve;
                });

            matrix_entries.push(block_entries);
        });

    // setup connections: set ghost points
    mesh.edges.iter().for_each(|edge| {
        if let BlockBoundary::Connection(connection) = edge {
            let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

            connection.donor.iter().for_each(|(i, j)| {
                // ghost point index
                let i_gp = i as isize + ghost_point_modifyer.0;
                let j_gp = j as isize + ghost_point_modifyer.1;

                // ghost point index in receiver block (physical index)
                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i_gp, j_gp]);
                let index_rec = matrix_entries[connection.receiver.block]
                    [[(i_rec + 1) as usize, (j_rec + 1) as usize]]
                .index;

                // set ghost point to point to receiver matrix entry
                matrix_entries[connection.donor.block]
                    [[(i_gp + 1) as usize, (j_gp + 1) as usize]]
                .index = index_rec;

                // set connection point to be solved
                matrix_entries[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]].prop =
                    PointProps::Solve;
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
        // TODO remove code duplication
        else if let BlockBoundary::PeriodicConnection(periodic_connection) = edge {
            let connection = &periodic_connection.connection;
            let translation = &periodic_connection.translation;

            let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

            connection.donor.iter().for_each(|(i, j)| {
                // ghost point index
                let i_gp = i as isize + ghost_point_modifyer.0;
                let j_gp = j as isize + ghost_point_modifyer.1;

                // ghost point index in receiver block (physical index)
                let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i_gp, j_gp]);
                let index_rec = matrix_entries[connection.receiver.block]
                    [[(i_rec + 1) as usize, (j_rec + 1) as usize]]
                .index;

                // set ghost point to point to receiver matrix entry
                matrix_entries[connection.donor.block]
                    [[(i_gp + 1) as usize, (j_gp + 1) as usize]]
                .index = index_rec;

                // set connection point to be solved
                matrix_entries[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]].prop =
                    PointProps::Solve;
            });

            connection.donor.iter().for_each(|(i, j)| {
                let (i_rec, j_rec) =
                    connection.get_index_in_receiver_block(&[i as isize, j as isize]);

                let donor_index = matrix_entries[connection.donor.block]
                    [[(i + 1) as usize, (j + 1) as usize]]
                .index;

                matrix_entries[connection.receiver.block]
                    [[(i_rec + 1) as usize, (j_rec + 1) as usize]]
                .prop = PointProps::ConnectPeriodic {
                    donor_index,
                    translation: *translation,
                };
            });
        }
    });

    // set all other edges to be fixed
    mesh.edges.iter().for_each(|edge| {
        let mut set_fixed_closure = |range: &BlockBoundaryRange| {
            range.iter().for_each(|(i, j)| {
                // to ghost point index
                let i = i + 1;
                let j = j + 1;

                matrix_entries[range.block][[i, j]].prop = PointProps::Fix;
            });
        };

        match edge {
            BlockBoundary::Connection(_) => (),
            BlockBoundary::PeriodicConnection { .. } => {}
            BlockBoundary::Inlet(range) => set_fixed_closure(range),
            BlockBoundary::Outlet(range) => set_fixed_closure(range),
            BlockBoundary::Wall(range) => set_fixed_closure(range),
        };
    });

    // recheck connection start and end points if one is to be fixed
    mesh.edges.iter().for_each(|edge| {
        if let BlockBoundary::Connection(connection) = edge {
            {
                let (i, j) = connection.donor.first();
                let (i_rec, j_rec) = connection.receiver.first();

                assert!(
                    (i_rec as isize, j_rec as isize)
                        == connection.get_index_in_receiver_block(&[i as isize, j as isize])
                );

                if matrix_entries[connection.donor.block][[i + 1, j + 1]].prop == PointProps::Fix
                    || matrix_entries[connection.receiver.block][[i_rec + 1, j_rec + 1]].prop
                        == PointProps::Fix
                {
                    matrix_entries[connection.donor.block][[i + 1, j + 1]].prop = PointProps::Fix;
                    matrix_entries[connection.receiver.block][[i_rec + 1, j_rec + 1]].prop =
                        PointProps::Fix;
                }
            }

            {
                let (i, j) = connection.donor.last();
                let (i_rec, j_rec) = connection.receiver.last();

                assert!(
                    (i_rec as isize, j_rec as isize)
                        == connection.get_index_in_receiver_block(&[i as isize, j as isize])
                );

                if matrix_entries[connection.donor.block][[i + 1, j + 1]].prop == PointProps::Fix
                    || matrix_entries[connection.receiver.block][[i_rec + 1, j_rec + 1]].prop
                        == PointProps::Fix
                {
                    matrix_entries[connection.donor.block][[i + 1, j + 1]].prop = PointProps::Fix;
                    matrix_entries[connection.receiver.block][[i_rec + 1, j_rec + 1]].prop =
                        PointProps::Fix;
                }
            }
        }
    });

    matrix_entries
}

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    let iterations = 20;

    // fields with ghost points
    let mut coords = coordinates_with_ghost_points(&mesh);
    let entries = matrix_entries(&mesh);

    // allocations
    let dof = mesh.points();
    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;
    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    // iterate
    for n in 0..iterations {
        print!("  iteration\t{n}");

        lhs.reset();
        rhs_x.fill(0.0);
        rhs_y.fill(0.0);

        // assemble LHS and RHS

        entries
            .iter()
            .enumerate()
            .for_each(|(block_index, block_points)| {
                let dim = block_points.dim();

                block_points
                    .slice(s![1..dim.0 - 1, 1..dim.1 - 1])
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

                                assert!(a_i_j.is_finite());
                                assert!(a_ip1_j.is_finite());
                                assert!(a_im1_j.is_finite());
                                assert!(a_i_jp1.is_finite());
                                assert!(a_i_jm1.is_finite());
                                assert!(a_ip1_jp1.is_finite());
                                assert!(a_ip1_jm1.is_finite());
                                assert!(a_im1_jp1.is_finite());
                                assert!(a_im1_jm1.is_finite());

                                lhs.put(index, index, a_i_j).unwrap();
                                lhs.put(index, entries[block_index][[i - 1, j]].index, a_im1_j)
                                    .unwrap();
                                lhs.put(index, entries[block_index][[i + 1, j]].index, a_ip1_j)
                                    .unwrap();
                                lhs.put(index, entries[block_index][[i, j - 1]].index, a_i_jm1)
                                    .unwrap();
                                lhs.put(index, entries[block_index][[i, j + 1]].index, a_i_jp1)
                                    .unwrap();
                                lhs.put(
                                    index,
                                    entries[block_index][[i + 1, j + 1]].index,
                                    a_ip1_jp1,
                                )
                                .unwrap();
                                lhs.put(
                                    index,
                                    entries[block_index][[i + 1, j - 1]].index,
                                    a_ip1_jm1,
                                )
                                .unwrap();
                                lhs.put(
                                    index,
                                    entries[block_index][[i - 1, j + 1]].index,
                                    a_im1_jp1,
                                )
                                .unwrap();
                                lhs.put(
                                    index,
                                    entries[block_index][[i - 1, j - 1]].index,
                                    a_im1_jm1,
                                )
                                .unwrap();
                            }
                            PointProps::Fix => {
                                let Vec2d(x, y) = coords[block_index][[i, j]];
                                lhs.put(index, index, 1.0).unwrap();
                                rhs_x[index] = x;
                                rhs_y[index] = y;
                            }
                            PointProps::Connect { donor_index } => {
                                lhs.put(index, index, 1.0).unwrap();
                                lhs.put(index, donor_index, -1.0).unwrap();
                            }
                            PointProps::ConnectPeriodic {
                                donor_index,
                                translation,
                            } => {
                                lhs.put(index, index, 1.0).unwrap();
                                lhs.put(index, donor_index, -1.0).unwrap();
                                rhs_x[index] = translation.0;
                                rhs_y[index] = translation.1;
                            }
                        }
                    });
            });

        // write matrix to ascii file for debugging
        // let (m, n) = lhs.dims();
        // let mut a = Matrix::new(m, n);
        // lhs.to_matrix(&mut a)?;
        // std::fs::write("matrix.txt", format!("{}", a)).unwrap();

        // solve
        let config = ConfigSolver::new();

        let mut x_new = Vector::new(dof);
        let mut y_new = Vector::new(dof);

        let mut solver = Solver::new(config)?;

        solver.initialize(&lhs)?;
        solver.factorize()?;
        solver.solve(&mut x_new, &rhs_x)?;
        solver.solve(&mut y_new, &rhs_y)?;

        // update coords field with solution and compute norm

        let mut x_norm_sqr = 0.0;
        let mut y_norm_sqr = 0.0;

        entries
            .iter()
            .zip(coords.iter_mut())
            .for_each(|(block_entries, block_coords)| {
                assert!(block_entries.dim() == block_coords.dim());
                let dim = block_entries.dim();
                let internal = s![1..dim.0 - 1, 1..dim.1 - 1];

                block_coords
                    .slice_mut(internal)
                    .iter_mut()
                    .zip(block_entries.slice(internal).iter())
                    .for_each(|(x, info)| {
                        let i = info.index;

                        let x_new = Vec2d(x_new[i], y_new[i]);

                        let dx = x_new - *x;

                        x_norm_sqr += dx.0 * dx.0;
                        y_norm_sqr += dx.1 * dx.1;

                        *x = x_new;
                    });
            });

        let norm = (x_norm_sqr + y_norm_sqr).sqrt();

        println!("\tresidual {norm}");

        // TODO check if this is necessary (for periodic it is necessary)
        // fill ghost points
        mesh.edges.iter().for_each(|edge| {
            if let BlockBoundary::Connection(connection) = edge {
                let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

                connection.donor.iter().for_each(|(i, j)| {
                    let i = i as isize + ghost_point_modifyer.0;
                    let j = j as isize + ghost_point_modifyer.1;

                    let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                    let x = mesh.blocks[connection.receiver.block].coords
                        [[i_rec as usize, j_rec as usize]];

                    coords[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]] = x;
                });
            }
            // TODO reduce code duplication
            else if let BlockBoundary::PeriodicConnection(periodic_connection) = edge {
                let connection = &periodic_connection.connection;
                let translation = &periodic_connection.translation;
                let ghost_point_modifyer = connection.donor.get_ghost_layer_modifyer();

                connection.donor.iter().for_each(|(i, j)| {
                    let i = i as isize + ghost_point_modifyer.0;
                    let j = j as isize + ghost_point_modifyer.1;

                    let (i_rec, j_rec) = connection.get_index_in_receiver_block(&[i, j]);

                    let x = mesh.blocks[connection.receiver.block].coords
                        [[i_rec as usize, j_rec as usize]];

                    coords[connection.donor.block][[(i + 1) as usize, (j + 1) as usize]] =
                        x + *translation;
                });
            }
        });
    }

    // update mesh coordinates
    coords
        .iter()
        .zip(mesh.blocks.iter_mut())
        .for_each(|(block_coords_ref, block)| {
            let dim = block_coords_ref.dim();
            block
                .coords
                .iter_mut()
                .zip(
                    block_coords_ref
                        .slice(s![1..dim.0 - 1, 1..dim.1 - 1])
                        .iter(),
                )
                .for_each(|(x, x_ref)| {
                    *x = *x_ref;
                });
        });

    Ok(())
}
