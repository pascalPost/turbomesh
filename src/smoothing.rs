// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

pub mod block_boundary_props;

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    smoothing::block_boundary_props::{block_boundary_points_solver_props, ConnectionData},
    types::{BlockBoundary, BlockConnection, PeriodicBlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use ndarray::{s, Array2};
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::error::Error;

use self::block_boundary_props::{BlockBoundaryPointSolverProp, BlockPointIndex};

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

    coords
}

fn fill_ghost_points(
    mesh: &Mesh,
    coords: &mut Vec<ndarray::Array2<Vec2d>>,
    // _corners: &Vec<(MeshPoint, MeshPoint)>,
) {
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
                    x + *translation;
            });
        }
    });

    // // fill corners
    // corners.iter().for_each(|(target, source)| {
    //     let x = coords[source.block][[source.point.0, source.point.1]];
    //     coords[target.block][[target.point.0, target.point.1]] = x;
    // });
}

fn try_set_ghost_point_data(
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    mesh: &Mesh,
    block_idx: usize,
    ghost_point: &(usize, usize),
    connection: &BlockConnection,
) -> bool {
    let potential_donor = connection
        .get_index_in_receiver_block(&[ghost_point.0 as isize - 1, ghost_point.1 as isize - 1]);

    // check if potential donor is a
    // valid point
    let rec_dim = mesh.blocks[connection.receiver.block].points();
    if (potential_donor.0 >= 0 || potential_donor.0 < rec_dim[0] as isize)
        && (potential_donor.1 >= 0 || potential_donor.1 < rec_dim[1] as isize)
    {
        let index = matrix_entries[connection.receiver.block][[
            potential_donor.0 as usize + 1,
            potential_donor.1 as usize + 1,
        ]]
        .index;
        matrix_entries[block_idx][[ghost_point.0, ghost_point.1]].index = index;

        true
    } else {
        false
    }
}

/// returns fields containing the unique matrix/vector indices for every point,
/// property how to treat every point (fixed, solved for, connected) and target
/// indices for every needed ghost point
fn matrix_entries(mesh: &Mesh) -> Vec<Array2<MatrixEntry>> {
    // TODO create mechanism to work for multi point conenctions (4 block
    // meeting point)
    //
    // 3 | 0
    // -- --
    // 2 | 1
    //

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

    let mut matrix_entries_new: Vec<Array2<MatrixEntry>> = Vec::with_capacity(mesh.blocks.len());

    // intialization
    mesh.blocks
        .iter()
        .zip(start_indices.iter())
        .for_each(|(block, start_idx)| {
            let dim = block.points();

            // initialize all block points including ghost  to fixed and usize::MAX as matrix index
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

            matrix_entries_new.push(block_entries);
        });

    // boundary indices
    let boundary_props = block_boundary_points_solver_props(mesh);

    boundary_props
        .iter()
        .enumerate()
        .for_each(|(block_idx, block_boundary_props)| {
            block_boundary_props
                .iter()
                .enumerate()
                .for_each(|(point_idx, point_prop)| {
                    let point = block_boundary_props.point_index(point_idx).unwrap();

                    // ghost point index
                    let i = point.0 + 1;
                    let j = point.1 + 1;

                    match point_prop {
                        BlockBoundaryPointSolverProp::Undefined => {
                            panic!("undefined block boundary point encountered")
                        }
                        BlockBoundaryPointSolverProp::Fixed => {
                            matrix_entries_new[block_idx][[i, j]].prop = PointProps::Fix
                        }
                        BlockBoundaryPointSolverProp::Solved(connections) => {
                            matrix_entries_new[block_idx][[i, j]].prop = PointProps::Solve;

                            // fill ghost points

                            // collect all needed ghost points
                            let dim = block_boundary_props.dim;
                            let ghost_points = if i == 1 {
                                [(0, j - 1), (0, j), (0, j + 1)]
                            } else if j == 1 {
                                [(i - 1, 0), (i, 0), (i + 1, 0)]
                            } else if i == dim[0] {
                                [(dim[0] + 1, j - 1), (dim[0] + 1, j), (dim[0] + 1, j + 1)]
                            } else if j == dim[1] {
                                [(i - 1, dim[1] + 1), (i, dim[1] + 1), (i + 1, dim[1] + 1)]
                            } else {
                                panic!("ghost point index out of bounds")
                            };

                            // search for ghost points donor in all connections to this point
                            for ghost_point in ghost_points.iter() {
                                let mut found = false;

                                for connection in connections.iter() {
                                    let edge = &mesh.edges[connection.edge];

                                    if connection.donor {
                                        if let BlockBoundary::Connection(connection)
                                        | BlockBoundary::PeriodicConnection(
                                            PeriodicBlockConnection { connection, .. },
                                        ) = &edge
                                        {
                                            found = try_set_ghost_point_data(
                                                &mut matrix_entries_new,
                                                mesh,
                                                block_idx,
                                                ghost_point,
                                                connection,
                                            );

                                            if found {
                                                break;
                                            }
                                        }
                                    }
                                }

                                if !found {
                                    panic!("no donor found for ghost point");
                                }
                            }
                        }
                        BlockBoundaryPointSolverProp::Connected(donor) => {
                            let donor_index = matrix_entries_new[donor.index.block]
                                [[donor.index.point.0 + 1, donor.index.point.1 + 1]]
                            .index;
                            matrix_entries_new[block_idx][[i, j]].prop =
                                PointProps::Connect { donor_index }
                        }
                        BlockBoundaryPointSolverProp::ConnectedPeriodic { donor, translation } => {
                            let donor_index = matrix_entries_new[donor.index.block]
                                [[donor.index.point.0 + 1, donor.index.point.1 + 1]]
                            .index;
                            matrix_entries_new[block_idx][[i, j]].prop =
                                PointProps::ConnectPeriodic {
                                    donor_index,
                                    translation: *translation,
                                }
                        }
                    }
                });
        });

    // check if block ghost point corners need to be filled
    mesh.blocks
        .iter()
        .enumerate()
        .for_each(|(block_idx, block)| {
            let dim = block.points();

            BlockCorner::all_corners().into_iter().for_each(|corner| {
                let point_idx_internal = BlockCorner::corner_index(corner, dim);
                // transform to contain ghost layer
                let point_idx = (point_idx_internal.0 + 1, point_idx_internal.1 + 1);

                // the ghost point corner that potentionaly needs to be set
                let ghost_point_corner =
                    BlockCorner::corner_index(corner, [dim[0] + 2, dim[1] + 2]);

                // the ghost point corner needs to be set if the corner point is
                // to be solved and it is not yet set (still usize::MAX)
                if matrix_entries_new[block_idx][[point_idx.0, point_idx.1]].prop
                    == PointProps::Solve
                    && matrix_entries_new[block_idx][[ghost_point_corner.0, ghost_point_corner.1]]
                        .index
                        == usize::MAX
                {
                    // loop over all connections to the corner point
                    let mut found = false;
                    if let BlockBoundaryPointSolverProp::Solved(connections) = boundary_props
                        [block_idx]
                        .get(point_idx_internal)
                        .expect("corner to be a boundary point")
                    {
                        for connection in connections {
                            match &mesh.edges[connection.edge] {
                                BlockBoundary::Connection(connection) => {
                                    if let Some(ghost_cell_corner_donor) = find_ghost_point_corner(
                                        mesh,
                                        block_idx,
                                        connection,
                                        ghost_point_corner,
                                    ) {
                                        let index = matrix_entries_new
                                            [ghost_cell_corner_donor.block][[
                                            ghost_cell_corner_donor.point.0,
                                            ghost_cell_corner_donor.point.1,
                                        ]]
                                        .index;
                                        matrix_entries_new[block_idx]
                                            [[ghost_point_corner.0, ghost_point_corner.1]]
                                        .index = index;
                                        matrix_entries_new[block_idx]
                                            [[ghost_point_corner.0, ghost_point_corner.1]]
                                        .prop = PointProps::Connect { donor_index: index };

                                        todo!("Set the corner in the matched block as well.");

                                        found = true;
                                        break;
                                    }
                                }
                                BlockBoundary::PeriodicConnection(PeriodicBlockConnection {
                                    connection,
                                    translation,
                                }) => {
                                    if let Some(ghost_cell_corner_donor) = find_ghost_point_corner(
                                        mesh,
                                        block_idx,
                                        connection,
                                        ghost_point_corner,
                                    ) {
                                        let index = matrix_entries_new
                                            [ghost_cell_corner_donor.block][[
                                            ghost_cell_corner_donor.point.0,
                                            ghost_cell_corner_donor.point.1,
                                        ]]
                                        .index;
                                        matrix_entries_new[block_idx]
                                            [[ghost_point_corner.0, ghost_point_corner.1]]
                                        .index = index;
                                        matrix_entries_new[block_idx]
                                            [[ghost_point_corner.0, ghost_point_corner.1]]
                                        .prop = PointProps::ConnectPeriodic {
                                            donor_index: index,
                                            translation: *translation,
                                        };

                                        todo!("Set the corner in the matched block as well.");

                                        found = true;
                                        break;
                                    }
                                }
                                _ => (),
                            }
                        }

                        if !found {
                            for connection in connections {
                                let connection_2nd_level = boundary_props[connection.index.block]
                                    .get(connection.index.point)
                                    .expect("corner to be a boundary point");

                                println!("{:?}", connection_2nd_level);

                                // if the point is not identical to the corner
                                // point
                                let conn_to_check = match connection_2nd_level {
                                    BlockBoundaryPointSolverProp::Undefined => todo!(),
                                    BlockBoundaryPointSolverProp::Fixed => todo!(),
                                    BlockBoundaryPointSolverProp::Solved(_) => todo!(),
                                    BlockBoundaryPointSolverProp::Connected(conn) => conn,
                                    BlockBoundaryPointSolverProp::ConnectedPeriodic {
                                        donor,
                                        ..
                                    } => donor,
                                };

                                match &mesh.edges[conn_to_check.edge] {
                                    BlockBoundary::Connection(connection) => {
                                        if let Some(ghost_cell_corner_donor) =
                                            find_ghost_point_corner(
                                                mesh,
                                                block_idx,
                                                connection,
                                                ghost_point_corner,
                                            )
                                        {
                                            let index = matrix_entries_new
                                                [ghost_cell_corner_donor.block][[
                                                ghost_cell_corner_donor.point.0,
                                                ghost_cell_corner_donor.point.1,
                                            ]]
                                            .index;
                                            matrix_entries_new[block_idx]
                                                [[ghost_point_corner.0, ghost_point_corner.1]]
                                            .index = index;
                                            matrix_entries_new[block_idx]
                                                [[ghost_point_corner.0, ghost_point_corner.1]]
                                            .prop = PointProps::Connect { donor_index: index };

                                            todo!("Set the corner in the matched block as well.");

                                            found = true;
                                            break;
                                        }
                                    }
                                    BlockBoundary::PeriodicConnection(
                                        PeriodicBlockConnection {
                                            connection,
                                            translation,
                                        },
                                    ) => {
                                        if let Some(ghost_cell_corner_donor) =
                                            find_ghost_point_corner(
                                                mesh,
                                                block_idx,
                                                connection,
                                                ghost_point_corner,
                                            )
                                        {
                                            let index = matrix_entries_new
                                                [ghost_cell_corner_donor.block][[
                                                ghost_cell_corner_donor.point.0,
                                                ghost_cell_corner_donor.point.1,
                                            ]]
                                            .index;
                                            matrix_entries_new[block_idx]
                                                [[ghost_point_corner.0, ghost_point_corner.1]]
                                            .index = index;
                                            matrix_entries_new[block_idx]
                                                [[ghost_point_corner.0, ghost_point_corner.1]]
                                            .prop = PointProps::ConnectPeriodic {
                                                donor_index: index,
                                                translation: *translation,
                                            };

                                            todo!("Set the corner in the matched block as well.");

                                            found = true;
                                            break;
                                        }
                                    }
                                    _ => (),
                                }
                            }
                        }

                        if !found {
                            todo!("check 3rd level connection to corner point");
                        }
                    } else {
                        panic!("corner point not solved");
                    }

                    // check both neighbor points to connect to the block that
                    // the ghost cell corner is in
                    // for neighbor in BlockCorner::corner_neighbors(corner, dim).into_iter() {
                    //     let connection = boundary_props[block_idx]
                    //         .get(neighbor)
                    //         .expect("neighbor to be a boundary point");

                    //     let mut connections = match connection {
                    //         BlockBoundaryPointSolverProp::Undefined => &[],
                    //         BlockBoundaryPointSolverProp::Fixed => &[],
                    //         BlockBoundaryPointSolverProp::Solved(connections) => {
                    //             connections.as_slice()
                    //         }
                    //         BlockBoundaryPointSolverProp::Connected(connection) => {
                    //             std::slice::from_ref(connection)
                    //         }
                    //         BlockBoundaryPointSolverProp::ConnectedPeriodic { donor, .. } => {
                    //             std::slice::from_ref(donor)
                    //         }
                    //     }
                    //     .to_vec();

                    //     // add 2nd level connections (connections of connections to neighbors)
                    //     let mut connections_2nd_level: Vec<_> = connections
                    //         .iter()
                    //         .map(|connection: &ConnectionData| -> Vec<ConnectionData> {
                    //             let conns = boundary_props[connection.index.block]
                    //                 .get(connection.index.point)
                    //                 .expect("connection point to neighbor to be a boundary point");

                    //             match conns {
                    //                 BlockBoundaryPointSolverProp::Undefined => vec![],
                    //                 BlockBoundaryPointSolverProp::Fixed => vec![],
                    //                 BlockBoundaryPointSolverProp::Solved(connections) => {
                    //                     connections.clone()
                    //                 }
                    //                 BlockBoundaryPointSolverProp::Connected(connection) => {
                    //                     vec![connection.clone()]
                    //                 }
                    //                 BlockBoundaryPointSolverProp::ConnectedPeriodic {
                    //                     donor,
                    //                     ..
                    //                 } => vec![donor.clone()],
                    //             }
                    //         })
                    //         .collect();

                    //     connections_2nd_level
                    //         .iter_mut()
                    //         .for_each(|conn| connections.append(conn));

                    //     println!(
                    //         "block {} corner {:?} neighbor {:?} connections: {:?}",
                    //         block_idx, corner, neighbor, connections
                    //     );

                    //     for connection in connections {
                    //         match &mesh.edges[connection.edge] {
                    //             BlockBoundary::Connection(connection) => {
                    //                 if let Some(ghost_cell_corner_donor) = find_ghost_point_corner(
                    //                     mesh,
                    //                     block_idx,
                    //                     connection,
                    //                     ghost_point_corner,
                    //                 ) {
                    //                     let index = matrix_entries_new
                    //                         [ghost_cell_corner_donor.block][[
                    //                         ghost_cell_corner_donor.point.0,
                    //                         ghost_cell_corner_donor.point.1,
                    //                     ]]
                    //                     .index;
                    //                     matrix_entries_new[block_idx]
                    //                         [[ghost_point_corner.0, ghost_point_corner.1]]
                    //                     .index = index;
                    //                     matrix_entries_new[block_idx]
                    //                         [[ghost_point_corner.0, ghost_point_corner.1]]
                    //                     .prop = PointProps::Connect { donor_index: index };

                    //                     break;
                    //                 }
                    //             }
                    //             BlockBoundary::PeriodicConnection(PeriodicBlockConnection {
                    //                 connection,
                    //                 translation,
                    //             }) => {
                    //                 if let Some(ghost_cell_corner_donor) = find_ghost_point_corner(
                    //                     mesh,
                    //                     block_idx,
                    //                     connection,
                    //                     ghost_point_corner,
                    //                 ) {
                    //                     let index = matrix_entries_new
                    //                         [ghost_cell_corner_donor.block][[
                    //                         ghost_cell_corner_donor.point.0,
                    //                         ghost_cell_corner_donor.point.1,
                    //                     ]]
                    //                     .index;
                    //                     matrix_entries_new[block_idx]
                    //                         [[ghost_point_corner.0, ghost_point_corner.1]]
                    //                     .index = index;
                    //                     matrix_entries_new[block_idx]
                    //                         [[ghost_point_corner.0, ghost_point_corner.1]]
                    //                     .prop = PointProps::ConnectPeriodic {
                    //                         donor_index: index,
                    //                         translation: *translation,
                    //                     };

                    //                     break;
                    //                 }
                    //             }
                    //             _ => (),
                    //         }
                    //     }
                    // }

                    assert!(
                        matrix_entries_new[block_idx][[ghost_point_corner.0, ghost_point_corner.1]]
                            .index
                            != usize::MAX,
                        "ghost cell corner not found in block {} ghost point {:?}",
                        block_idx,
                        ghost_point_corner
                    );
                }
            });
        });

    println!("1 (20, 10) {:?}", boundary_props[1].get((20, 10)));
    println!("2 ( 0, 10) {:?}", boundary_props[2].get((0, 10)));
    println!("6 (20,  0) {:?}", boundary_props[6].get((20, 0)));

    // write matrix_entries to file
    {
        use std::io::Write;
        let mut file = std::fs::File::create("matrix_entries_new.txt").unwrap();
        for (block_idx, block) in matrix_entries_new.iter().enumerate() {
            writeln!(file, "block {}", block_idx).unwrap();
            for ((i, j), entry) in block.indexed_iter() {
                writeln!(file, "{}, ({}, {}) = {:?}", block_idx, i, j, entry).unwrap();
            }
        }
    }

    matrix_entries_new
}

fn find_ghost_point_corner(
    mesh: &Mesh,
    block_idx: usize,
    connection: &BlockConnection,
    ghost_cell_corner: (usize, usize),
) -> Option<BlockPointIndex> {
    // check of donor_point exists and is thus the target
    // point_data

    let owner_block = if connection.donor.block == block_idx {
        connection.receiver.block
    } else {
        connection.donor.block
    };

    let potential_donor_point = connection
        .get_index_in_receiver_block(&[ghost_cell_corner.0 as isize, ghost_cell_corner.1 as isize]);

    let owner_dim = mesh.blocks[owner_block].points();

    if potential_donor_point.0 >= 0
        && potential_donor_point.0 < owner_dim[0] as isize
        && potential_donor_point.1 >= 0
        && potential_donor_point.1 < owner_dim[1] as isize
    {
        Some(BlockPointIndex::new(
            owner_block,
            (
                potential_donor_point.0 as usize,
                potential_donor_point.1 as usize,
            ),
        ))
    } else {
        None
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum BlockCorner {
    BottomLeft,
    BottomRight,
    TopLeft,
    TopRight,
}

impl BlockCorner {
    /// returns all possible block corners for iteration
    fn all_corners() -> [Self; 4] {
        [
            BlockCorner::BottomLeft,
            BlockCorner::BottomRight,
            BlockCorner::TopLeft,
            BlockCorner::TopRight,
        ]
    }

    /// returns the index of the corner given the block dimensions
    fn corner_index(corner: BlockCorner, dim: [usize; 2]) -> (usize, usize) {
        match corner {
            BlockCorner::BottomLeft => (0, 0),
            BlockCorner::BottomRight => (dim[0] - 1, 0),
            BlockCorner::TopLeft => (0, dim[1] - 1),
            BlockCorner::TopRight => (dim[0] - 1, dim[1] - 1),
        }
    }

    /// returns the two neighbor points of the corner
    fn corner_neighbors(corner: BlockCorner, dim: [usize; 2]) -> [(usize, usize); 2] {
        match corner {
            BlockCorner::BottomLeft => [(1, 0), (0, 1)],
            BlockCorner::BottomRight => [(dim[0] - 2, 0), (dim[0] - 1, 1)],
            BlockCorner::TopLeft => [(0, dim[1] - 2), (1, dim[1] - 1)],
            BlockCorner::TopRight => [(dim[0] - 1, dim[1] - 2), (dim[0] - 2, dim[1] - 1)],
        }
    }
}

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    let iterations = 20;

    // fields with ghost points
    let entries = matrix_entries(&mesh);
    let mut coords = coordinates_with_ghost_points(&mesh);

    // allocations
    let dof = mesh.points();
    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;
    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    // iterate
    for n in 0..iterations {
        print!("  iteration\t{n}");

        fill_ghost_points(mesh, &mut coords /* , &corners */);

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

                                let Vec2d(x_ip1_j, y_ip1_j) = coords[block_index][[i + 1, j]];
                                let Vec2d(x_im1_j, y_im1_j) = coords[block_index][[i - 1, j]];
                                let Vec2d(x_i_jp1, y_i_jp1) = coords[block_index][[i, j + 1]];
                                let Vec2d(x_i_jm1, y_i_jm1) = coords[block_index][[i, j - 1]];

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
