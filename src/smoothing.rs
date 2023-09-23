// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

pub mod block_boundary_props;

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use crate::{
    smoothing::block_boundary_props::{
        block_boundary_points_solver_props, BlockCorner, ConnectionData,
    },
    types::{BlockBoundary, BlockConnection, PeriodicBlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use ndarray::{s, Array2};
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::error::Error;

use self::block_boundary_props::{
    BlockBoundaryArray, BlockBoundaryPointProp, BlockBoundaryPointSolverProp, BlockPointIndex,
    BoundaryProps,
};

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
    // TODO index can be removed and the global matrix index can be computed given the block sizes
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

fn update_ghost_point_locations(
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

                // TODO check if this is still correct for complicated
                // connection. Perhaps we have to save a transform matrix.
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
    translation: Option<Vec2d>,
) -> bool {
    let potential_donor = connection
        .get_index_in_receiver_block(&[ghost_point.0 as isize - 1, ghost_point.1 as isize - 1]);

    // check if potential donor is a
    // valid point
    let rec_dim = mesh.blocks[connection.receiver.block].points();
    if (potential_donor.0 >= 0 && potential_donor.0 < rec_dim[0] as isize)
        && (potential_donor.1 >= 0 && potential_donor.1 < rec_dim[1] as isize)
    {
        let index = matrix_entries[connection.receiver.block][[
            potential_donor.0 as usize + 1,
            potential_donor.1 as usize + 1,
        ]]
        .index;

        assert!(index != usize::MAX);

        matrix_entries[block_idx][[ghost_point.0, ghost_point.1]].index = index;
        matrix_entries[block_idx][[ghost_point.0, ghost_point.1]].prop = translation.map_or_else(
            || PointProps::Connect { donor_index: index },
            |trans| PointProps::ConnectPeriodic {
                donor_index: index,
                translation: trans,
            },
        );

        println!(
            " . . . . . . . Found GhostPoint donor for block {}, ghost point ({}, {}) in block {}, point ({}, {}) (excluding ghost layer), matrix index {}",
            block_idx, ghost_point.0,ghost_point.1,
            connection.receiver.block, potential_donor.0, potential_donor.1, index
        );

        true
    } else {
        false
    }
}

/// returns true if the given point is to be fixed
fn is_point_fixed(
    block_id: usize,
    boundary_point_id: usize,
    boundary_props: &BoundaryProps,
    checked_edges: &mut Vec<usize>,
) -> bool {
    for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
        // println!(" . . . . . Checking {:?}", bc);

        match bc {
            BlockBoundaryPointProp::Undefined => {
                panic!("undefined block boundary point encountered.")
            }

            BlockBoundaryPointProp::Fixed => return true,

            BlockBoundaryPointProp::Connected(data) | BlockBoundaryPointProp::Periodic(data) => {
                // check every connection that has not yet been checked
                let edge_id = data.edge;

                if !checked_edges.contains(&edge_id) {
                    let con_idx = data.index;
                    let con_block = &boundary_props.data[con_idx.block];
                    let con_boundary_point_id =
                        con_block.boundary_point_index(con_idx.point).unwrap();
                    checked_edges.push(edge_id);
                    if is_point_fixed(
                        con_idx.block,
                        con_boundary_point_id,
                        boundary_props,
                        checked_edges,
                    ) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

// /// sets the given point and all its connected points to fixed
// fn set_point_fixed(
//     block_id: usize,
//     boundary_point_id: usize,
//     boundary_props: &BoundaryProps,
//     matrix_entries: &mut Vec<Array2<MatrixEntry>>,
// ) {
//     // set point to fixed
//     let point_id = boundary_props.data[block_id]
//         .point_index(boundary_point_id)
//         .unwrap();

//     println!(
//         " . . . . Setting block {} point {} / ({}, {}) to fixed.",
//         block_id, boundary_point_id, point_id.0, point_id.1
//     );

//     debug_assert_eq!(
//         matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].prop,
//         PointProps::Fix
//     );

//     matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].prop = PointProps::Fix;

//     // check for connected points
//     for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
//         match bc {
//             BlockBoundaryPointProp::Undefined => todo!(),
//             BlockBoundaryPointProp::Fixed => {
//                 // already set
//             }
//             BlockBoundaryPointProp::Connected(data) | BlockBoundaryPointProp::Periodic(data) => {
//                 let con_idx = data.index;
//                 let con_block = &boundary_props.data[con_idx.block];
//                 let con_point_id = con_block.boundary_point_index(con_idx.point).unwrap();

//                 if data.donor {
//                     // call recursively to set donor and all connected points to fixed
//                     set_point_fixed(con_idx.block, con_point_id, boundary_props, matrix_entries);
//                 } else {
//                     // set non donor to fixed

//                     println!(
//                         " . . . . Setting block {} point {} / ({}, {}) to fixed.",
//                         block_id, boundary_point_id, point_id.0, point_id.1
//                     );

//                     debug_assert_eq!(
//                         matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
//                             .prop,
//                         PointProps::Fix
//                     );

//                     matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
//                         .prop = PointProps::Fix;
//                 }
//             }
//         }
//     }
// }

/// sets the given point to solve and returns its index
fn set_solve(
    block_id: usize,
    boundary_point_id: usize,
    boundary_props: &BoundaryProps,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
) -> usize {
    let point_id = boundary_props.data[block_id]
        .point_index(boundary_point_id)
        .unwrap();

    println!(
        " . . . . Setting block {} boundaryPoint {} / ({}, {}) to solve.",
        block_id, boundary_point_id, point_id.0, point_id.1
    );

    matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].prop = PointProps::Solve;

    println!(
        " . . . . . matrix entry {}, ({}, {}) {:?}",
        block_id,
        point_id.0 + 1,
        point_id.1 + 1,
        matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]]
    );

    matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].index
}

fn connect_points(
    block_id: usize,
    boundary_point_id: usize,
    boundary_props: &BoundaryProps,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    donor_index: usize,
    mesh: &Mesh,
) {
    // let point_id = boundary_props.data[block_id]
    //     .point_index(boundary_point_id)
    //     .unwrap();

    // check for connected points
    for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
        match bc {
            BlockBoundaryPointProp::Undefined => todo!(),
            BlockBoundaryPointProp::Fixed => {
                panic!("Fixed point must not be encountered when setting connected points")
            }
            BlockBoundaryPointProp::Connected(data) => {
                let con_idx = data.index;
                let con_block = &boundary_props.data[con_idx.block];
                let con_point_id = con_block.boundary_point_index(con_idx.point).unwrap();

                if data.donor {
                    println!(
                        " . . . . Setting block {} boundaryPoint {} / ({}, {}) to connect.",
                        con_idx.block, con_point_id, con_idx.point.0, con_idx.point.1
                    );

                    matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
                        .prop = PointProps::Connect { donor_index };

                    println!(
                        " . . . . . matrix entry {}, ({}, {}) {:?}",
                        con_idx.block,
                        con_idx.point.0 + 1,
                        con_idx.point.1 + 1,
                        matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
                    );

                    // for donor points, call recursively to set all connected point
                    connect_points(
                        con_idx.block,
                        con_point_id,
                        boundary_props,
                        matrix_entries,
                        donor_index,
                        mesh,
                    );
                }
            }
            BlockBoundaryPointProp::Periodic(data) => {
                let con_idx = data.index;
                let con_block = &boundary_props.data[con_idx.block];
                let con_point_id = con_block.boundary_point_index(con_idx.point).unwrap();

                if data.donor {
                    println!(
                        " . . . . Setting block {} point {} / ({}, {}) to ConnectPeriodic.",
                        con_idx.block, con_point_id, con_idx.point.0, con_idx.point.1
                    );

                    let translation =
                        if let BlockBoundary::PeriodicConnection(data) = &mesh.edges[data.edge] {
                            data.translation
                        } else {
                            panic!("Wrong boundary type encountered.");
                        };

                    matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
                        .prop = PointProps::ConnectPeriodic {
                        donor_index,
                        translation,
                    };

                    println!(
                        " . . . . . matrix entry {}, ({}, {}) {:?}",
                        con_idx.block,
                        con_idx.point.0 + 1,
                        con_idx.point.1 + 1,
                        matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
                    );

                    // for donor points, call recursively to set all connected point
                    connect_points(
                        con_idx.block,
                        con_point_id,
                        boundary_props,
                        matrix_entries,
                        donor_index,
                        mesh,
                    );
                }
            }
        }
    }
}

/// set the matrix entry for the given point and all its connected points
fn set_matrix_entry(
    block_id: usize,
    boundary_point_id: usize,
    boundary_props: &BoundaryProps,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    mesh: &Mesh,
) {
    // check if the matrix entry is already set to connect or connect periodic
    let point_id = boundary_props.data[block_id]
        .point_index(boundary_point_id)
        .unwrap();

    if let PointProps::Connect { .. } | PointProps::ConnectPeriodic { .. } =
        matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].prop
    {
        println!(
            " . . . . . matrix entry {}, ({}, {}) {:?}",
            block_id,
            point_id.0 + 1,
            point_id.1 + 1,
            matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]]
        );
        return;
    }

    // (1.1) check tree for fixed points
    // println!(" . . . . Checking for fixed points");
    let mut checked_edges = Vec::<usize>::with_capacity(10);
    if is_point_fixed(
        block_id,
        boundary_point_id,
        boundary_props,
        &mut checked_edges,
    ) {
        // println!(" . . . . Point is to be fixed.");
        println!(" . . . . Point remains fixed.");

        // (1.2) set all connected points to fixed

        // TODO does not need to be done as matrix entries are by default fixed
        // set_point_fixed(block_id, boundary_point_id, boundary_props,
        // matrix_entries);

        println!(
            " . . . . . matrix entry {}, ({}, {}) {:?}",
            block_id,
            point_id.0 + 1,
            point_id.1 + 1,
            matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]]
        );
    } else {
        // (2.1) set minimal connected point to solve
        let donor_index = set_solve(block_id, boundary_point_id, boundary_props, matrix_entries);

        // (2.2) set all other connected points to connected or periodic connected
        connect_points(
            block_id,
            boundary_point_id,
            boundary_props,
            matrix_entries,
            donor_index,
            mesh,
        );
    }
}

// fn check_entry(
//     entry: &BlockBoundaryPointProp,
//     block_id: usize,
//     point_id: usize,
//     block_boundary_props: &BlockBoundaryArray<Vec<BlockBoundaryPointProp>>,
//     matrix_entries: &mut Vec<Array2<MatrixEntry>>,
// ) -> Option<PointProps> {
//     match entry {
//         block_boundary_props::BlockBoundaryPointProp::Undefined => todo!(),
//         block_boundary_props::BlockBoundaryPointProp::Fixed => {
//             // set point to fixed
//             // let point_index = block_boundary_props.point_index(point_id).unwrap();
//             // matrix_entries[block_id][[point_index.0, point_index.1]].prop = PointProps::Fix;
//             // println!(
//             //     " . . . Final Matrix entry: {:?}",
//             //     matrix_entries[block_id][[point_index.0, point_index.1]].prop
//             // );
//             return Some(PointProps::Fix);
//         }
//         block_boundary_props::BlockBoundaryPointProp::Connected(connected_point) => {
//             // println!("{:?}", connected_point);

//             // only check if donor is true

//             println!(
//                 " . . . . . Connected Point BCs: {:?}",
//                 boundary_props[&connected_point.index]
//             );
//             // connected_point

//             // search all connected points.

//             // if the connected point is smaller, check
//             // for solver properties

//             todo!()
//         }
//         block_boundary_props::BlockBoundaryPointProp::Periodic(_) => {
//             todo!()
//         }
//     }
// }

struct MatrixEntries {
    /// matrix entries as vector of blocks containing the matrix entries for
    /// every 2d block point
    pub data: Vec<Array2<MatrixEntry>>,
}

impl MatrixEntries {
    fn new(mesh: &Mesh) -> MatrixEntries {
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

        let mut matrix_entries: Vec<Array2<MatrixEntry>> = Vec::with_capacity(mesh.blocks.len());

        // intialization
        mesh.blocks
            .iter()
            .zip(start_indices.iter())
            .for_each(|(block, start_idx)| {
                let dim = block.points();

                // initialize all block points including ghost to fixed and usize::MAX as matrix index
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

        // boundary points
        let boundary_props = BoundaryProps::new(mesh);

        println!("Setting block boundary point matirx properties.");

        boundary_props
            .data
            .iter()
            .enumerate()
            .for_each(|(block_id, block_boundary_props)| {
                println!(" . Block {}", block_id);

                block_boundary_props.iter().enumerate().for_each(
                    |(boundary_point_id, point_props)| {
                        println!(
                            " . . Block {} BoundaryPoint {} / {:?}",
                            block_id,
                            boundary_point_id,
                            block_boundary_props.point_index(boundary_point_id).unwrap()
                        );
                        println!(" . . . BCs: {:?}", point_props);

                        set_matrix_entry(
                            block_id,
                            boundary_point_id,
                            &boundary_props,
                            &mut matrix_entries,
                            mesh,
                        );
                    },
                );
            });

        // let mut checked_edges = Vec::<usize>::with_capacity(10);
        // assert!(!is_point_fixed(7, 119, &boundary_props,&mut checked_edges));

        // println!("2, (40,9) / (41,10): {:?}", matrix_entries[2][[41, 10]]);
        // println!("7, (0,1) / (1,2): {:?}", matrix_entries[7][[1, 2]]);

        // println!("0, (0,0) / (1,1): {:?}", matrix_entries[0][[1, 1]]);
        // println!("0, (1,0) / (2,1): {:?}", matrix_entries[0][[2, 1]]);
        // println!("0, (1,-1) / (2,0): {:?}", matrix_entries[0][[2, 0]]);
        // println!("");
        // println!("0, (10,0) / (11,1): {:?}", matrix_entries[0][[11, 1]]);
        // println!("0, (10,1) / (11,2): {:?}", matrix_entries[0][[11, 2]]);
        // println!("1, (0,1) / (1,2): {:?}", matrix_entries[1][[1, 2]]);
        // println!("");
        // println!("0, (10,9) / (11,10): {:?}", matrix_entries[0][[11, 10]]);
        // println!("0, (10,10) / (11,11): {:?}", matrix_entries[0][[11, 11]]);

        println!("5, (41,1) / (42,2): {:?}", matrix_entries[5][[42, 2]]);
        println!("6, (19,30) / (20,31): {:?}", matrix_entries[6][[20, 31]]);
        println!("6, (19,31) / (20,32): {:?}", matrix_entries[6][[20, 32]]);

        // std::process::exit(0);

        println!("Searching ghost points for solution points.");

        boundary_props
            .data
            .iter()
            .enumerate()
            .for_each(|(block_id, block_boundary_props)| {
                println!(" . Block {}", block_id);

                block_boundary_props.iter().enumerate().for_each(
                    |(boundary_point_id, point_props)| {
                        let (i, j) = block_boundary_props.point_index(boundary_point_id).unwrap();

                        println!(
                            " . . Block {} BoundaryPoint {} / ({}, {}) or ({}, {}) including ghost point layers",
                            block_id, boundary_point_id, i, j, i + 1, j + 1
                        );

                        println!(
                            " . . . Prop: {:?}",
                            matrix_entries[block_id][[i + 1, j + 1]]
                        );
                        println!(" . . . BCs: {:?}", point_props);

                        if let PointProps::Solve = matrix_entries[block_id][[i + 1, j + 1]].prop {
                            fill_ghost_points(mesh, &boundary_props, block_id, boundary_point_id, &mut matrix_entries);
                        }
                    },
                );
            });

        MatrixEntries {
            data: matrix_entries,
        }
    }
}

/// function to recursively find the donor of the given ghost point
// fn search_ghost_point_donor() {}

/// fill the ghost points for the given boundary point in ghost layer including
/// point index
fn fill_ghost_points(
    mesh: &Mesh,
    boundary_props: &BoundaryProps,
    block_id: usize,
    boundary_point_id: usize,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
) {
    let point_id = boundary_props.data[block_id]
        .point_index(boundary_point_id)
        .unwrap();

    // add ghost point layer to index
    let (i, j) = (point_id.0 + 1, point_id.1 + 1);

    // collect all needed ghost points
    let dim = mesh.blocks[block_id].point_size;
    let ghost_points = if i == 1 {
        [(0, j - 1), (0, j), (0, j + 1)]
    } else if j == 1 {
        [(i - 1, 0), (i, 0), (i + 1, 0)]
    } else if i == dim.0 {
        [(dim.0 + 1, j - 1), (dim.0 + 1, j), (dim.0 + 1, j + 1)]
    } else if j == dim.1 {
        [(i - 1, dim.1 + 1), (i, dim.1 + 1), (i + 1, dim.1 + 1)]
    } else {
        panic!("ghost point index out of bounds")
    };

    println!(
        " . . . Needed GhostPoints (index space including ghost point layes): {:?}",
        ghost_points
    );

    // search for ghost points donor in all connections to this point
    for ghost_point in ghost_points.iter() {
        println!(" . . . . Searching GhostPoint {:?}", ghost_point);

        let mut found = false;

        for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
            println!(" . . . . . Checking BC {:?}", bc);

            match bc {
                BlockBoundaryPointProp::Undefined => todo!(),
                BlockBoundaryPointProp::Fixed => todo!(),
                BlockBoundaryPointProp::Connected(connection)
                | BlockBoundaryPointProp::Periodic(connection) => {
                    println!(" . . . . . . ConnectionData {:?}", connection);
                    // println!(" . . . . . . EdgeData {:?}",
                    // mesh.edges[connection.edge]);

                    if !connection.donor {
                        println!(" . . . . . . Skipping as not a donor connection");
                        continue;
                    }

                    let edge_id = connection.edge;
                    let edge = &mesh.edges[edge_id];

                    match edge {
                        BlockBoundary::Connection(connection) => {
                            println!(
                                " . . . . . . Checking connection to block {} via edge {}",
                                connection.receiver.block, edge_id
                            );

                            assert!(connection.donor.block == block_id);
                            // TODO assert if point is in donor range

                            if try_set_ghost_point_data(
                                matrix_entries,
                                mesh,
                                block_id,
                                ghost_point,
                                connection,
                                None,
                            ) {
                                found = true;
                                break;
                            }
                        }
                        BlockBoundary::PeriodicConnection(PeriodicBlockConnection {
                            connection,
                            translation,
                        }) => {
                            println!(
                                " . . . . . . Checking periodic connection to block {} via edge {}",
                                connection.receiver.block, edge_id
                            );

                            assert!(connection.donor.block == block_id);
                            // TODO assert if point is in donor range

                            if try_set_ghost_point_data(
                                matrix_entries,
                                mesh,
                                block_id,
                                ghost_point,
                                connection,
                                Some(translation.clone()),
                            ) {
                                found = true;
                                break;
                            }
                        }
                        BlockBoundary::Inlet(_) => todo!(),
                        BlockBoundary::Outlet(_) => todo!(),
                        BlockBoundary::Wall(_) => todo!(),
                    }
                }
            }
        }

        if !found {
            panic!("no donor found for ghost point");
        }
    }
}

/// returns fields containing the unique matrix/vector indices for every point,
/// property how to treat every point (fixed, solved for, connected) and target
/// indices for every needed ghost point
// fn matrix_entries(mesh: &Mesh) -> Vec<Array2<MatrixEntry>> {
//     // TODO create mechanism to work for multi point conenctions (4 block
//     // meeting point)
//     //
//     // 3 | 0
//     // -- --
//     // 2 | 1
//     //

//     let start_indices: Vec<usize> = mesh
//         .blocks
//         .iter()
//         .map(|block| block.coords.len())
//         .scan(0 as usize, |state, internal_size| {
//             let idx = *state;
//             *state += internal_size;
//             Some(idx)
//         })
//         .collect();

//     let mut matrix_entries_new: Vec<Array2<MatrixEntry>> = Vec::with_capacity(mesh.blocks.len());

//     // intialization
//     mesh.blocks
//         .iter()
//         .zip(start_indices.iter())
//         .for_each(|(block, start_idx)| {
//             let dim = block.points();

//             // initialize all block points including ghost  to fixed and usize::MAX as matrix index
//             let mut block_entries = Array2::<MatrixEntry>::from_elem(
//                 [dim[0] + 2, dim[1] + 2],
//                 MatrixEntry::new(usize::MAX, PointProps::Fix),
//             );

//             // set index for all block points
//             block_entries
//                 .slice_mut(s![1..dim[0] + 1, 1..dim[1] + 1])
//                 .iter_mut()
//                 .enumerate()
//                 .for_each(|(counter, e)| {
//                     e.index = start_idx + counter;
//                 });

//             // set all internal block points to solve
//             block_entries
//                 .slice_mut(s![2..dim[0], 2..dim[1]])
//                 .iter_mut()
//                 .for_each(|e| {
//                     e.prop = PointProps::Solve;
//                 });

//             matrix_entries_new.push(block_entries);
//         });

//     // boundary indices
//     let boundary_props = block_boundary_points_solver_props(mesh);

//     // println!("boundary props: {:#?}", boundary_props);

//     // boundary_props
//     //     .iter()
//     //     .enumerate()
//     //     .for_each(|(block_idx, block)| {
//     //         println!("block: {block_idx}");
//     //         block.iter().enumerate().for_each(|(point_idx, point)| {
//     //             println!(" . point: {point_idx}");
//     //             println!(" . . index: {:?}", block.point_index(point_idx).unwrap());
//     //             println!(" . . data: {:#?}", point);
//     //         });
//     //         // println!(" . props: {:#?}", block);
//     //     });

//     // exit(0);

//     boundary_props
//         .iter()
//         .enumerate()
//         .for_each(|(block_idx, block_boundary_props)| {
//             block_boundary_props
//                 .iter()
//                 .enumerate()
//                 .for_each(|(point_idx, point_prop)| {
//                     let point = block_boundary_props.point_index(point_idx).unwrap();

//                     // ghost point index
//                     let i = point.0 + 1;
//                     let j = point.1 + 1;

//                     match point_prop {
//                         BlockBoundaryPointSolverProp::Undefined => {
//                             panic!("undefined block boundary point encountered")
//                         }
//                         BlockBoundaryPointSolverProp::Fixed => {
//                             matrix_entries_new[block_idx][[i, j]].prop = PointProps::Fix
//                         }
//                         BlockBoundaryPointSolverProp::Solved(connections) => {
//                             matrix_entries_new[block_idx][[i, j]].prop = PointProps::Solve;

//                             // fill ghost points

//                             // collect all needed ghost points
//                             let dim = block_boundary_props.dim;
//                             let ghost_points = if i == 1 {
//                                 [(0, j - 1), (0, j), (0, j + 1)]
//                             } else if j == 1 {
//                                 [(i - 1, 0), (i, 0), (i + 1, 0)]
//                             } else if i == dim[0] {
//                                 [(dim[0] + 1, j - 1), (dim[0] + 1, j), (dim[0] + 1, j + 1)]
//                             } else if j == dim[1] {
//                                 [(i - 1, dim[1] + 1), (i, dim[1] + 1), (i + 1, dim[1] + 1)]
//                             } else {
//                                 panic!("ghost point index out of bounds")
//                             };

//                             // search for ghost points donor in all connections to this point
//                             for ghost_point in ghost_points.iter() {
//                                 let mut found = false;

//                                 for connection in connections.iter() {
//                                     let edge = &mesh.edges[connection.edge];

//                                     if connection.donor {
//                                         if let BlockBoundary::Connection(connection)
//                                         | BlockBoundary::PeriodicConnection(
//                                             PeriodicBlockConnection { connection, .. },
//                                         ) = &edge
//                                         {
//                                             found = try_set_ghost_point_data(
//                                                 &mut matrix_entries_new,
//                                                 mesh,
//                                                 block_idx,
//                                                 ghost_point,
//                                                 connection,
//                                             );

//                                             if found {
//                                                 break;
//                                             }
//                                         }
//                                     }
//                                 }

//                                 if !found {
//                                     panic!("no donor found for ghost point");
//                                 }
//                             }
//                         }
//                         BlockBoundaryPointSolverProp::Connected(donor) => {
//                             let donor_index = matrix_entries_new[donor.index.block]
//                                 [[donor.index.point.0 + 1, donor.index.point.1 + 1]]
//                             .index;
//                             matrix_entries_new[block_idx][[i, j]].prop =
//                                 PointProps::Connect { donor_index }
//                         }
//                         BlockBoundaryPointSolverProp::ConnectedPeriodic { donor, translation } => {
//                             let donor_index = matrix_entries_new[donor.index.block]
//                                 [[donor.index.point.0 + 1, donor.index.point.1 + 1]]
//                             .index;
//                             matrix_entries_new[block_idx][[i, j]].prop =
//                                 PointProps::ConnectPeriodic {
//                                     donor_index,
//                                     translation: *translation,
//                                 }
//                         }
//                     }
//                 });
//         });

//     println!(
//         "1 (20, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[1]
//             .get((20, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[1][[21, 11]]
//     );

//     println!(
//         "5 (40, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[5]
//             .get((40, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[5][[41, 11]]
//     );

//     println!(
//         "2 (0, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[2]
//             .get((0, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[2][[1, 11]]
//     );

//     println!(
//         "6 (20, 40) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[6]
//             .get((20, 40))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[6][[21, 41]]
//     );

//     println!("\n\n\n");

//     // check if block ghost point corners need to be filled
//     mesh.blocks
//         .iter()
//         .enumerate()
//         .for_each(|(block_idx, block)| {
//             let dim = block.points();

//             BlockCorner::all_corners().into_iter().for_each(|corner| {
//                 let point_idx_internal = BlockCorner::corner_index(corner, dim);
//                 // transform to contain ghost layer
//                 let point_idx = (point_idx_internal.0 + 1, point_idx_internal.1 + 1);

//                 // the ghost point corner that potentionaly needs to be set
//                 let ghost_point_corner =
//                     BlockCorner::corner_index(corner, [dim[0] + 2, dim[1] + 2]);

//                 // the ghost point corner needs to be set if the corner point is
//                 // to be solved
//                 if matrix_entries_new[block_idx][[point_idx.0, point_idx.1]].prop
//                     == PointProps::Solve
//                 {
//                     // loop over all connections to the corner point
//                     let mut found = false;
//                     if let BlockBoundaryPointSolverProp::Solved(corner_connections) = boundary_props
//                         [block_idx]
//                         .get(point_idx_internal)
//                         .expect("corner to be a boundary point")
//                     {
//                         // check matrix properties of the corner connection
//                         // (must be connect or connect periodic)
//                         for corner_connection in corner_connections {
//                             match matrix_entries_new[corner_connection.index.block][[
//                                 corner_connection.index.point.0 + 1,
//                                 corner_connection.index.point.1 + 1,
//                             ]]
//                             .prop
//                             {
//                                 PointProps::Connect { .. } => (),
//                                 PointProps::ConnectPeriodic { .. } => (),
//                                 _ => {
//                                     panic!("corner connection must be connect or connect periodic")
//                                 }
//                             }

//                             match &mesh.edges[corner_connection.edge] {
//                                 BlockBoundary::Connection(corner_connection_data) => {
//                                     if find_ghost_point_corner(
//                                         mesh,
//                                         block_idx,
//                                         corner_connection_data,
//                                         ghost_point_corner,
//                                         &mut matrix_entries_new,
//                                         None,
//                                         &boundary_props,
//                                         &corner_connection.index,
//                                         point_idx,
//                                     ) {
//                                         found = true;
//                                         break;
//                                     }
//                                 }
//                                 BlockBoundary::PeriodicConnection(
//                                     periodic_corner_connection_data,
//                                 ) => {
//                                     if find_ghost_point_corner(
//                                         mesh,
//                                         block_idx,
//                                         &periodic_corner_connection_data.connection,
//                                         ghost_point_corner,
//                                         &mut matrix_entries_new,
//                                         Some(periodic_corner_connection_data.translation),
//                                         &boundary_props,
//                                         &corner_connection.index,
//                                         point_idx,
//                                     ) {
//                                         found = true;
//                                         break;
//                                     }
//                                 }
//                                 _ => (),
//                             }
//                         }

//                         // if !found {
//                         //     // search the connections to the connections (2nd level)

//                         //     for corner_connection in corner_connections {
//                         //         // match &mesh.edges[corner_connection.edge] {
//                         //         //     BlockBoundary::Connection(corner_connection_data) => {
//                         //         //         if find_ghost_point_corner(
//                         //         //             mesh,
//                         //         //             block_idx,
//                         //         //             corner_connection_data,
//                         //         //             ghost_point_corner,
//                         //         //             &mut matrix_entries_new,
//                         //         //             None,
//                         //         //         ) {
//                         //         //             found = true;
//                         //         //             break;
//                         //         //         }
//                         //         //     }
//                         //         //     BlockBoundary::PeriodicConnection(
//                         //         //         periodic_corner_connection_data,
//                         //         //     ) => {
//                         //         //         if find_ghost_point_corner(
//                         //         //             mesh,
//                         //         //             block_idx,
//                         //         //             &periodic_corner_connection_data.connection,
//                         //         //             ghost_point_corner,
//                         //         //             &mut matrix_entries_new,
//                         //         //             Some(periodic_corner_connection_data.translation),
//                         //         //         ) {
//                         //         //             found = true;
//                         //         //             break;
//                         //         //         }
//                         //         //     }
//                         //         //     _ => (),
//                         //         // }

//                         //         let corner_connection_connection = boundary_props
//                         //             [corner_connection.index.block]
//                         //             .get(corner_connection.index.point)
//                         //             .expect("corner to be a boundary point");

//                         //         // if the point is not identical to the corner
//                         //         // point
//                         //         let conn_to_check = match corner_connection_connection {
//                         //             BlockBoundaryPointSolverProp::Undefined => todo!(),
//                         //             BlockBoundaryPointSolverProp::Fixed => todo!(),
//                         //             BlockBoundaryPointSolverProp::Solved(_) => todo!(),
//                         //             BlockBoundaryPointSolverProp::Connected(conn) => conn,
//                         //             BlockBoundaryPointSolverProp::ConnectedPeriodic {
//                         //                 donor,
//                         //                 ..
//                         //             } => donor,
//                         //         };

//                         //         let donor_block = if corner_connection.index.block == block_idx {
//                         //             connection.receiver.block
//                         //         } else {
//                         //             connection.donor.block
//                         //         };

//                         //         let potential_donor_point =
//                         //             connection.get_index_in_receiver_block(&[
//                         //                 ghost_point_corner.0 as isize,
//                         //                 ghost_point_corner.1 as isize,
//                         //             ]);

//                         //         match &mesh.edges[conn_to_check.edge] {
//                         //             BlockBoundary::Connection(connection) => {
//                         //                 if let Some(ghost_cell_corner_donor) =
//                         //                     find_ghost_point_corner(
//                         //                         mesh,
//                         //                         block_idx,
//                         //                         connection,
//                         //                         ghost_point_corner,
//                         //                     )
//                         //                 {
//                         //                     let index = matrix_entries_new
//                         //                         [ghost_cell_corner_donor.block][[
//                         //                         ghost_cell_corner_donor.point.0,
//                         //                         ghost_cell_corner_donor.point.1,
//                         //                     ]]
//                         //                     .index;
//                         //                     matrix_entries_new[block_idx]
//                         //                         [[ghost_point_corner.0, ghost_point_corner.1]]
//                         //                     .index = index;
//                         //                     matrix_entries_new[block_idx]
//                         //                         [[ghost_point_corner.0, ghost_point_corner.1]]
//                         //                     .prop = PointProps::Connect { donor_index: index };

//                         //                     found = true;
//                         //                     break;
//                         //                 }
//                         //             }
//                         //             BlockBoundary::PeriodicConnection(
//                         //                 PeriodicBlockConnection {
//                         //                     connection,
//                         //                     translation,
//                         //                 },
//                         //             ) => {
//                         //                 if let Some(ghost_cell_corner_donor) =
//                         //                     find_ghost_point_corner(
//                         //                         mesh,
//                         //                         block_idx,
//                         //                         connection,
//                         //                         ghost_point_corner,
//                         //                     )
//                         //                 {
//                         //                     let index = matrix_entries_new
//                         //                         [ghost_cell_corner_donor.block][[
//                         //                         ghost_cell_corner_donor.point.0,
//                         //                         ghost_cell_corner_donor.point.1,
//                         //                     ]]
//                         //                     .index;
//                         //                     matrix_entries_new[block_idx]
//                         //                         [[ghost_point_corner.0, ghost_point_corner.1]]
//                         //                     .index = index;
//                         //                     matrix_entries_new[block_idx]
//                         //                         [[ghost_point_corner.0, ghost_point_corner.1]]
//                         //                     .prop = PointProps::ConnectPeriodic {
//                         //                         donor_index: index,
//                         //                         translation: *translation,
//                         //                     };

//                         //                     found = true;
//                         //                     break;
//                         //                 }
//                         //             }
//                         //             _ => (),
//                         //         }
//                         //     }
//                         // }
//                     } else {
//                         panic!("corner point not solved");
//                     }

//                     // assert!(
//                     //     matrix_entries_new[block_idx][[ghost_point_corner.0, ghost_point_corner.1]]
//                     //         .index
//                     //         != usize::MAX,
//                     //     "ghost cell corner not found in block {} ghost point {:?}",
//                     //     block_idx,
//                     //     ghost_point_corner
//                     // );
//                 }
//             });
//         });

//     println!(
//         "1 (20, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[1]
//             .get((20, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[1][[21, 11]]
//     );

//     println!(
//         "5 (40, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[5]
//             .get((40, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[5][[41, 11]]
//     );

//     println!(
//         "2 (0, 10) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[2]
//             .get((0, 10))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[2][[1, 11]]
//     );

//     println!(
//         "6 (20, 0) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[6]
//             .get((20, 0))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[6][[21, 41]]
//     );

//     println!(
//         "6 (20, 40) boundary_props: {:?} matrix_entry: {:?}",
//         boundary_props[6]
//             .get((20, 40))
//             .expect("corner to be a boundary point"),
//         matrix_entries_new[6][[21, 41]]
//     );

//     std::process::exit(0);

//     // write matrix_entries to file
//     {
//         use std::io::Write;
//         let mut file = std::fs::File::create("matrix_entries_new.txt").unwrap();
//         for (block_idx, block) in matrix_entries_new.iter().enumerate() {
//             writeln!(file, "block {}", block_idx).unwrap();
//             for ((i, j), entry) in block.indexed_iter() {
//                 writeln!(file, "{}, ({}, {}) = {:?}", block_idx, i, j, entry).unwrap();
//             }
//         }
//     }

//     matrix_entries_new
// }

// fn find_ghost_point_corner(
//     mesh: &Mesh,
//     block_idx: usize,
//     connection: &BlockConnection,
//     ghost_point_corner: (usize, usize),
//     matrix_entries: &mut Vec<Array2<MatrixEntry>>,
//     translation: Option<Vec2d>,
//     boundary_props: &Vec<BlockBoundaryArray<BlockBoundaryPointSolverProp>>,
//     corner_connection: &BlockPointIndex,
//     corner_index: (usize, usize),
// ) -> bool {
//     // check of donor_point exists and is thus the target
//     // point_data

//     let potential_donor_point = connection.get_index_in_receiver_block(&[
//         ghost_point_corner.0 as isize,
//         ghost_point_corner.1 as isize,
//     ]);

//     let donor_block_dim = mesh.blocks[corner_connection.block].points();

//     if potential_donor_point.0 >= 0
//         && potential_donor_point.0 < donor_block_dim[0] as isize
//         && potential_donor_point.1 >= 0
//         && potential_donor_point.1 < donor_block_dim[1] as isize
//     {
//         let donor_point = potential_donor_point;

//         let index = matrix_entries[corner_connection.block]
//             [[donor_point.0 as usize, donor_point.1 as usize]]
//         .index;

//         matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].index = index;

//         if translation.is_none() {
//             matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].prop =
//                 PointProps::Connect { donor_index: index };
//         } else {
//             matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].prop =
//                 PointProps::ConnectPeriodic {
//                     donor_index: index,
//                     translation: translation.unwrap(),
//                 };
//         }

//         true
//     } else {
//         // search connections to connection (2nd level)
//         let corner_connection_connection = boundary_props[corner_connection.block]
//             .get(corner_connection.point)
//             .expect("corner to be a boundary point");

//         let corner_connection_connection_data = match corner_connection_connection {
//             BlockBoundaryPointSolverProp::Connected(conn) => conn,
//             BlockBoundaryPointSolverProp::ConnectedPeriodic { donor, .. } => donor,
//             _ => panic!("unexpected corner connection connection type"),
//         };

//         let data = match &mesh.edges[corner_connection_connection_data.edge] {
//             BlockBoundary::Connection(con) => (con, None),
//             BlockBoundary::PeriodicConnection(per_con) => {
//                 (&per_con.connection, Some(per_con.translation))
//             }
//             _ => {
//                 panic!("unexpected corner connection connection type")
//             }
//         };

//         let potential_donor_point = data
//             .0
//             .get_index_in_receiver_block(&[potential_donor_point.0, potential_donor_point.1]);

//         let donor_block_dim = mesh.blocks[corner_connection_connection_data.index.block].points();

//         if potential_donor_point.0 >= 0
//             && potential_donor_point.0 < donor_block_dim[0] as isize
//             && potential_donor_point.1 >= 0
//             && potential_donor_point.1 < donor_block_dim[1] as isize
//         {
//             // set target
//             let donor_point = potential_donor_point;

//             let index = matrix_entries[corner_connection_connection_data.index.block]
//                 [[donor_point.0 as usize, donor_point.1 as usize]]
//             .index;

//             matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].index = index;

//             if data.1.is_none() {
//                 matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].prop =
//                     PointProps::Connect { donor_index: index };
//             } else {
//                 matrix_entries[block_idx][[ghost_point_corner.0, ghost_point_corner.1]].prop =
//                     PointProps::ConnectPeriodic {
//                         donor_index: index,
//                         translation: data.1.unwrap(),
//                     };
//             }

//             // check donor solver settings
//             let index = matrix_entries[block_idx][[corner_index.0, corner_index.1]].index;

//             if data.1.is_none() {
//                 matrix_entries[corner_connection_connection_data.index.block]
//                     [[donor_point.0 as usize, donor_point.1 as usize]]
//                 .prop = PointProps::Connect { donor_index: index };
//             } else {
//                 matrix_entries[corner_connection_connection_data.index.block]
//                     [[donor_point.0 as usize, donor_point.1 as usize]]
//                 .prop = PointProps::ConnectPeriodic {
//                     donor_index: index,
//                     translation: data.1.unwrap(),
//                 };
//             }

//             true
//         } else {
//             false
//         }
//     }
// }

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    let iterations = 20;

    // fields with ghost points
    let entries = MatrixEntries::new(mesh).data;

    // TODO remove
    // write matrix_entries to file
    {
        use std::io::Write;
        let mut file = std::fs::File::create("matrix_entries_new_algo.txt").unwrap();
        for (block_idx, block) in entries.iter().enumerate() {
            writeln!(file, "block {}", block_idx).unwrap();
            for ((i, j), entry) in block.indexed_iter() {
                writeln!(file, "{}, ({}, {}) = {:?}", block_idx, i, j, entry).unwrap();
            }
        }
    }

    // let entries = matrix_entries(&mesh);
    let mut coords = coordinates_with_ghost_points(&mesh);

    // allocations
    let dof = mesh.points();
    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;
    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    // iterate
    for n in 0..iterations {
        print!("  iteration\t{n}");

        update_ghost_point_locations(mesh, &mut coords /* , &corners */);

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
