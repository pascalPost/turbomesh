// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

pub mod block_boundary_props;

// for mumps compilation, see also https://github.com/scivision/mumps

// runs based on slightly modified repo https://github.com/cpmech/russell

use self::block_boundary_props::{BlockBoundaryPointProp, BlockPointIndex, BoundaryProps};
use crate::{
    types::{BlockBoundary, BlockConnection, PeriodicBlockConnection},
    Block2d, Mesh, Scalar, Vec2d,
};
use log::{debug, log_enabled, Level};
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
    coords: &mut Vec<ndarray::Array2<Vec2d>>,
    updates: &Vec<(BlockPointIndex, BlockPointIndex, Option<Vec2d>)>,
) {
    updates.iter().for_each(|update| {
        let target = update.0;
        let source = update.1;

        let x = coords[source.block][[source.point.0, source.point.1]];

        coords[target.block][[target.point.0, target.point.1]] =
            update.2.map_or_else(|| x, |translation| x + translation);
    });
}

/// returns true if the given point is to be fixed
fn is_point_fixed(
    block_id: usize,
    boundary_point_id: usize,
    boundary_props: &BoundaryProps,
    checked_edges: &mut Vec<usize>,
) -> bool {
    for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
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

    debug!(
        " . . . . Setting block {} boundaryPoint {} / ({}, {}) to solve.",
        block_id, boundary_point_id, point_id.0, point_id.1
    );

    matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]].prop = PointProps::Solve;

    debug!(
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
                    debug!(
                        " . . . . Setting block {} boundaryPoint {} / ({}, {}) to connect.",
                        con_idx.block, con_point_id, con_idx.point.0, con_idx.point.1
                    );

                    matrix_entries[con_idx.block][[con_idx.point.0 + 1, con_idx.point.1 + 1]]
                        .prop = PointProps::Connect { donor_index };

                    debug!(
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
                    debug!(
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

                    debug!(
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
        debug!(
            " . . . . . matrix entry {}, ({}, {}) {:?}",
            block_id,
            point_id.0 + 1,
            point_id.1 + 1,
            matrix_entries[block_id][[point_id.0 + 1, point_id.1 + 1]]
        );
        return;
    }

    // (1.1) check tree for fixed points
    let mut checked_edges = Vec::<usize>::with_capacity(10);
    if is_point_fixed(
        block_id,
        boundary_point_id,
        boundary_props,
        &mut checked_edges,
    ) {
        debug!(" . . . . Point remains fixed.");

        // (1.2) set all connected points to fixed

        debug!(
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

struct MatrixEntries {
    /// matrix entries as vector of blocks containing the matrix entries for
    /// every 2d block point
    pub data: Vec<Array2<MatrixEntry>>,

    /// number of degrees of freedom (number of mesh points + periodic ghost points)
    pub dof: usize,
}

impl MatrixEntries {
    fn new(mesh: &Mesh) -> MatrixEntries {
        let mut dof = mesh.points();

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

        debug!("Setting block boundary point matrix properties.");

        boundary_props
            .data
            .iter()
            .enumerate()
            .for_each(|(block_id, block_boundary_props)| {
                debug!(" . Block {}", block_id);

                block_boundary_props.iter().enumerate().for_each(
                    |(boundary_point_id, point_props)| {
                        debug!(
                            " . . Block {} BoundaryPoint {} / {:?}",
                            block_id,
                            boundary_point_id,
                            block_boundary_props.point_index(boundary_point_id).unwrap()
                        );
                        debug!(" . . . BCs: {:?}", point_props);

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

        debug!("Searching ghost points for solution points.");

        boundary_props
            .data
            .iter()
            .enumerate()
            .for_each(|(block_id, block_boundary_props)| {
                debug!(" . Block {}", block_id);

                block_boundary_props.iter().enumerate().for_each(
                    |(boundary_point_id, point_props)| {
                        let (i, j) = block_boundary_props.point_index(boundary_point_id).unwrap();

                        debug!(
                            " . . Block {} BoundaryPoint {} / ({}, {}) or ({}, {}) including ghost point layers",
                            block_id, boundary_point_id, i, j, i + 1, j + 1
                        );

                        debug!(
                            " . . . Prop: {:?}",
                            matrix_entries[block_id][[i + 1, j + 1]]
                        );
                        debug!(" . . . BCs: {:?}", point_props);

                        if let PointProps::Solve = matrix_entries[block_id][[i + 1, j + 1]].prop {
                            fill_ghost_points(mesh, &boundary_props, block_id, boundary_point_id, &mut matrix_entries, &mut dof);
                        }
                    },
                );
            });

        MatrixEntries {
            data: matrix_entries,
            dof,
        }
    }
}

enum GhostPointCheckResult {
    Found { index: usize },
    NotFound { potential_donor: BlockPointIndex },
}

/// if the ghost point is not found, the potential donor is returned for subsequent searches
fn check_ghost_point(
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    mesh: &Mesh,
    block_idx: usize,
    ghost_point: &(usize, usize),
    connection: &BlockConnection,
) -> GhostPointCheckResult {
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

        debug!(
            " . . . . . . . Found GhostPoint donor for block {}, ghost point ({}, {}) in block {}, point ({}, {}) (excluding ghost layer), matrix index {}",
            block_idx, ghost_point.0,ghost_point.1,
            connection.receiver.block, potential_donor.0, potential_donor.1, index
        );

        GhostPointCheckResult::Found { index }
    } else {
        GhostPointCheckResult::NotFound {
            potential_donor: BlockPointIndex {
                block: connection.receiver.block,
                point: (
                    (potential_donor.0 + 1).try_into().unwrap(),
                    (potential_donor.1 + 1).try_into().unwrap(),
                ),
            },
        }
    }
}

/// function to recursively find the donor of the given ghost point
fn recursive_search_and_set_ghost_point(
    block_id: usize,
    boundary_point_id: usize,
    ghost_point: &(usize, usize),
    mesh: &Mesh,
    boundary_props: &BoundaryProps,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    checked_edges: &mut Vec<usize>,
) -> Option<(usize, Option<Vec2d>)> {
    debug!(
        " . . . . . Checking BCs of block {} boundaryPoint {}: {:?}",
        block_id, boundary_point_id, boundary_props.data[block_id].points[boundary_point_id]
    );

    for bc in boundary_props.data[block_id].points[boundary_point_id].iter() {
        match bc {
            BlockBoundaryPointProp::Connected(connection)
            | BlockBoundaryPointProp::Periodic(connection) => {
                debug!(" . . . . . . ConnectionData {:?}", connection);

                let edge_id = connection.edge;
                let edge = &mesh.edges[edge_id];

                if !checked_edges.contains(&edge_id) {
                    match edge {
                        BlockBoundary::Connection(connection) => {
                            debug!(
                                " . . . . . . Checking connection to block {} via edge {}",
                                connection.receiver.block, edge_id
                            );

                            // the transformation does only work from donor to receiver
                            assert!(connection.donor.block == block_id);
                            // TODO assert if point is in donor range

                            match check_ghost_point(
                                matrix_entries,
                                mesh,
                                block_id,
                                ghost_point,
                                connection,
                            ) {
                                GhostPointCheckResult::Found { index } => {
                                    return Some((index, None));
                                }
                                GhostPointCheckResult::NotFound { potential_donor } => {
                                    checked_edges.push(edge_id);

                                    // convert donor boundary point to
                                    // overlapping boundary point on receiver
                                    // block

                                    let donor_bc_point = boundary_props.data[block_id]
                                        .point_index(boundary_point_id)
                                        .unwrap();

                                    let receiver_bc_point =
                                        connection.get_index_in_receiver_block(&[
                                            donor_bc_point.0 as isize,
                                            donor_bc_point.1 as isize,
                                        ]);

                                    let receiver_boundary_point_id = boundary_props.data
                                        [connection.receiver.block]
                                        .boundary_point_index((
                                            receiver_bc_point.0.try_into().unwrap(),
                                            receiver_bc_point.1.try_into().unwrap(),
                                        ))
                                        .unwrap();

                                    let res = recursive_search_and_set_ghost_point(
                                        potential_donor.block,
                                        receiver_boundary_point_id,
                                        &potential_donor.point,
                                        mesh,
                                        boundary_props,
                                        matrix_entries,
                                        checked_edges,
                                    );

                                    if res.is_some() {
                                        return res;
                                    }
                                }
                            }
                        }
                        BlockBoundary::PeriodicConnection(PeriodicBlockConnection {
                            connection,
                            translation,
                        }) => {
                            debug!(
                                " . . . . . . Checking periodic connection to block {} via edge {}",
                                connection.receiver.block, edge_id
                            );

                            match check_ghost_point(
                                matrix_entries,
                                mesh,
                                block_id,
                                ghost_point,
                                connection,
                            ) {
                                GhostPointCheckResult::Found { index } => {
                                    return Some((index, Some(translation.clone())));
                                }
                                GhostPointCheckResult::NotFound { potential_donor: _ } => {}
                            }
                        }
                        _ => todo!(),
                    }
                } else {
                    debug!(" . . . . . . Skipping already checked edge {}", edge_id);
                }
            }
            _ => todo!(),
        }
    }

    None
}

/// fill the ghost points for the given boundary point in ghost layer including
/// point index
fn fill_ghost_points(
    mesh: &Mesh,
    boundary_props: &BoundaryProps,
    block_id: usize,
    boundary_point_id: usize,
    matrix_entries: &mut Vec<Array2<MatrixEntry>>,
    dof: &mut usize,
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

    debug!(
        " . . . Needed GhostPoints (index space including ghost point layes): {:?}",
        ghost_points
    );

    // search for ghost points donor in all connections to this point
    for ghost_point in ghost_points.iter() {
        debug!(" . . . . Searching GhostPoint {:?}", ghost_point);

        // check if ghost point is already set
        if let PointProps::Connect { .. } | PointProps::ConnectPeriodic { .. } =
            matrix_entries[block_id][[ghost_point.0, ghost_point.1]].prop
        {
            debug!(
                " . . . . . ghost point block {} ({}, {}) already set: {:?}.",
                block_id,
                ghost_point.0,
                ghost_point.1,
                matrix_entries[block_id][[ghost_point.0, ghost_point.1]]
            );
            continue;
        }

        let mut checked_edges = Vec::<usize>::with_capacity(10);
        let donor = recursive_search_and_set_ghost_point(
            block_id,
            boundary_point_id,
            ghost_point,
            mesh,
            boundary_props,
            matrix_entries,
            &mut checked_edges,
        );

        if donor.is_none() {
            panic!(
                "no donor found for ghost point block {} ({}, {}) [index including ghost layers]",
                block_id, ghost_point.0, ghost_point.1
            );
        } else {
            let (index, translation) = donor.unwrap();

            if translation.is_none() {
                matrix_entries[block_id][[ghost_point.0, ghost_point.1]].index = index;
                matrix_entries[block_id][[ghost_point.0, ghost_point.1]].prop =
                    PointProps::Connect { donor_index: index };
            } else {
                // introduce new matrix index
                matrix_entries[block_id][[ghost_point.0, ghost_point.1]].index = *dof;
                matrix_entries[block_id][[ghost_point.0, ghost_point.1]].prop =
                    PointProps::ConnectPeriodic {
                        donor_index: index,
                        translation: -translation.unwrap(),
                    };

                *dof += 1;
            }
        }
    }
}

struct SingleIndexConverter {
    block_sizes: Vec<(usize, usize)>,
}

impl SingleIndexConverter {
    fn new(mesh: &Mesh) -> Self {
        let mut block_sizes = Vec::<(usize, usize)>::with_capacity(mesh.blocks.len());

        mesh.blocks.iter().for_each(|block| {
            block_sizes.push(block.point_size);
        });

        Self { block_sizes }
    }

    /// returns the block point index for the given index
    fn get_block_point_index(&self, index: usize) -> BlockPointIndex {
        let mut buf = 0;

        for (block_id, block_dim) in self.block_sizes.iter().enumerate() {
            let block_size = block_dim.0 * block_dim.1;
            if index < buf + block_size {
                let block_index = index - buf;
                let j: usize = block_index % block_dim.1;
                let i: usize = block_index / block_dim.1;
                return BlockPointIndex::new(block_id, (i + 1, j + 1));
            }
            buf += block_size;
        }

        panic!("Could not convert the given index.");
    }
}

fn compute_ghost_point_updates(
    mesh: &Mesh,
    matrix_entries: &Vec<Array2<MatrixEntry>>,
) -> Vec<(BlockPointIndex, BlockPointIndex, Option<Vec2d>)> {
    let max_size = matrix_entries
        .iter()
        .map(|block| {
            let dim = block.dim();
            2 * dim.0 + 2 * (dim.1 - 2)
        })
        .sum();

    let mut updates =
        Vec::<(BlockPointIndex, BlockPointIndex, Option<Vec2d>)>::with_capacity(max_size);

    let single_index_converter = SingleIndexConverter::new(mesh);

    matrix_entries
        .iter()
        .enumerate()
        .for_each(|(block_id, block_matrix_entries)| {
            let dim = block_matrix_entries.dim();

            // TODO optimize traversal based on matrix storage format for coord
            // access down the line

            let mut func = |target_index: BlockPointIndex, prop: &PointProps| match prop {
                PointProps::Connect { donor_index } => {
                    let source_block_point_index =
                        single_index_converter.get_block_point_index(*donor_index);
                    updates.push((target_index, source_block_point_index, None));
                }
                PointProps::ConnectPeriodic {
                    donor_index,
                    translation,
                } => {
                    let source_block_point_index =
                        single_index_converter.get_block_point_index(*donor_index);
                    updates.push((
                        target_index,
                        source_block_point_index,
                        Some(translation.clone()),
                    ));
                }
                _ => (),
            };

            // loop left
            for i in 0..dim.1 {
                func(
                    BlockPointIndex::new(block_id, (0, i)),
                    &block_matrix_entries[[0, i]].prop,
                );
            }

            // loop right
            for i in 0..dim.1 {
                func(
                    BlockPointIndex::new(block_id, (dim.0 - 1, i)),
                    &block_matrix_entries[[dim.0 - 1, i]].prop,
                );
            }

            // loop top
            for i in 1..(dim.0 - 1) {
                func(
                    BlockPointIndex::new(block_id, (i, 0)),
                    &block_matrix_entries[[i, 0]].prop,
                );
            }

            // loop bottom
            for i in 1..(dim.0 - 1) {
                func(
                    BlockPointIndex::new(block_id, (i, dim.1 - 1)),
                    &block_matrix_entries[[i, dim.1 - 1]].prop,
                );
            }
        });

    updates.shrink_to_fit();
    updates
}

/// represents the info needed to add the periodic ghost point to the matrix
#[derive(Debug)]
struct PeriodicGhostPointMatrixData {
    ghost_point_index: usize,
    donor_index: usize,
    translation: Vec2d,
}

impl PeriodicGhostPointMatrixData {
    fn new(
        ghost_point_index: usize,
        donor_index: usize,
        translation: Vec2d,
    ) -> PeriodicGhostPointMatrixData {
        PeriodicGhostPointMatrixData {
            ghost_point_index,
            donor_index,
            translation,
        }
    }
}

fn collect_periodic_ghost_points(
    matrix_entries: &Vec<Array2<MatrixEntry>>,
    dof: usize,
    num_mesh_points: usize,
) -> Vec<PeriodicGhostPointMatrixData> {
    let size = dof - num_mesh_points;
    let mut periodic_ghost_points = Vec::<PeriodicGhostPointMatrixData>::with_capacity(size);

    matrix_entries
        .iter()
        .enumerate()
        .for_each(|(block_id, block_matrix_entries)| {
            let dim = block_matrix_entries.dim();

            let mut func = |i, j| {
                if matrix_entries[block_id][[i, j]].index >= num_mesh_points
                    && matrix_entries[block_id][[i, j]].index < usize::MAX
                {
                    if let PointProps::ConnectPeriodic {
                        donor_index,
                        translation,
                    } = matrix_entries[block_id][[i, j]].prop
                    {
                        periodic_ghost_points.push(PeriodicGhostPointMatrixData::new(
                            matrix_entries[block_id][[i, j]].index,
                            donor_index,
                            translation,
                        ));
                    } else {
                        panic!("Ghost point to add to matrix is not periodic.");
                    }
                }
            };

            // loop left
            for i in 0..dim.1 {
                func(0, i);
            }

            // loop right
            for i in 0..dim.1 {
                func(dim.0 - 1, i);
            }

            // loop top
            for i in 1..(dim.0 - 1) {
                func(i, 0);
            }

            // loop bottom
            for i in 1..(dim.0 - 1) {
                func(i, dim.1 - 1);
            }
        });

    assert!(periodic_ghost_points.len() == size);
    periodic_ghost_points
}

pub fn smooth_mesh(mesh: &mut Mesh) -> Result<(), Box<dyn Error>> {
    let iterations = 20;

    // fields with ghost points
    let matrix_data = MatrixEntries::new(mesh);
    let entries = &matrix_data.data;

    if log_enabled!(Level::Trace) {
        // write matrix_entries to file

        use std::io::Write;
        let mut file = std::fs::File::create("matrix_entries.txt").unwrap();
        for (block_idx, block) in entries.iter().enumerate() {
            writeln!(file, "block {}", block_idx).unwrap();
            for ((i, j), entry) in block.indexed_iter() {
                writeln!(file, "{}, ({}, {}) = {:?}", block_idx, i, j, entry).unwrap();
            }
        }
    }

    let peridic_ghost_points =
        collect_periodic_ghost_points(&entries, matrix_data.dof, mesh.points());

    let mut coords = coordinates_with_ghost_points(&mesh);
    let ghost_point_coord_updates = compute_ghost_point_updates(mesh, &entries);

    // allocations
    let dof = matrix_data.dof;
    let mut lhs = SparseTriplet::new(dof, dof, dof * 9, Symmetry::No)?;
    let mut rhs_x = Vector::new(dof);
    let mut rhs_y = Vector::new(dof);

    // iterate
    for n in 0..iterations {
        print!("  iteration\t{n}");

        update_ghost_point_locations(&mut coords, &ghost_point_coord_updates);

        if log_enabled!(Level::Trace) {
            // write coords to file

            use std::io::Write;
            let mut file = std::fs::File::create("coords.txt").unwrap();
            for (block_idx, block) in coords.iter().enumerate() {
                writeln!(file, "block {}", block_idx).unwrap();
                for ((i, j), entry) in block.indexed_iter() {
                    writeln!(file, "{}, ({}, {}) = {:?}", block_idx, i, j, entry).unwrap();
                }
            }

            // write periodic ghost point entries to file

            let mut file = std::fs::File::create("periodic_ghost_points.txt").unwrap();
            peridic_ghost_points.iter().for_each(|line| {
                writeln!(file, "{:?}", line).unwrap();
            });
        }

        lhs.reset();
        rhs_x.fill(0.0);
        rhs_y.fill(0.0);

        // assemble LHS and RHS

        // add all mesh points
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

        // add all periodic ghost points
        peridic_ghost_points.iter().for_each(
            |PeriodicGhostPointMatrixData {
                 ghost_point_index,
                 donor_index,
                 translation,
             }| {
                let index = *ghost_point_index;

                lhs.put(index, index, 1.0).unwrap();
                lhs.put(index, *donor_index, -1.0).unwrap();
                rhs_x[index] = translation.0;
                rhs_y[index] = translation.1;
            },
        );

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

    // update (internal) mesh coordinates
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
