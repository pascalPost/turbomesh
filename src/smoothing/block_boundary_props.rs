// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use std::cell::Cell;

use crate::types::{BlockBoundary, BlockBoundaryRange, BlockConnection};
use crate::{Block2d, Mesh};

#[derive(Debug, Clone, PartialEq)]
pub struct ConnectionData {
    block: usize,
    point: (usize, usize),
}

impl ConnectionData {
    pub fn new(block: usize, point: (usize, usize)) -> Self {
        Self { block, point }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BlockBoundaryPointProp {
    Undefined,
    Fixed,
    Connected(ConnectionData),
}

/// stores the block boundary point edge mapping for all blocks
pub struct BoundaryProps {
    pub blocks: Vec<BlockBoundaryProps<Vec<(BlockBoundaryPointProp, usize)>>>,
}

impl BoundaryProps {
    pub fn new(mesh: &Mesh) -> Self {
        let mut blocks =
            Vec::<BlockBoundaryProps<Vec<(BlockBoundaryPointProp, usize)>>>::with_capacity(
                mesh.blocks.len(),
            );

        mesh.blocks.iter().for_each(|block| {
            blocks.push(
                BlockBoundaryProps::<Vec<(BlockBoundaryPointProp, usize)>>::new(block, || {
                    Vec::<(BlockBoundaryPointProp, usize)>::with_capacity(4)
                }),
            );
        });

        mesh.edges
            .iter()
            .enumerate()
            .for_each(|(edge_idx, edge)| match edge {
                BlockBoundary::Connection(connection) => {
                    add_connected_points(&mut blocks, connection, edge_idx)
                }
                BlockBoundary::PeriodicConnection(per) => {
                    add_connected_points(&mut blocks, &per.connection, edge_idx)
                }
                BlockBoundary::Inlet(range)
                | BlockBoundary::Outlet(range)
                | BlockBoundary::Wall(range) => add_fixed_points(&mut blocks, &range, edge_idx),
            });

        blocks.iter_mut().for_each(|bound| bound.shrink_to_fit());

        Self { blocks }
    }
}

fn add_fixed_points(
    blocks: &mut Vec<BlockBoundaryProps<Vec<(BlockBoundaryPointProp, usize)>>>,
    range: &BlockBoundaryRange,
    edge_idx: usize,
) {
    let block_boundary = &mut blocks[range.block];
    range.iter().for_each(|p| {
        block_boundary
            .get_mut(p)
            .push((BlockBoundaryPointProp::Fixed, edge_idx));
    });
}

fn add_connected_points(
    blocks: &mut Vec<BlockBoundaryProps<Vec<(BlockBoundaryPointProp, usize)>>>,
    connection: &BlockConnection,
    edge_idx: usize,
) {
    connection
        .donor
        .iter()
        .zip(connection.receiver.iter())
        .for_each(|(donor, rec)| {
            blocks[connection.donor.block].get_mut(donor).push((
                BlockBoundaryPointProp::Connected(ConnectionData::new(
                    connection.receiver.block,
                    rec,
                )),
                edge_idx,
            ));

            blocks[connection.receiver.block].get_mut(rec).push((
                BlockBoundaryPointProp::Connected(ConnectionData::new(
                    connection.donor.block,
                    donor,
                )),
                edge_idx,
            ));
        });
}

/// stores for every block boundary point the edges it is contained in
pub struct BlockBoundaryProps<T> {
    dim: [usize; 2],

    // saves for every point on the block boundary the mesh edges the point is
    // contained in
    pub points: Vec<T>,
}

impl<T> BlockBoundaryProps<T> {
    fn new<F>(block: &Block2d, init: F) -> Self
    where
        F: FnMut() -> T,
    {
        let dim = block.points();
        let mut points = Vec::<T>::new();
        points.resize_with(dim[0] * 2 + dim[1] * 2 - 4, init);
        Self { dim, points }
    }

    pub fn get(&self, point: (usize, usize)) -> &T {
        let point_idx = boundary_point_index(point, (self.dim[0], self.dim[1])).unwrap();
        &self.points[point_idx]
    }

    pub fn get_mut(&mut self, point: (usize, usize)) -> &mut T {
        let point_idx = boundary_point_index(point, (self.dim[0], self.dim[1])).unwrap();
        &mut self.points[point_idx]
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.points.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.points.iter_mut()
    }

    pub fn point_index(&self, boundary_point_index: usize) -> Option<(usize, usize)> {
        boundary_point_to_point(boundary_point_index, (self.dim[0], self.dim[1]))
    }

    fn shrink_to_fit(&mut self) {
        self.points.shrink_to_fit();
    }
}

/// returns a flat array boundary index for the given point on the block of
/// boudary dimensions dim
///
/// example block of size 6x4:
///          (0,3) (1,3) (2,3) (3,3) (4,3) (5,3)
///          13    12    11    10    09    08
/// (0,3)    x     x     x     x     x     x     (5,3)
/// (0,2) 14 x                             x  07 (5,2)
/// (0,1) 15 x                             x  06 (5,1)
/// (0,0)    x     x     x     x     x     x     (5,0)
///          00    01    02    03    04    05
///         (0,0) (1,0) (2,0) (3,0) (4,0) (5,0)
///
/// ```
/// use turbomesh::smoothing::block_boundary_props::boundary_point_index;
///
/// let dim = (6,4);
/// assert_eq!(boundary_point_index((0, 0), dim), Some(0));
/// assert_eq!(boundary_point_index((1, 0), dim), Some(1));
/// assert_eq!(boundary_point_index((2, 0), dim), Some(2));
/// assert_eq!(boundary_point_index((3, 0), dim), Some(3));
/// assert_eq!(boundary_point_index((4, 0), dim), Some(4));
/// assert_eq!(boundary_point_index((5, 0), dim), Some(5));
///
/// assert_eq!(boundary_point_index((5, 0), dim), Some(5));
/// assert_eq!(boundary_point_index((5, 1), dim), Some(6));
/// assert_eq!(boundary_point_index((5, 2), dim), Some(7));
/// assert_eq!(boundary_point_index((5, 3), dim), Some(8));
///
/// assert_eq!(boundary_point_index((5, 3), dim), Some(8));
/// assert_eq!(boundary_point_index((4, 3), dim), Some(9));
/// assert_eq!(boundary_point_index((3, 3), dim), Some(10));
/// assert_eq!(boundary_point_index((2, 3), dim), Some(11));
/// assert_eq!(boundary_point_index((1, 3), dim), Some(12));
/// assert_eq!(boundary_point_index((0, 3), dim), Some(13));
///
/// assert_eq!(boundary_point_index((0, 3), dim), Some(13));
/// assert_eq!(boundary_point_index((0, 2), dim), Some(14));
/// assert_eq!(boundary_point_index((0, 1), dim), Some(15));
/// ```
pub fn boundary_point_index(point: (usize, usize), dim: (usize, usize)) -> Option<usize> {
    if point.1 == 0 {
        if point.0 < dim.0 {
            Some(point.0)
        } else {
            None
        }
    } else if point.0 == dim.0 - 1 {
        if point.1 < dim.1 {
            Some(point.0 + point.1)
        } else {
            None
        }
    } else if point.1 == dim.1 - 1 {
        if point.0 < dim.0 {
            Some(2 * (dim.0 - 1) - point.0 + point.1)
        } else {
            None
        }
    } else if point.0 == 0 {
        if point.1 < dim.1 {
            Some((2 * dim.0 + 2 * dim.1 - 4) - point.1)
        } else {
            None
        }
    } else {
        None
    }
}

/// returns the point index given the flat array boundary index for the block of
/// boudary dimensions dim
///
/// example block of size 6x4:
///          (0,3) (1,3) (2,3) (3,3) (4,3) (5,3)
///          13    12    11    10    09    08
/// (0,3)    x     x     x     x     x     x     (5,3)
/// (0,2) 14 x                             x  07 (5,2)
/// (0,1) 15 x                             x  06 (5,1)
/// (0,0)    x     x     x     x     x     x     (5,0)
///          00    01    02    03    04    05
///         (0,0) (1,0) (2,0) (3,0) (4,0) (5,0)
///
/// ```
/// use turbomesh::smoothing::block_boundary_props::boundary_point_to_point;
///
/// let dim = (6,4);
/// assert_eq!(boundary_point_to_point(0, dim), Some((0, 0)));
/// assert_eq!(boundary_point_to_point(1, dim), Some((1, 0)));
/// assert_eq!(boundary_point_to_point(2, dim), Some((2, 0)));
/// assert_eq!(boundary_point_to_point(3, dim), Some((3, 0)));
/// assert_eq!(boundary_point_to_point(4, dim), Some((4, 0)));
/// assert_eq!(boundary_point_to_point(5, dim), Some((5, 0)));
///
/// assert_eq!(boundary_point_to_point(5, dim), Some((5, 0)));
/// assert_eq!(boundary_point_to_point(6, dim), Some((5, 1)));
/// assert_eq!(boundary_point_to_point(7, dim), Some((5, 2)));
/// assert_eq!(boundary_point_to_point(8, dim), Some((5, 3)));
///
/// assert_eq!(boundary_point_to_point(8, dim), Some((5, 3)));
/// assert_eq!(boundary_point_to_point(9, dim), Some((4, 3)));
/// assert_eq!(boundary_point_to_point(10, dim), Some((3, 3)));
/// assert_eq!(boundary_point_to_point(11, dim), Some((2, 3)));
/// assert_eq!(boundary_point_to_point(12, dim), Some((1, 3)));
/// assert_eq!(boundary_point_to_point(13, dim), Some((0, 3)));
///
/// assert_eq!(boundary_point_to_point(13, dim), Some((0, 3)));
/// assert_eq!(boundary_point_to_point(14, dim), Some((0, 2)));
/// assert_eq!(boundary_point_to_point(15, dim), Some((0, 1)));
/// assert_eq!(boundary_point_to_point(16, dim), None);
/// ```
pub fn boundary_point_to_point(
    boundary_point_index: usize,
    dim: (usize, usize),
) -> Option<(usize, usize)> {
    if boundary_point_index < dim.0 {
        Some((boundary_point_index, 0))
    } else if boundary_point_index < dim.0 + dim.1 - 1 {
        Some((dim.0 - 1, boundary_point_index - dim.0 + 1))
    } else if boundary_point_index < 2 * dim.0 + dim.1 - 2 {
        Some((2 * dim.0 + dim.1 - 3 - boundary_point_index, dim.1 - 1))
    } else if boundary_point_index < 2 * dim.0 + 2 * dim.1 - 4 {
        Some((0, 2 * dim.0 + 2 * dim.1 - 4 - boundary_point_index))
    } else {
        None
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BlockBoundaryPointSolverProp {
    Undefined,
    Fixed,
    Solved(Vec<ConnectionData>),
    Connected(ConnectionData),
}

/// compute the matrix treatment for every mesh point
pub fn block_boundary_points_solver_props(
    mesh: &Mesh,
) -> Vec<BlockBoundaryProps<BlockBoundaryPointSolverProp>> {
    let boundary_props = BoundaryProps::new(mesh);

    // init all points to undefined
    let mut data = mesh
        .blocks
        .iter()
        .map(|block| BlockBoundaryProps::new(block, || BlockBoundaryPointSolverProp::Undefined))
        .collect::<Vec<_>>();

    // set all fixed points to fixed
    boundary_props.blocks.iter().zip(data.iter_mut()).for_each(
        |(block_data, block_solver_data)| {
            block_data
                .iter()
                .zip(block_solver_data.iter_mut())
                .for_each(|(point_data, point_solver_data)| {
                    for p in point_data.iter() {
                        assert!(p.0 != BlockBoundaryPointProp::Undefined);
                        if let BlockBoundaryPointProp::Fixed = p.0 {
                            *point_solver_data = BlockBoundaryPointSolverProp::Fixed;
                            break;
                        }
                    }
                });
        },
    );

    for (block_idx, block) in boundary_props.blocks.iter().enumerate() {
        for (point_idx, point) in block.iter().enumerate() {
            // check if point props need to be defined
            if data[block_idx].points[point_idx] == BlockBoundaryPointSolverProp::Undefined {
                // loop over all connected points and check for a fixed point. If
                // there is a fixed point, set all the current point and all
                // connected points to fixed
                let mut is_fixed = false;
                for p in point.iter() {
                    if let BlockBoundaryPointProp::Connected(con) = &p.0 {
                        if *data[con.block].get(con.point) == BlockBoundaryPointSolverProp::Fixed {
                            is_fixed = true;
                            break;
                        }
                    }
                }
                if is_fixed {
                    data[block_idx].points[point_idx] = BlockBoundaryPointSolverProp::Fixed;
                    for p in point.iter() {
                        if let BlockBoundaryPointProp::Connected(con) = &p.0 {
                            *data[con.block].get_mut(con.point) =
                                BlockBoundaryPointSolverProp::Fixed;
                        }
                    }
                } else {
                    // set the current point to be solved and the connected
                    // points to be connected
                    let mut connected_points = Vec::<ConnectionData>::with_capacity(4);
                    for p in point.iter() {
                        if let BlockBoundaryPointProp::Connected(con) = &p.0 {
                            *data[con.block].get_mut(con.point) =
                                BlockBoundaryPointSolverProp::Connected(ConnectionData::new(
                                    block_idx,
                                    block.point_index(point_idx).unwrap(),
                                ));

                            connected_points.push(ConnectionData::new(con.block, con.point));
                        } else {
                            panic!("point is not connected");
                        }
                    }

                    connected_points.shrink_to_fit();

                    data[block_idx].points[point_idx] =
                        BlockBoundaryPointSolverProp::Solved(connected_points);
                }
            }
        }
    }

    // boundary_props
    //     .blocks
    //     .iter()
    //     .enumerate()
    //     .for_each(|(block_idx, block)| {
    //         block.iter().enumerate().for_each(|(point_idx, point)| {
    //             // loop over all connected points

    //             // check if the current point is to be fixed (if one connected point is fixed)
    //             let mut is_fixed = false;

    //             for p in point.iter() {
    //                 if let BlockBoundaryPointProp::Connected(con) = &p.0 {
    //                     if *data[con.block].get(con.point) == BlockBoundaryPointSolverProp::Fixed {
    //                         is_fixed = true;
    //                         break;
    //                     }
    //                 }
    //             }

    //             if is_fixed {
    //                 data[block_idx].points[point_idx] = BlockBoundaryPointSolverProp::Fixed;
    //             } else {
    //                 // the point is to be solved

    //                 // check if one of the points it already set to be solved

    //                 let mut is_already_solved = false;

    //                 for p in point.iter() {
    //                     if let BlockBoundaryPointProp::Connected(con) = &p.0 {
    //                         if *data[con.block].get(con.point)
    //                             == BlockBoundaryPointSolverProp::Solved
    //                         {
    //                             is_already_solved = true;
    //                             break;
    //                         }
    //                     }
    //                 }

    //                 if !is_already_solved {
    //                     // set the current point to be solved and all other
    //                     // points to be connected
    //                 }
    //             }
    //         })
    //     });

    // // check all connected points for fixed points and set these to fixed
    // {
    //     // TODO preallocate to gain speed
    //     let mut changes = Vec::<(usize, usize)>::new();

    //     boundary_props
    //         .blocks
    //         .iter()
    //         .zip(data.iter())
    //         .enumerate()
    //         .for_each(|(block_idx, (block_data, block_solver_data))| {
    //             block_data
    //                 .iter()
    //                 .zip(block_solver_data.iter())
    //                 .enumerate()
    //                 .for_each(|(point_idx, (point_data, point_solver_data))| {
    //                     if *point_solver_data != BlockBoundaryPointSolverProp::Fixed {
    //                         for p in point_data.iter() {
    //                             if let BlockBoundaryPointProp::Connected(con) = &p.0 {
    //                                 // check if connected point is fixed
    //                                 if *data[con.block].get(con.point)
    //                                     == BlockBoundaryPointSolverProp::Fixed
    //                                 {
    //                                     changes.push((block_idx, point_idx));
    //                                     break;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 });
    //         });

    //     changes.into_iter().for_each(|(block_idx, point_idx)| {
    //         data[block_idx].points[point_idx] = BlockBoundaryPointSolverProp::Fixed;
    //     });
    // }

    data
}
