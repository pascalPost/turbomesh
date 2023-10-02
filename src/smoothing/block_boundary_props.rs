// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{BlockBoundary, BlockBoundaryRange, BlockConnection, PeriodicBlockConnection};
use crate::{Block2d, Mesh};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlockPointIndex {
    pub block: usize,
    pub point: (usize, usize),
}

impl BlockPointIndex {
    pub fn new(block: usize, point: (usize, usize)) -> Self {
        Self { block, point }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConnectionData {
    pub index: BlockPointIndex,
    pub edge: usize,
    pub donor: bool,
}

impl ConnectionData {
    fn new(index: BlockPointIndex, edge: usize, donor: bool) -> Self {
        Self { index, edge, donor }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BlockBoundaryPointProp {
    Undefined,
    Fixed,
    Connected(ConnectionData),
    Periodic(ConnectionData),
}

/// stores the block boundary point edge mapping for all blocks
#[derive(Debug, Clone)]
pub struct BoundaryProps {
    pub data: Vec<BlockBoundaryArray<Vec<BlockBoundaryPointProp>>>,
}

impl BoundaryProps {
    pub fn new(mesh: &Mesh) -> Self {
        let mut blocks = Vec::<BlockBoundaryArray<Vec<BlockBoundaryPointProp>>>::with_capacity(
            mesh.blocks.len(),
        );

        mesh.blocks.iter().for_each(|block| {
            blocks.push(BlockBoundaryArray::<Vec<BlockBoundaryPointProp>>::new(
                block,
                || Vec::<BlockBoundaryPointProp>::with_capacity(4),
            ));
        });

        mesh.edges
            .iter()
            .enumerate()
            .for_each(|(edge_idx, edge)| match edge {
                BlockBoundary::Connection(connection) => {
                    add_connected_points(&mut blocks, connection, edge_idx)
                }
                BlockBoundary::PeriodicConnection(per) => {
                    add_periodic_points(&mut blocks, &per, edge_idx)
                }
                BlockBoundary::Inlet(range)
                | BlockBoundary::Outlet(range)
                | BlockBoundary::Wall(range) => add_fixed_points(&mut blocks, &range),
            });

        blocks.iter_mut().for_each(|bound| bound.shrink_to_fit());

        Self { data: blocks }
    }
}

impl std::ops::Index<&BlockPointIndex> for BoundaryProps {
    type Output = Vec<BlockBoundaryPointProp>;
    fn index<'a>(&'a self, i: &BlockPointIndex) -> &'a Self::Output {
        &self.data[i.block].get(i.point).unwrap()
    }
}

impl fmt::Display for BoundaryProps {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.data
            .iter()
            .enumerate()
            .try_for_each(|(block_idx, block)| {
                writeln!(f, "block: {block_idx}")?;
                block
                    .iter()
                    .enumerate()
                    .try_for_each(|(point_idx, point)| -> Result<(), _> {
                        writeln!(f, " . point: {point_idx}")?;
                        writeln!(f, " . . index: {:?}", block.point_index(point_idx).unwrap())?;
                        writeln!(f, " . . data: {:#?}", point)
                        // ControlFlow::Continue(())
                        // Ok(())
                    })
                // println!(" . props: {:#?}", block);
                // Ok(())
            })

        // write!(f, "Point [{} {}]", self.x, self.y)
        // f.debug_struct("Point")
        //  .field("x", &self.x)
        //  .field("y", &self.y)
        //  .finish()
    }
}

/// adds the fixed points to the block boundary point edge mapping
fn add_fixed_points(
    blocks: &mut Vec<BlockBoundaryArray<Vec<BlockBoundaryPointProp>>>,
    range: &BlockBoundaryRange,
) {
    let block_boundary = &mut blocks[range.block];
    range.iter().for_each(|p| {
        block_boundary
            .get_mut(p)
            .unwrap()
            .push(BlockBoundaryPointProp::Fixed);
    });
}

/// adds the connected points to the block boundary point edge mapping
fn add_connected_points(
    blocks: &mut Vec<BlockBoundaryArray<Vec<BlockBoundaryPointProp>>>,
    connection: &BlockConnection,
    edge_idx: usize,
) {
    connection
        .donor
        .iter()
        .zip(connection.receiver.iter())
        .for_each(|(donor, rec)| {
            blocks[connection.donor.block].get_mut(donor).unwrap().push(
                BlockBoundaryPointProp::Connected(ConnectionData::new(
                    BlockPointIndex::new(connection.receiver.block, rec),
                    edge_idx,
                    true,
                )),
            );

            blocks[connection.receiver.block]
                .get_mut(rec)
                .unwrap()
                .push(BlockBoundaryPointProp::Connected(ConnectionData::new(
                    BlockPointIndex::new(connection.donor.block, donor),
                    edge_idx,
                    false,
                )));
        });
}

/// adds the periodic points to the block boundary point edge mapping, where the
/// donor is set to periodic and the receiver is set to connected
fn add_periodic_points(
    blocks: &mut Vec<BlockBoundaryArray<Vec<BlockBoundaryPointProp>>>,
    per: &PeriodicBlockConnection,
    edge_idx: usize,
) {
    let connection = &per.connection;

    connection
        .donor
        .iter()
        .zip(connection.receiver.iter())
        .for_each(|(donor, rec)| {
            blocks[connection.donor.block].get_mut(donor).unwrap().push(
                BlockBoundaryPointProp::Periodic(ConnectionData::new(
                    BlockPointIndex::new(connection.receiver.block, rec),
                    edge_idx,
                    true,
                )),
            );

            blocks[connection.receiver.block]
                .get_mut(rec)
                .unwrap()
                .push(BlockBoundaryPointProp::Periodic(ConnectionData::new(
                    BlockPointIndex::new(connection.donor.block, donor),
                    edge_idx,
                    false,
                )));
        });
}

pub struct MeshBoundaryArray<T> {
    pub blocks: Vec<BlockBoundaryArray<T>>,
}

impl<T> MeshBoundaryArray<T> {
    pub fn new<F>(mesh: &Mesh, init: F) -> Self
    where
        F: FnMut() -> T + Copy,
    {
        let mut blocks = Vec::<BlockBoundaryArray<T>>::with_capacity(mesh.blocks.len());

        mesh.blocks.iter().for_each(|block| {
            blocks.push(BlockBoundaryArray::<T>::new(block, init));
        });

        Self { blocks }
    }
}

/// block boundary point array
#[derive(Debug, Clone)]
pub struct BlockBoundaryArray<T> {
    /// block dimension
    pub dim: [usize; 2],

    /// vector containing every block boundary point
    pub points: Vec<T>,
}

impl<T> BlockBoundaryArray<T> {
    /// creates a new block boundary array given a block and an initializer function
    fn new<F>(block: &Block2d, init: F) -> Self
    where
        F: FnMut() -> T,
    {
        let dim = block.points();
        let mut points = Vec::<T>::new();
        points.resize_with(dim[0] * 2 + dim[1] * 2 - 4, init);
        Self { dim, points }
    }

    pub fn get(&self, point: (usize, usize)) -> Option<&T> {
        let point_idx = boundary_point_index(point, (self.dim[0], self.dim[1]))?;
        Some(&self.points[point_idx])
    }

    pub fn get_mut(&mut self, point: (usize, usize)) -> Option<&mut T> {
        let point_idx = boundary_point_index(point, (self.dim[0], self.dim[1]))?;
        Some(&mut self.points[point_idx])
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

    pub fn boundary_point_index(&self, point: (usize, usize)) -> Option<usize> {
        boundary_point_index(point, (self.dim[0], self.dim[1]))
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
