// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{EdgeIndex, SegmentFunction};
use crate::{Mesh, Scalar};
use std::slice::SliceIndex;
use subslice_index::subslice_index;

#[derive(Debug)]
pub struct BlockBoundaryRange {
    pub block: usize,
    pub edge: EdgeIndex,
    pub start: usize,
    pub end: usize,
}

impl BlockBoundaryRange {
    pub fn new<Index>(mesh: &Mesh, block: usize, edge: EdgeIndex, segments: Index) -> Self
    where
        Index: SliceIndex<
            [Box<(dyn SegmentFunction + 'static)>],
            Output = [Box<(dyn SegmentFunction)>],
        >,
    {
        let edge_segments = mesh.blocks[block].segments(edge);
        let segments = &edge_segments[segments];

        let start = subslice_index(edge_segments.as_slice(), segments);

        let mut start_idx = 0;
        for i in 0..=start {
            start_idx += edge_segments[i].len();
        }

        let len: usize = segments.iter().map(|seg| seg.len()).sum();
        let end_idx = start_idx + len - 1;

        Self {
            block,
            edge,
            start: start_idx,
            end: end_idx,
        }
    }
}

// TODO implement non-fixed inlet and outlet BCs
#[derive(Debug)]
pub enum BlockBoundary {
    Connection(BlockConnection),
    PeriodicConnection {
        connection: (BlockBoundaryRange, BlockBoundaryRange),
        translation: Scalar,
    },
    Inlet(BlockBoundaryRange),
    Outlet(BlockBoundaryRange),
    Wall(BlockBoundaryRange),
}

#[derive(Debug)]
pub struct BlockConnection(pub BlockBoundaryRange, pub BlockBoundaryRange);

impl BlockConnection {
    pub fn new(mesh: &Mesh, ranges: (BlockBoundaryRange, BlockBoundaryRange)) -> Self {
        // check input
        &mesh.blocks[ranges.0.block];

        Self(ranges.0, ranges.1)
    }
}
