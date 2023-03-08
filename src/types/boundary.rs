// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{EdgeIndex, SegmentFunction};
use crate::{Mesh, Scalar, Vec2d};
use float_cmp::approx_eq;
use ndarray::{s, Array2, ArrayView1, ArrayViewMut1};
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

        let start_idx: usize = edge_segments[0..start]
            .iter()
            .map(|seg| seg.len())
            .sum::<usize>()
            - edge_segments[0..start].len();

        let len: usize = segments.iter().map(|seg| seg.len()).sum();
        let end_idx = start_idx + len - segments.len();

        Self {
            block,
            edge,
            start: start_idx,
            end: end_idx,
        }
    }

    pub fn reverse(mut self) -> Self {
        std::mem::swap(&mut self.start, &mut self.end);
        self
    }

    pub fn iter(&self) -> Box<dyn Iterator<Item = usize>> {
        if self.end > self.start {
            Box::new((self.start..=self.end).into_iter())
        } else {
            Box::new((self.end..=self.start).rev().into_iter())
        }
    }

    pub fn iter_inner(&self) -> Box<dyn Iterator<Item = usize>> {
        if self.end > self.start {
            Box::new((self.start + 1..=self.end - 1).into_iter())
        } else {
            Box::new((self.end + 1..=self.start - 1).rev().into_iter())
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
        const EPSILON: f64 = 1e-10;

        // check input
        let block_0 = &mesh.blocks[ranges.0.block].coords;
        let block_1 = &mesh.blocks[ranges.1.block].coords;

        let edge_0 = edge_to_view(block_0, ranges.0.edge, ranges.0.start, ranges.0.end);
        let edge_1 = edge_to_view(block_1, ranges.1.edge, ranges.1.start, ranges.1.end);

        edge_0
            .indexed_iter()
            .zip(edge_1.indexed_iter())
            .for_each(|((i_0, x_0), (i_1, x_1))| {
                assert!(
                    approx_eq!(&Vec2d, x_0, x_1, epsilon = EPSILON),
                    "Non-matching Coordinates for Connection:\n
                    {:?} at index {} with {}\n
                    {:?} at index {} with {}",
                    ranges.0,
                    i_0,
                    x_0,
                    ranges.1,
                    i_1,
                    x_1
                )
            });

        Self(ranges.0, ranges.1)
    }
}

pub fn edge_view<'a, T>(
    arrays: &'a Vec<Array2<T>>,
    range: &BlockBoundaryRange,
) -> ArrayView1<'a, T> {
    let array = &arrays[range.block];
    let (size_i, size_j) = array.dim();
    let start = range.start;
    let end = range.end;
    match range.edge {
        EdgeIndex::IMin => array.slice(s![start..=end, 0]),
        EdgeIndex::IMax => array.slice(s![start..=end, size_j - 1]),
        EdgeIndex::JMin => array.slice(s![0, start..=end]),
        EdgeIndex::JMax => array.slice(s![size_i - 1, start..=end]),
    }
}

pub fn edge_view_mut<'a, T>(
    arrays: &'a mut Vec<Array2<T>>,
    range: &BlockBoundaryRange,
) -> ArrayViewMut1<'a, T> {
    let array = &mut arrays[range.block];
    let (size_i, size_j) = array.dim();
    let start = range.start;
    let end = range.end;
    match range.edge {
        EdgeIndex::IMin => array.slice_mut(s![start..=end, 0]),
        EdgeIndex::IMax => array.slice_mut(s![start..=end, size_j - 1]),
        EdgeIndex::JMin => array.slice_mut(s![0, start..=end]),
        EdgeIndex::JMax => array.slice_mut(s![size_i - 1, start..=end]),
    }
}

pub fn edge_to_view<T>(
    array: &Array2<T>,
    edge: EdgeIndex,
    start: usize,
    end: usize,
) -> ArrayView1<T> {
    let (size_i, size_j) = array.dim();
    match edge {
        EdgeIndex::IMin => array.slice(s![start..=end, 0]),
        EdgeIndex::IMax => array.slice(s![start..=end, size_j - 1]),
        EdgeIndex::JMin => array.slice(s![0, start..=end]),
        EdgeIndex::JMax => array.slice(s![size_i - 1, start..=end]),
    }
}

// pub fn edge_to_view_mut<T>(
//     array: &mut Array2<T>,
//     edge: EdgeIndex,
//     start: usize,
//     end: usize,
// ) -> ArrayViewMut1<T> {
//     let (size_i, size_j) = array.dim();
//     match edge {
//         EdgeIndex::IMin => array.slice_mut(s![start..=end, 0]),
//         EdgeIndex::IMax => array.slice_mut(s![start..=end, size_j - 1]),
//         EdgeIndex::JMin => array.slice_mut(s![0, start..=end]),
//         EdgeIndex::JMax => array.slice_mut(s![size_i - 1, start..=end]),
//     }
// }
