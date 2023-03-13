// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{EdgeIndex, SegmentFunction};
use crate::{Mesh, Scalar, Vec2d};
use float_cmp::approx_eq;
use ndarray::{s, Array2, ArrayViewMut1};
use std::slice::SliceIndex;
use subslice_index::subslice_index;

#[derive(Debug, Clone)]
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
            Box::new(self.start + 1..=self.end - 1)
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
pub struct BlockBoundaryRangeNew {
    pub block: usize,
    pub start: Array2<isize>,
    pub end: Array2<isize>,
}

impl BlockBoundaryRangeNew {
    pub fn new(mesh: &Mesh, range: BlockBoundaryRange) -> Self {
        let block = range.block;

        let mut start = Array2::<isize>::zeros((2, 1));
        let mut end = Array2::<isize>::zeros((2, 1));

        match range.edge {
            EdgeIndex::IMin => {
                // j = 0, set i for start and end
                start[[0, 0]] = range.start as isize;
                end[[0, 0]] = range.end as isize;
            }
            EdgeIndex::IMax => {
                // j = max , set i for start and end
                start[[0, 0]] = range.start as isize;
                start[[1, 0]] = (mesh.blocks[block].points()[1] - 1) as isize;

                end[[0, 0]] = range.end as isize;
                end[[1, 0]] = (mesh.blocks[block].points()[1] - 1) as isize;
            }
            EdgeIndex::JMin => {
                // i = 0, set j for start and end
                start[[1, 0]] = range.start as isize;
                end[[1, 0]] = range.end as isize;
            }
            EdgeIndex::JMax => {
                // i = max , set j for start and end
                start[[0, 0]] = (mesh.blocks[block].points()[0] - 1) as isize;
                start[[1, 0]] = range.start as isize;

                end[[0, 0]] = (mesh.blocks[block].points()[0] - 1) as isize;
                end[[1, 0]] = range.end as isize;
            }
        }

        Self { block, start, end }
    }

    pub fn iter(&self) -> BlockBoundaryRangeNewIter {
        BlockBoundaryRangeNewIter::new(self)
    }

    pub fn first(&self) -> (usize, usize) {
        (self.start[[0, 0]] as usize, self.start[[1, 0]] as usize)
    }

    pub fn last(&self) -> (usize, usize) {
        (self.end[[0, 0]] as usize, self.end[[1, 0]] as usize)
    }
}

pub struct BlockBoundaryRangeNewIter<'a> {
    range: &'a BlockBoundaryRangeNew,
    dim: usize,
    step: isize,
    steps: usize,
    count: usize,
}

impl<'a> BlockBoundaryRangeNewIter<'a> {
    fn new(range: &'a BlockBoundaryRangeNew) -> Self {
        let mut dim = 1;

        if range.start.first().unwrap() == range.end.first().unwrap() {
            dim = 0;
        }

        let mut step = 1;

        if range.start[[dim, 0]] > range.end[[dim, 0]] {
            step = -1;
        }

        let steps = if step < 0 {
            (range.start[[dim, 0]] - range.end[[dim, 0]]) as usize
        } else {
            (range.end[[dim, 0]] - range.start[[dim, 0]]) as usize
        };

        Self {
            range,
            dim,
            step,
            steps,
            count: 0,
        }
    }
}

impl<'a> Iterator for BlockBoundaryRangeNewIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.count += 1;

        if self.count <= self.steps {
            let mut index = self.range.start.clone();
            index[[self.dim, 0]] += self.step;

            Some((index[[0, 0]] as usize, index[[1, 0]] as usize))
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct BlockConnection {
    pub donor: BlockBoundaryRangeNew,
    pub receiver: BlockBoundaryRangeNew,

    // 2x2 transformation matrix, see
    // https://cgns.github.io/CGNS_docs_current/sids/cnct.html#Transform
    transform: ndarray::Array2<isize>,
}

impl BlockConnection {
    pub fn new(mesh: &Mesh, ranges: (BlockBoundaryRange, BlockBoundaryRange)) -> Self {
        const EPSILON: f64 = 1e-10;

        let donor = BlockBoundaryRangeNew::new(mesh, ranges.0.clone());
        let receiver = BlockBoundaryRangeNew::new(mesh, ranges.1.clone());

        // create transformation matrix
        let mut transform = Array2::<isize>::zeros((2, 2));

        let donor_parallel_sign = if ranges.0.start > ranges.0.end { -1 } else { 1 };
        let rec_parallel_sign = if ranges.1.start > ranges.1.end { -1 } else { 1 };

        let parallel_sign = donor_parallel_sign * rec_parallel_sign;

        let donor_orth_sign = match ranges.0.edge {
            EdgeIndex::IMin => -1,
            EdgeIndex::IMax => 1,
            EdgeIndex::JMin => -1,
            EdgeIndex::JMax => 1,
        };

        let rec_orth_sign = match ranges.1.edge {
            EdgeIndex::IMin => -1,
            EdgeIndex::IMax => 1,
            EdgeIndex::JMin => -1,
            EdgeIndex::JMax => 1,
        };

        let orth_sign = donor_orth_sign * rec_orth_sign;

        if ranges.1.edge == EdgeIndex::IMin || ranges.1.edge == EdgeIndex::IMax {
            transform[[0, 0]] = parallel_sign;
            transform[[1, 1]] = orth_sign;
        } else {
            transform[[0, 1]] = orth_sign;
            transform[[1, 0]] = parallel_sign;
        }

        let connection = Self {
            donor,
            receiver,
            transform,
        };

        // check overlap

        let block_0 = &mesh.blocks[ranges.0.block].coords;
        let block_1 = &mesh.blocks[ranges.1.block].coords;

        connection
            .donor
            .iter()
            .zip(connection.receiver.iter())
            .for_each(|(donor_idx, rec_idx)| {
                let x_0 = block_0[[donor_idx.0, donor_idx.1]];
                let x_1 = block_1[[rec_idx.0, rec_idx.1]];
                assert!(
                    approx_eq!(&Vec2d, &x_0, &x_1, epsilon = EPSILON),
                    "Non-matching Coordinates for Connection:\n
                    {:?} at index {:?} with {}\n
                    {:?} at index {:?} with {}",
                    ranges.0,
                    donor_idx,
                    x_0,
                    ranges.1,
                    rec_idx,
                    x_1
                )
            });

        // check transform from donor to receiver w/ ghost layer on donor, i.e.
        // first internal layer on receiver

        let add = match ranges.0.edge {
            EdgeIndex::IMin => (0, -1),
            EdgeIndex::IMax => (0, mesh.blocks[ranges.0.block].points()[1] as isize),
            EdgeIndex::JMin => (-1, 0),
            EdgeIndex::JMax => (mesh.blocks[ranges.0.block].points()[0] as isize, 0),
        };

        let add_rec: (isize, isize) = match ranges.1.edge {
            EdgeIndex::IMin => (0, 1),
            EdgeIndex::IMax => (0, -1),
            EdgeIndex::JMin => (1, 0),
            EdgeIndex::JMax => (-1, 0),
        };

        connection
            .donor
            .iter()
            .zip(connection.receiver.iter())
            .for_each(|((i, j), (i_rec, j_rec))| {
                let i = i as isize + add.0;
                let j = j as isize + add.1;

                let ghost_point_comp = connection.get_index_in_receiver_block(&[i, j]);

                let ghost_point = (i_rec as isize + add_rec.0, j_rec as isize + add_rec.1);

                assert_eq!(ghost_point_comp, ghost_point);
            });

        connection
    }

    pub fn get_index_in_receiver_block(&self, donor_index: &[isize; 2]) -> (isize, isize) {
        let donor_index = Array2::from_shape_vec((2, 1), donor_index.to_vec()).unwrap();

        let receiver_index = self.transform.clone() * (donor_index - self.donor.start.clone())
            + self.receiver.start.clone();

        (receiver_index[[0, 0]], receiver_index[[1, 0]])
    }
}

// pub fn edge_view<'a, T>(
//     arrays: &'a Vec<Array2<T>>,
//     range: &BlockBoundaryRange,
// ) -> ArrayView1<'a, T> {
//     let array = &arrays[range.block];
//     let (size_i, size_j) = array.dim();
//     let start = range.start;
//     let end = range.end;
//     match range.edge {
//         EdgeIndex::IMin => array.slice(s![start..=end, 0]),
//         EdgeIndex::IMax => array.slice(s![start..=end, size_j - 1]),
//         EdgeIndex::JMin => array.slice(s![0, start..=end]),
//         EdgeIndex::JMax => array.slice(s![size_i - 1, start..=end]),
//     }
// }

pub fn edge_view_mut<'a, T>(
    arrays: &'a mut Vec<Array2<T>>,
    range: &BlockBoundaryRange,
) -> ArrayViewMut1<'a, T> {
    let array = &mut arrays[range.block];
    let (size_i, size_j) = array.dim();
    let start = range.start;
    let end = range.end;

    if start < end {
        match range.edge {
            EdgeIndex::IMin => array.slice_mut(s![start..=end, 0]),
            EdgeIndex::IMax => array.slice_mut(s![start..=end, size_j - 1]),
            EdgeIndex::JMin => array.slice_mut(s![0, start..=end]),
            EdgeIndex::JMax => array.slice_mut(s![size_i - 1, start..=end]),
        }
    } else {
        match range.edge {
            EdgeIndex::IMin => array.slice_mut(s![start..=end, 0]),
            EdgeIndex::IMax => array.slice_mut(s![start..=end, size_j - 1]),
            EdgeIndex::JMin => array.slice_mut(s![0, start..=end]),
            EdgeIndex::JMax => array.slice_mut(s![size_i - 1, start..=end]),
        }
    }
}

// pub fn edge_to_view<T>(
//     array: &Array2<T>,
//     edge: EdgeIndex,
//     start: usize,
//     end: usize,
// ) -> ArrayView1<T> {
//     let (size_i, size_j) = array.dim();
//     match edge {
//         EdgeIndex::IMin => array.slice(s![start..=end, 0]),
//         EdgeIndex::IMax => array.slice(s![start..=end, size_j - 1]),
//         EdgeIndex::JMin => array.slice(s![0, start..=end]),
//         EdgeIndex::JMax => array.slice(s![size_i - 1, start..=end]),
//     }
// }

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
