// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{EdgeIndex, SegmentFunction};
use crate::{Mesh, Scalar, Vec2d};
use float_cmp::approx_eq;
use ndarray::Array2;
use std::slice::SliceIndex;
use subslice_index::subslice_index;

// TODO implement non-fixed inlet and outlet BCs
#[derive(Debug)]
pub enum BlockBoundary {
    Connection(BlockConnection),
    PeriodicConnection(PeriodicBlockConnection),
    Inlet(BlockBoundaryRange),
    Outlet(BlockBoundaryRange),
    Wall(BlockBoundaryRange),
}

#[derive(Debug)]
pub struct BlockBoundaryRange {
    pub block: usize,
    pub start: Array2<isize>,
    pub end: Array2<isize>,
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

        let mut start = Array2::<isize>::zeros((2, 1));
        let mut end = Array2::<isize>::zeros((2, 1));

        match edge {
            EdgeIndex::IMin => {
                // j = 0, set i for start and end
                start[[0, 0]] = start_idx as isize;
                end[[0, 0]] = end_idx as isize;
            }
            EdgeIndex::IMax => {
                // j = max , set i for start and end
                start[[0, 0]] = start_idx as isize;
                start[[1, 0]] = (mesh.blocks[block].points()[1] - 1) as isize;

                end[[0, 0]] = end_idx as isize;
                end[[1, 0]] = (mesh.blocks[block].points()[1] - 1) as isize;
            }
            EdgeIndex::JMin => {
                // i = 0, set j for start and end
                start[[1, 0]] = start_idx as isize;
                end[[1, 0]] = end_idx as isize;
            }
            EdgeIndex::JMax => {
                // i = max , set j for start and end
                start[[0, 0]] = (mesh.blocks[block].points()[0] - 1) as isize;
                start[[1, 0]] = start_idx as isize;

                end[[0, 0]] = (mesh.blocks[block].points()[0] - 1) as isize;
                end[[1, 0]] = end_idx as isize;
            }
        }

        Self { block, start, end }
    }

    pub fn edge_type(&self) -> EdgeIndex {
        if self.start[[0, 0]] == self.end[[0, 0]] {
            if self.start[[0, 0]] == 0 {
                EdgeIndex::JMin
            } else {
                EdgeIndex::JMax
            }
        } else {
            if self.start[[1, 0]] == 0 {
                EdgeIndex::IMin
            } else {
                EdgeIndex::IMax
            }
        }
    }

    pub fn is_increasing(&self) -> bool {
        if self.start[[0, 0]] == self.end[[0, 0]] {
            self.start[[1, 0]] < self.end[[1, 0]]
        } else {
            self.start[[0, 0]] < self.end[[0, 0]]
        }
    }

    pub fn reverse(mut self) -> Self {
        std::mem::swap(&mut self.start, &mut self.end);
        self
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

    pub fn get_ghost_layer_modifyer(&self) -> (isize, isize) {
        if self.start[[0, 0]] == self.end[[0, 0]] {
            if self.start[[0, 0]] == 0 {
                (-1, 0)
            } else {
                (1, 0)
            }
        } else {
            if self.start[[1, 0]] == 0 {
                (0, -1)
            } else {
                (0, 1)
            }
        }
    }
}

pub struct BlockBoundaryRangeNewIter<'a> {
    range: &'a BlockBoundaryRange,
    dim: usize,
    step: isize,
    steps: usize,
    count: usize,
}

impl<'a> BlockBoundaryRangeNewIter<'a> {
    fn new(range: &'a BlockBoundaryRange) -> Self {
        let mut dim = 0;

        if range.start.first().unwrap() == range.end.first().unwrap() {
            dim = 1;
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

        if self.count <= self.steps + 1 {
            let mut index = self.range.start.clone();
            index[[self.dim, 0]] += self.step * (self.count as isize - 1);

            Some((index[[0, 0]] as usize, index[[1, 0]] as usize))
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct PeriodicBlockConnection {
    pub connection: BlockConnection,
    pub translation: Vec2d,
}

impl PeriodicBlockConnection {
    pub fn new(
        mesh: &Mesh,
        ranges: (BlockBoundaryRange, BlockBoundaryRange),
        translation: Vec2d,
    ) -> Self {
        let connection = BlockConnection::new_unchecked(mesh, ranges);
        connection.check_overlap(mesh, translation);
        Self {
            connection,
            translation,
        }
    }
}

#[derive(Debug)]
pub struct BlockConnection {
    pub donor: BlockBoundaryRange,
    pub receiver: BlockBoundaryRange,

    // 2x2 transformation matrix, see
    // https://cgns.github.io/CGNS_docs_current/sids/cnct.html#Transform
    transform: ndarray::Array2<isize>,
}

impl BlockConnection {
    pub fn new(mesh: &Mesh, ranges: (BlockBoundaryRange, BlockBoundaryRange)) -> Self {
        let connection = BlockConnection::new_unchecked(mesh, ranges);
        connection.check_overlap(mesh, Vec2d(0.0, 0.0));
        connection
    }

    pub fn new_unchecked(mesh: &Mesh, ranges: (BlockBoundaryRange, BlockBoundaryRange)) -> Self {
        let donor = ranges.0;
        let receiver = ranges.1;

        // create transformation matrix
        let mut transform = Array2::<isize>::zeros((2, 2));

        let donor_parallel_sign = if donor.is_increasing() { 1 } else { -1 };
        let rec_parallel_sign = if receiver.is_increasing() { 1 } else { -1 };

        let parallel_sign = donor_parallel_sign * rec_parallel_sign;

        let donor_edge_type = donor.edge_type();
        let rec_edge_type = receiver.edge_type();

        let donor_normal_sign = match donor_edge_type {
            EdgeIndex::IMin => -1,
            EdgeIndex::IMax => 1,
            EdgeIndex::JMin => -1,
            EdgeIndex::JMax => 1,
        };

        let rec_normal_sign = match rec_edge_type {
            EdgeIndex::IMin => 1,
            EdgeIndex::IMax => -1,
            EdgeIndex::JMin => 1,
            EdgeIndex::JMax => -1,
        };

        let normal_sign = donor_normal_sign * rec_normal_sign;

        if rec_edge_type == EdgeIndex::JMin || rec_edge_type == EdgeIndex::JMax {
            // parallel: t21 or t22
            if donor_edge_type == EdgeIndex::JMin || donor_edge_type == EdgeIndex::JMax {
                // parallel: t22
                transform[[1, 1]] = parallel_sign;
                // normal: t11
                transform[[0, 0]] = normal_sign;
            } else {
                // parallel: t21
                transform[[1, 0]] = parallel_sign;
                // normal: t12
                transform[[0, 1]] = normal_sign;
            }
        } else {
            // parallel: t11 or t12
            if donor_edge_type == EdgeIndex::JMin || donor_edge_type == EdgeIndex::JMax {
                // parallel: t12
                transform[[0, 1]] = parallel_sign;
                // normal: t21
                transform[[1, 0]] = normal_sign;
            } else {
                // parallel: t11
                transform[[0, 0]] = parallel_sign;
                // normal: t22
                transform[[1, 1]] = normal_sign;
            }
        }

        Self {
            donor,
            receiver,
            transform,
        }
    }

    fn check_overlap(&self, mesh: &Mesh, translation: Vec2d) {
        const EPSILON: f64 = 1e-10;

        let block_0 = &mesh.blocks[self.donor.block].coords;
        let block_1 = &mesh.blocks[self.receiver.block].coords;

        self.donor
            .iter()
            .zip(self.receiver.iter())
            .for_each(|(donor_idx, rec_idx)| {
                let x_0 = block_0[[donor_idx.0, donor_idx.1]];
                let x_1 = block_1[[rec_idx.0, rec_idx.1]] - translation;
                assert!(
                    approx_eq!(&Vec2d, &x_0, &x_1, epsilon = EPSILON),
                    "Non-matching Coordinates for Connection:\n
                            {:?} at index {:?} with {}\n
                            {:?} at index {:?} with {}",
                    self.donor,
                    donor_idx,
                    x_0,
                    self.receiver,
                    rec_idx,
                    x_1
                )
            });

        // check transform from donor to receiver w/ ghost layer on donor, i.e.
        // first internal layer on receiver

        let donor_edge_type = self.donor.edge_type();
        let rec_edge_type = self.receiver.edge_type();

        let add = match donor_edge_type {
            EdgeIndex::IMin => (0, -1),
            EdgeIndex::IMax => (0, 1),
            EdgeIndex::JMin => (-1, 0),
            EdgeIndex::JMax => (1, 0),
        };

        let add_rec: (isize, isize) = match rec_edge_type {
            EdgeIndex::IMin => (0, 1),
            EdgeIndex::IMax => (0, -1),
            EdgeIndex::JMin => (1, 0),
            EdgeIndex::JMax => (-1, 0),
        };

        self.donor
            .iter()
            .zip(self.receiver.iter())
            .for_each(|((i, j), (i_rec, j_rec))| {
                let i = i as isize + add.0;
                let j = j as isize + add.1;

                let ghost_point_comp = self.get_index_in_receiver_block(&[i, j]);

                let ghost_point = (i_rec as isize + add_rec.0, j_rec as isize + add_rec.1);

                assert_eq!(ghost_point_comp, ghost_point);
            });
    }

    pub fn get_index_in_receiver_block(&self, donor_index: &[isize; 2]) -> (isize, isize) {
        let donor_index = Array2::from_shape_vec((2, 1), donor_index.to_vec()).unwrap();
        let delta = donor_index - self.donor.start.clone();
        let receiver_index = self.transform.dot(&delta) + self.receiver.start.clone();
        (receiver_index[[0, 0]], receiver_index[[1, 0]])
    }
}
