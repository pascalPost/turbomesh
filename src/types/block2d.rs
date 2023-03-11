// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use super::segment::SegmentFunction;
use crate::tfi::tfi_linear_2d;
use crate::types::{BlockEdgeData, Index, Scalar, Vec2d};
use ndarray::Array2;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeIndex {
    IMin,
    IMax,
    JMin,
    JMax,
}

/// block corner and edge order:
///        2
///   3 ------- 2
///   |         |
/// 3 |         | 1
///   |         |
///   0 ------- 1
///        0
// #[derive(Debug)]
pub struct Block2d {
    pub name: String,
    pub coords: Array2<Vec2d>,

    // TODO change the edge data to array of size 4
    pub edge_i_min: Vec<Box<dyn SegmentFunction>>,
    pub edge_i_max: Vec<Box<dyn SegmentFunction>>,
    pub edge_j_min: Vec<Box<dyn SegmentFunction>>,
    pub edge_j_max: Vec<Box<dyn SegmentFunction>>,
}

impl Block2d {
    /// new block given a name and the number of cells as the shape of the block.
    /// The number of cells in every direction needs to be greater than 0
    pub fn new(
        name: String,
        // n_cells: (Index, Index),
        edge_i_min: Vec<Box<dyn SegmentFunction>>,
        edge_i_max: Vec<Box<dyn SegmentFunction>>,
        edge_j_min: Vec<Box<dyn SegmentFunction>>,
        edge_j_max: Vec<Box<dyn SegmentFunction>>,
    ) -> Self {
        // compute block size from edges (accounting for overlapping points)

        let n_points_i_min =
            edge_i_min.iter().map(|seg| seg.len()).sum::<usize>() - (edge_i_min.len() - 1);
        let n_points_i_max =
            edge_i_max.iter().map(|seg| seg.len()).sum::<usize>() - (edge_i_max.len() - 1);
        assert!(
            n_points_i_min == n_points_i_max,
            "Number of points of edge i min ({}) and edge i max ({}) must be equal.",
            n_points_i_min,
            n_points_i_max
        );

        let n_points_j_min =
            edge_j_min.iter().map(|seg| seg.len()).sum::<usize>() - (edge_j_min.len() - 1);
        let n_points_j_max =
            edge_j_max.iter().map(|seg| seg.len()).sum::<usize>() - (edge_j_max.len() - 1);
        assert!(
            n_points_j_min == n_points_j_max,
            "Number of points of edge j min ({}) and edge j max ({}) must be equal.",
            n_points_j_min,
            n_points_j_max
        );

        let n_points = (n_points_i_min, n_points_j_min);

        // allocate 1d arrays representing the clustering in the intermediate
        // control domain with values 0<=u/v<=1
        let mut s1 = vec![Scalar::NAN; n_points.0]; // i_min edge
        let mut s2 = vec![Scalar::NAN; n_points.0]; // i_max edge
        let mut t1 = vec![Scalar::NAN; n_points.1]; // j_min edge
        let mut t2 = vec![Scalar::NAN; n_points.1]; // j_max edge

        // allocate 1d arrays representing the physical coordinates on the edges
        let mut x_i_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.0];
        let mut x_i_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.0];
        let mut x_j_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.1];
        let mut x_j_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.1];

        Self::apply_clustering_and_mapping(&edge_i_min, &mut s1, &mut x_i_min);
        Self::apply_clustering_and_mapping(&edge_i_max, &mut s2, &mut x_i_max);
        Self::apply_clustering_and_mapping(&edge_j_min, &mut t1, &mut x_j_min);
        Self::apply_clustering_and_mapping(&edge_j_max, &mut t2, &mut x_j_max);

        // allocate 2d array and run tfi on it
        let mut x = Array2::<Vec2d>::zeros((n_points.0, n_points.1));
        tfi_linear_2d(
            &s1, &s2, &t1, &t2, &x_i_min, &x_i_max, &x_j_min, &x_j_max, &mut x,
        );

        Block2d {
            name,
            coords: x,
            edge_i_min,
            edge_i_max,
            edge_j_min,
            edge_j_max,
        }
    }

    pub fn segments(&self, edge: EdgeIndex) -> &Vec<Box<dyn SegmentFunction>> {
        match edge {
            EdgeIndex::IMin => &self.edge_i_min,
            EdgeIndex::IMax => &self.edge_i_max,
            EdgeIndex::JMin => &self.edge_j_min,
            EdgeIndex::JMax => &self.edge_j_max,
        }
    }

    /// returns the discrete edge data
    pub fn edge_data(&self, edge: EdgeIndex) -> BlockEdgeData {
        let n_points = self.points();
        match edge {
            EdgeIndex::IMin => {
                let mut s1 = vec![Scalar::NAN; n_points[0]];
                let mut x_i_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points[0]];
                Self::apply_clustering_and_mapping(&self.edge_i_min, &mut s1, &mut x_i_min);
                BlockEdgeData::new(s1, x_i_min)
            }
            EdgeIndex::IMax => {
                let mut s2 = vec![Scalar::NAN; n_points[0]];
                let mut x_i_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points[0]];
                Self::apply_clustering_and_mapping(&self.edge_i_max, &mut s2, &mut x_i_max);
                BlockEdgeData::new(s2, x_i_max)
            }
            EdgeIndex::JMin => {
                let mut t1 = vec![Scalar::NAN; n_points[1]];
                let mut x_j_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points[1]];
                Self::apply_clustering_and_mapping(&self.edge_j_min, &mut t1, &mut x_j_min);
                BlockEdgeData::new(t1, x_j_min)
            }
            EdgeIndex::JMax => {
                let mut t2 = vec![Scalar::NAN; n_points[1]];
                let mut x_j_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points[1]];
                Self::apply_clustering_and_mapping(&self.edge_j_max, &mut t2, &mut x_j_max);
                BlockEdgeData::new(t2, x_j_max)
            }
        }
    }

    pub fn edge_segment(&self, edge: EdgeIndex, segment_index: usize) -> BlockEdgeData {
        let segment = &self.segments(edge)[segment_index];

        let n_points = segment.len();
        let mut u = vec![Scalar::NAN; n_points];
        let mut x = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points];
        Self::apply_clustering_and_mapping(std::slice::from_ref(segment), &mut u, &mut x);
        BlockEdgeData::new(u, x)
    }

    pub fn edge_segments<Idx>(
        &self,
        edge: EdgeIndex,
        segments: impl std::slice::SliceIndex<
            [std::boxed::Box<(dyn SegmentFunction + 'static)>],
            Output = [Box<dyn SegmentFunction>],
        >,
    ) -> BlockEdgeData {
        let segments = &self.segments(edge)[segments];

        let n_points = segments.iter().map(|seg| seg.len()).sum();

        let mut u = vec![Scalar::NAN; n_points];
        let mut x = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points];

        Self::apply_clustering_and_mapping(segments, &mut u, &mut x);
        BlockEdgeData::new(u, x)
    }

    /// apply segment clustering and mapping to the given edga data
    fn apply_clustering_and_mapping(
        edge: &[Box<dyn SegmentFunction>],
        u: &mut Vec<Scalar>,
        x: &mut Vec<Vec2d>,
    ) {
        let mut start = 0;
        edge.iter().for_each(|segment| {
            let end = start + segment.len();
            segment.add_segment_to_edge(&mut u[start..end], &mut x[start..end]);
            start = end - 1;
        });
    }

    /// returns the number of cells in every direction
    pub fn cells(&self) -> [Index; 2] {
        let [i_points, j_points] = self.points();
        [i_points - 1, j_points - 1]
    }

    /// returns the numnber of points in every direction
    pub fn points(&self) -> [Index; 2] {
        let shape = self.coords.shape();
        [shape[0], shape[1]]
    }
}
