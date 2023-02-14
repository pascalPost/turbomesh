// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::tfi::tfi_linear_2d;
use crate::types::{Array2d, Index, Scalar, Vec2d};

use super::segment::SegmentFunction;

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
    pub coords: Array2d<Vec2d>,
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
        let mut x = Array2d::<Vec2d>::new([n_points.0, n_points.1]);
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

    fn apply_clustering_and_mapping(
        edge: &Vec<Box<dyn SegmentFunction>>,
        u: &mut Vec<Scalar>,
        x: &mut Vec<Vec2d>,
    ) {
        let mut start = 0;
        edge.iter().for_each(|segment| {
            let end = start + segment.len();
            segment.apply_clustering(&mut u[start..end]);
            start = end - 1;
        });

        // renormalization to 0<=u<=1
        let u_max = *u.last().unwrap();
        u.iter_mut().for_each(|u| *u /= u_max);

        let mut start = 0;
        edge.iter().for_each(|segment| {
            let end = start + segment.len();
            segment.computational_to_physical(&u[start..end], &mut x[start..end]);
            start = end - 1;
        });
    }

    /// returns the number of cells in every direction
    pub fn cells(&self) -> [Index; 2] {
        let [i_points, j_points] = self.points();
        [i_points - 1, j_points - 1]
    }

    pub fn points(&self) -> [Index; 2] {
        self.coords.shape
    }
}
