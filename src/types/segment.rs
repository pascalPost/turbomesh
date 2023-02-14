// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Scalar, Vec2d};

/// trait representing the function that performs the clustering on an edge segment
pub trait ClusteringFunction {
    fn apply_clustering(&self, u: &mut [Scalar]);
}

/// trait representing the function to map the computational 0<=u<=1 to the
/// physical space
pub trait MappingFunction {
    fn computational_to_physical(&self, u: &[Scalar], x: &mut [Vec2d]);

    fn computational_to_physical_vec(&self, u: &[Scalar]) -> Vec<Vec2d> {
        let mut x = vec![Vec2d(0.0, 0.0); u.len()];
        self.computational_to_physical(u, x.as_mut_slice());
        x
    }

    fn computational_to_physical_val(&self, u: Scalar) -> Vec2d {
        let mut x = [Vec2d(0.0, 0.0); 1];
        self.computational_to_physical(&[u], &mut x);
        x[0]
    }
}

pub trait SegmentFunction: ClusteringFunction + MappingFunction {
    fn len(&self) -> usize;
}

/// representing needed edge properties for tfi
pub struct Segment<C: ClusteringFunction, M: MappingFunction> {
    points: usize,
    pub clustering: C,
    pub mapping: M,
}

impl<C: ClusteringFunction, M: MappingFunction> Segment<C, M> {
    pub fn new(points: usize, clustering: C, mapping: M) -> Self {
        Self {
            points,
            clustering,
            mapping,
        }
    }

    pub fn apply_clustering(&self, u: &mut [Scalar]) {
        self.clustering.apply_clustering(u);
    }

    pub fn apply_mapping(&self, u: &[Scalar], x: &mut [Vec2d]) {
        self.mapping.computational_to_physical(&u, x);
    }
}

impl<C: ClusteringFunction, M: MappingFunction> ClusteringFunction for Segment<C, M> {
    fn apply_clustering(&self, u: &mut [Scalar]) {
        self.clustering.apply_clustering(u);
    }
}

impl<C: ClusteringFunction, M: MappingFunction> MappingFunction for Segment<C, M> {
    fn computational_to_physical(&self, u: &[Scalar], x: &mut [Vec2d]) {
        self.mapping.computational_to_physical(u, x);
    }
}

impl<C: ClusteringFunction, M: MappingFunction> SegmentFunction for Segment<C, M> {
    fn len(&self) -> usize {
        self.points
    }
}
