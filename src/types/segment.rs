// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Scalar, Vec2d};

/// trait representing the function that performs the clustering on an edge segment
pub trait ClusteringFunction {
    fn get_clustering(&self, points: usize) -> Vec<Scalar>;
}

pub trait DiscretizableCurve {
    fn discretize_curve(
        &self,
        clustering: &impl ClusteringFunction,
        u: &mut [Scalar],
        x: &mut [Vec2d],
    );

    fn arclength(&self) -> Scalar;
}

// TODO rename to EdgeSegment
pub trait SegmentFunction {
    fn add_segment_to_edge(&self, u: &mut [Scalar], x: &mut [Vec2d]);
    fn len(&self) -> usize;
}

/// representing needed edge properties for tfi
pub struct Segment<CF: ClusteringFunction, DC: DiscretizableCurve> {
    points: usize,
    pub clustering: CF,
    pub curve: DC,
}

impl<CF: ClusteringFunction, DC: DiscretizableCurve> Segment<CF, DC> {
    pub fn new(points: usize, clustering: CF, curve: DC) -> Self {
        Self {
            points,
            clustering,
            curve,
        }
    }
}

impl<CF: ClusteringFunction, DC: DiscretizableCurve> SegmentFunction for Segment<CF, DC> {
    fn add_segment_to_edge(&self, u: &mut [Scalar], x: &mut [Vec2d]) {
        self.curve.discretize_curve(&self.clustering, u, x);
    }

    fn len(&self) -> usize {
        self.points
    }
}
