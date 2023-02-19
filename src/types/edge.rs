// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Index, Scalar};
use crate::{Segment, Vec2d};
use std::rc::Rc;

use super::segment::SegmentFunction;
use super::{ClusteringFunction, MappingFunction};

/// descritzed edge (1D block)
///   0 ------- 1
#[derive(Debug)]
pub struct Edge {
    pub name: String,

    /// clustering points with u in [0,1]
    pub u: Vec<Scalar>,

    /// physical coordinates corresponding to clustering u
    pub coords: Vec<Vec2d>,
}

impl Edge {
    /// new edge given a name, the number of cells on the edge and the
    pub fn new<C: ClusteringFunction, M: MappingFunction>(
        name: String,
        segments: &Vec<Segment<C, M>>,
    ) -> Self {
        let n_points = segments.iter().map(|seg| seg.len()).sum();

        // allocate 1d array representing the clustering in the intermediate
        // control domain with values 0<=u/v<=1
        let mut u = vec![Scalar::NAN; n_points]; // i_min edge

        // allocate 1d arrays representing the physical coordinates on the edges
        let mut x = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points];

        // apply input to the clustering and mappings
        segments.iter().for_each(|segment| {
            segment.apply_clustering(&mut u);
            segment.apply_mapping(&u, &mut x);
        });

        Edge { name, coords: x, u }
    }

    /// returns the number of cells
    pub fn cells(&self) -> Index {
        self.points() + 1
    }

    /// returns the number of points
    pub fn points(&self) -> Index {
        self.coords.len()
    }
}

/// represents a view into an edge
#[derive(Debug)]
pub struct EdgeView {
    edge: Rc<Edge>,
    pub start: Index,
    pub end: Index,
}

impl EdgeView {
    pub fn new(edge: Edge) -> Self {
        let end = edge.coords.len();
        Self {
            edge: Rc::new(edge),
            start: 0,
            end,
        }
    }

    pub fn edge(&self) -> &Rc<Edge> {
        &self.edge
    }

    pub fn len(&self) -> Index {
        if self.start > self.end {
            self.start - self.end + 1
        } else {
            self.end - self.start + 1
        }
    }

    pub fn split_at(&self, point: Index) -> (Self, Self) {
        assert!(
            point > 0 && self.start + point < self.end,
            "Point must be with in range ({},{})",
            0,
            self.len() - 1
        );

        (
            Self {
                edge: self.edge.clone(),
                start: self.start,
                end: self.start + point,
            },
            Self {
                edge: self.edge.clone(),
                start: self.start + point,
                end: self.end,
            },
        )
    }

    /// return a reversed version of the given EdgeView
    pub fn rev(&self) -> Self {
        Self {
            edge: self.edge.clone(),
            start: self.end,
            end: self.start,
        }
    }

    pub fn point_coord(&self, index: usize) -> Vec2d {
        if self.start > self.end {
            assert!(
                index >= self.end && index <= self.start,
                "Given index {} is not in [{},{}]",
                index,
                self.end,
                self.start
            );
        } else {
            assert!(
                index >= self.start && index <= self.end,
                "Given index {} is not in [{},{}]",
                index,
                self.start,
                self.end
            );
        }

        self.edge.coords[index]
    }
}

impl ClusteringFunction for EdgeView {
    fn apply_clustering(&self, u: &mut [Scalar]) {
        // check if u needs to be initialized (nan) or if it contains value from
        // other segment
        if u[0].is_nan() {
            u[0] = 0.0;
        }

        if self.start > self.end {
            self.edge.u[self.end..=self.start]
                .windows(2)
                .rev()
                .zip(u[1..].iter_mut())
                .for_each(|(u_ref, u)| *u = u_ref[1] - u_ref[0]);
        } else {
            self.edge.u[self.start..=self.end]
                .windows(2)
                .zip(u[1..].iter_mut())
                .for_each(|(u_ref, u)| *u = u_ref[1] - u_ref[0]);
        };

        // cum sum
        let u_max = u.iter_mut().fold(0.0, |acc, x| {
            *x += acc;
            *x
        });
    }
}

impl MappingFunction for EdgeView {
    fn computational_to_physical(&self, _u: &[Scalar], x: &mut [Vec2d]) {
        assert!(
            self.len() == x.len(),
            "EdgeView of size {} cannot be applied to Edge of size {}",
            self.edge.coords.len(),
            x.len()
        );

        if self.start > self.end {
            self.edge.coords[self.end..=self.start]
                .iter()
                .rev()
                .zip(x.iter_mut())
                .for_each(|(x_ref, x)| *x = *x_ref);
        } else {
            self.edge.coords[self.start..=self.end]
                .iter()
                .zip(x.iter_mut())
                .for_each(|(x_ref, x)| *x = *x_ref);
        };
    }
}

impl SegmentFunction for EdgeView {
    fn len(&self) -> usize {
        self.len()
    }
}

pub struct BlockEdgeData {
    pub u: Vec<Scalar>,
    pub x: Vec<Vec2d>,
}

impl BlockEdgeData {
    pub fn new(u: Vec<Scalar>, x: Vec<Vec2d>) -> Self {
        Self { u, x }
    }

    pub fn len(&self) -> usize {
        assert!(self.u.len() == self.x.len());
        self.u.len()
    }
}

impl ClusteringFunction for BlockEdgeData {
    fn apply_clustering(&self, u: &mut [Scalar]) {
        assert!(
            u.len() == self.u.len(),
            "BlockEdgeData of length {} cannot be applied to edge of length {}",
            self.len(),
            u.len()
        );

        self.u
            .iter()
            .zip(u.iter_mut())
            .for_each(|(u_ref, u)| *u = *u_ref);
    }
}

impl MappingFunction for BlockEdgeData {
    fn computational_to_physical(&self, _u: &[Scalar], x: &mut [Vec2d]) {
        assert!(
            self.len() == x.len(),
            "BlockEdgeData of length {} cannot be applied to edge of length {}",
            self.len(),
            x.len()
        );

        self.x
            .iter()
            .zip(x.iter_mut())
            .for_each(|(x_ref, x)| *x = *x_ref);
    }
}

impl SegmentFunction for BlockEdgeData {
    fn len(&self) -> usize {
        self.len()
    }
}
