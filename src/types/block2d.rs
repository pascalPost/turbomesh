// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::tfi::{boundary_blended_control_function, tfi_linear_2d};
use crate::types::{Array2d, Index, Scalar, Vec2d};

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
    // /// edge mapping from computational to physical space (x, y) = f(u in [0, 1])
    // pub edge_mappings: Vec<Box<dyn Fn(&[Scalar]) -> Box<[Vec2d]>>>,
}

#[derive(Copy, Clone)]
pub enum Edge {
    IMin,
    IMax,
    JMin,
    JMax,
}

pub trait EdgeFunction {
    fn i_min(&self, coords: &mut Array2d<Vec2d>);
    fn i_max(&self, coords: &mut Array2d<Vec2d>);
    fn j_min(&self, coords: &mut Array2d<Vec2d>);
    fn j_max(&self, coords: &mut Array2d<Vec2d>);
}

struct LinearEdge {}

impl LinearEdge {
    /// map intermediate control domain edge to physical
    /// assumes the corners to be well defined
    /// input: 0 <= rho <= 1
    fn linear_edge(rho: Scalar, x_0: Vec2d, x_1: Vec2d) -> Vec2d {
        x_0 + rho * (x_1 - x_0)
    }
}

impl EdgeFunction for LinearEdge {
    fn i_min(&self, coords: &mut Array2d<Vec2d>) {
        let [i_points, _] = coords.shape;

        let x_0 = coords[[0, 0]];
        let x_1 = coords[[i_points - 1, 0]];

        for i in 1..i_points - 1 {
            let rho = coords[[i, 0]].0;
            coords[[i, 0]] = Self::linear_edge(rho, x_0, x_1);
        }
    }

    fn i_max(&self, coords: &mut Array2d<Vec2d>) {
        let [i_points, j_points] = coords.shape;

        let x_0 = coords[[0, j_points - 1]];
        let x_1 = coords[[i_points - 1, j_points - 1]];

        for i in 1..i_points - 1 {
            let rho = coords[[i, j_points - 1]].0;
            coords[[i, j_points - 1]] = Self::linear_edge(rho, x_0, x_1);
        }
    }

    fn j_min(&self, coords: &mut Array2d<Vec2d>) {
        let [_, j_points] = coords.shape;

        let x_0 = coords[[0, 0]];
        let x_1 = coords[[0, j_points - 1]];

        for j in 1..j_points - 1 {
            let rho = coords[[0, j]].1;
            coords[[0, j]] = Self::linear_edge(rho, x_0, x_1);
        }
    }

    fn j_max(&self, coords: &mut Array2d<Vec2d>) {
        let [i_points, j_points] = coords.shape;

        let x_0 = coords[[i_points - 1, 0]];
        let x_1 = coords[[i_points - 1, j_points - 1]];

        for j in 1..j_points - 1 {
            let rho = coords[[i_points - 1, j]].1;
            coords[[i_points - 1, j]] = Self::linear_edge(rho, x_0, x_1);
        }
    }
}

fn equidistant_clustering(i: Index, n_points: Index) -> Scalar {
    i as Scalar / (n_points - 1) as Scalar
}

// pub struct EdgeSegmentProps<I: Iterator<Item = Index>> {
//     pub edge: Edge,
//     // pub range: Box<dyn Iterator<Item = Index>>,
//     pub range: I,
//     pub clustering: Box<dyn Fn(&[Index])>,
//     pub mapping: Box<dyn Fn(&[Scalar]) -> Box<[Vec2d]>>,
// }

pub struct EdgeSegmentProps {
    pub edge: Edge,
    pub start: Index, // TODO try inserting a type of range
    pub end: Index,
    pub clustering: Box<dyn Fn(&mut [Scalar])>,
    pub mapping: Box<dyn Fn(&[Scalar], &mut [Vec2d])>,
}

impl Block2d {
    /// new block given a name and the number of cells as the shape of the block.
    /// The number of cells in every direction needs to be greater than 0
    pub fn new(name: String, n_cells: (Index, Index), edges: Vec<EdgeSegmentProps>) -> Self {
        // check input
        assert!(
            n_cells.0 > 0 && n_cells.1 > 0,
            "size in every direction needs to be greater than 0."
        );

        // TODO check that all edges are properly defined
        // TODO check start end end indices for edges

        // allocate & initialize memory
        // initialize in intermediate control domain to give a sain default
        let mut block = Block2d {
            name,
            coords: Array2d::new([n_cells.0 + 1, n_cells.1 + 1]),
            // TODO save clusterings and mappings to struct
        };

        let [n_points_i, n_points_j] = block.coords.shape;

        // allocate 1d arrays representing the clustering in the intermediate
        // control domain with values 0<=u/v<=1
        let mut s1 = vec![Scalar::NAN; n_points_i]; // i_min edge
        let mut s2 = vec![Scalar::NAN; n_points_i]; // i_max edge

        let mut t1 = vec![Scalar::NAN; n_points_j]; // j_min edge
        let mut t2 = vec![Scalar::NAN; n_points_j]; // j_max edge

        // allocate 1d arrays representing the physical coordinates on the edges
        let mut coords_i_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points_i];
        let mut coords_i_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points_i];
        let mut coords_j_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points_j];
        let mut coords_j_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points_j];

        for edge in edges.iter() {
            match edge.edge {
                Edge::IMin => {
                    (edge.clustering)(&mut s1[edge.start..edge.end]);
                    (edge.mapping)(
                        &s1[edge.start..edge.end],
                        &mut coords_i_min[edge.start..edge.end],
                    );
                }
                Edge::IMax => {
                    (edge.clustering)(&mut s2[edge.start..edge.end]);
                    (edge.mapping)(
                        &s2[edge.start..edge.end],
                        &mut coords_i_max[edge.start..edge.end],
                    );
                }
                Edge::JMin => {
                    (edge.clustering)(&mut t1[edge.start..edge.end]);
                    (edge.mapping)(
                        &t1[edge.start..edge.end],
                        &mut coords_j_min[edge.start..edge.end],
                    );
                }
                Edge::JMax => {
                    (edge.clustering)(&mut t2[edge.start..edge.end]);
                    (edge.mapping)(
                        &t2[edge.start..edge.end],
                        &mut coords_j_max[edge.start..edge.end],
                    );
                }
            }
        }

        // check if all is not NAN any more to see if all values are set
        // TODO add better error message
        s1.iter().for_each(|u| assert!(!u.is_nan()));
        s2.iter().for_each(|u| assert!(!u.is_nan()));
        t1.iter().for_each(|u| assert!(!u.is_nan()));
        t2.iter().for_each(|u| assert!(!u.is_nan()));
        coords_i_min.iter().for_each(|Vec2d(x, y)| {
            assert!(!x.is_nan());
            assert!(!y.is_nan());
        });
        coords_i_max.iter().for_each(|Vec2d(x, y)| {
            assert!(!x.is_nan());
            assert!(!y.is_nan());
        });
        coords_j_min.iter().for_each(|Vec2d(x, y)| {
            assert!(!x.is_nan());
            assert!(!y.is_nan());
        });
        coords_j_max.iter().for_each(|Vec2d(x, y)| {
            assert!(!x.is_nan());
            assert!(!y.is_nan());
        });

        // TODO check that corners are 0 or 1 in computational space

        let [i_len, j_len] = block.coords.shape;

        // check equality of corners in physical space
        let x_0_0 = coords_j_min[0];
        assert_eq!(x_0_0, coords_i_min[0]);

        let x_1_0 = coords_i_min[i_len - 1];
        assert_eq!(x_1_0, coords_j_max[0]);

        let x_0_1 = coords_i_max[0];
        assert_eq!(x_0_1, coords_j_min[j_len - 1]);

        let x_1_1 = coords_j_max[j_len - 1];
        assert_eq!(x_1_1, coords_i_max[i_len - 1]);

        // TODO switch to iterator syntax to get independent from mem layout and
        // better suited for parallel exec

        for j in 0..j_len {
            let t1 = t1[j];
            let t2 = t2[j];

            let x_0_eta = coords_j_min[j];
            let x_1_eta = coords_j_max[j];

            for i in 0..i_len {
                let s1 = s1[i];
                let s2 = s2[i];

                let x_xi_0 = coords_i_min[i];
                let x_xi_1 = coords_i_max[i];

                let xi = ((1.0 - t1) * s1 + t1 * s2) / (1.0 - (s2 - s1) * (t2 - t1));
                let eta = ((1.0 - s1) * t1 + s1 * t2) / (1.0 - (t2 - t1) * (s2 - s1));

                let u_ij = (1.0 - xi) * x_0_eta + xi * x_1_eta;
                let v_ij = (1.0 - eta) * x_xi_0 + eta * x_xi_1;
                let uv_ij = xi * eta * x_1_1
                    + xi * (1.0 - eta) * x_1_0
                    + (1.0 - xi) * eta * x_0_1
                    + (1.0 - xi) * (1.0 - eta) * x_0_0;

                block.coords[[i, j]] = u_ij + v_ij - uv_ij
            }
        }

        block
    }

    /// returns the number of cells in every direction
    pub fn cells(&self) -> [Index; 2] {
        let [i_points, j_points] = self.points();
        [i_points - 1, j_points - 1]
    }

    pub fn points(&self) -> [Index; 2] {
        self.coords.shape
    }

    pub fn corners(mut self, c0: Vec2d, c1: Vec2d, c2: Vec2d, c3: Vec2d) -> Self {
        let [i_points, j_points] = self.coords.shape;
        self.coords[[0, 0]] = c0;
        self.coords[[i_points - 1, 0]] = c1;
        self.coords[[i_points - 1, j_points - 1]] = c2;
        self.coords[[0, j_points - 1]] = c3;

        self = self
            .mod_edge(Edge::IMin, LinearEdge {})
            .mod_edge(Edge::IMax, LinearEdge {})
            .mod_edge(Edge::JMin, LinearEdge {})
            .mod_edge(Edge::JMax, LinearEdge {});

        tfi_linear_2d(&mut self);
        self
    }

    pub fn mod_edge_control_domain<F>(&mut self, edge: Edge, start: Index, end: Index, f: F)
    where
        F: Fn(Index, Index) -> Scalar,
    {
        match edge {
            Edge::IMin => {
                let i_points = self.points()[0];
                for i in start..end {
                    self.coords[[i, 0]] = Vec2d(f(i, i_points), 0.0);
                }
            }
            Edge::IMax => {
                let [i_points, j_points] = self.points();
                for i in start..end {
                    self.coords[[i, j_points - 1]] = Vec2d(f(i, i_points), 1.0);
                }
            }
            Edge::JMin => {
                let j_points = self.points()[1];
                for j in start..end {
                    self.coords[[0, j]] = Vec2d(0.0, f(j, j_points));
                }
            }
            Edge::JMax => {
                let [i_points, j_points] = self.points();
                for j in start..end {
                    self.coords[[i_points - 1, j]] = Vec2d(1.0, f(j, j_points));
                }
            }
        }
    }

    pub fn mod_edge(mut self, edge: Edge, f: impl EdgeFunction) -> Self {
        match edge {
            Edge::IMin => f.i_min(&mut self.coords),
            Edge::IMax => f.i_max(&mut self.coords),
            Edge::JMin => f.j_min(&mut self.coords),
            Edge::JMax => f.j_max(&mut self.coords),
        }

        self
    }
}
