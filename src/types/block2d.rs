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
#[derive(Debug)]
pub struct Block2d {
    pub name: String,
    pub coords: Array2d<Vec2d>,
}

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

impl Block2d {
    /// new block given a name and the number of cells as the shape of the block.
    /// The number of cells in every direction needs to be greater than 0
    pub fn new(name: String, n_cells: (Index, Index)) -> Self {
        // check input
        assert!(
            n_cells.0 > 0 && n_cells.1 > 0,
            "size in every direction needs to be greater than 0."
        );

        // allocate & initialize memory
        // initialize in intermediate control domain to give a sain default
        let mut block = Block2d {
            name,
            coords: Array2d::new([n_cells.0 + 1, n_cells.1 + 1]),
        }
        .mod_edge_control_domain(Edge::IMin, equidistant_clustering)
        .mod_edge_control_domain(Edge::IMax, equidistant_clustering)
        .mod_edge_control_domain(Edge::JMin, equidistant_clustering)
        .mod_edge_control_domain(Edge::JMax, equidistant_clustering);

        boundary_blended_control_function(&mut block.coords);

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

    pub fn mod_edge_control_domain<F>(mut self, edge: Edge, f: F) -> Self
    where
        F: Fn(Index, Index) -> Scalar,
    {
        match edge {
            Edge::IMin => {
                let i_points = self.points()[0];
                for i in 0..i_points {
                    self.coords[[i, 0]] = Vec2d(f(i, i_points), 0.0);
                }
            }
            Edge::IMax => {
                let [i_points, j_points] = self.points();
                for i in 0..i_points {
                    self.coords[[i, j_points - 1]] = Vec2d(f(i, i_points), 1.0);
                }
            }
            Edge::JMin => {
                let j_points = self.points()[1];
                for j in 0..j_points {
                    self.coords[[0, j]] = Vec2d(0.0, f(j, j_points));
                }
            }
            Edge::JMax => {
                let [i_points, j_points] = self.points();
                for j in 0..j_points {
                    self.coords[[i_points - 1, j]] = Vec2d(1.0, f(j, j_points));
                }
            }
        }

        self
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
