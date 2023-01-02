// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::tfi::tfi_linear_2d;
use crate::types::{Array2d, Index, Scalar, Vec2d};

// trait Clustering {

// }

/// trait representing the function to map the computational 0<=u<=1 to the
/// physical space
pub trait Mapping {
    fn computational_to_physical(&self, u: &[Scalar], x: &mut [Vec2d]);
}

/// representing needed edge properties for tfi
pub struct Edge {
    pub clustering: Box<dyn Fn(&mut [Scalar])>,
    pub mapping: Box<dyn Mapping>,
}

impl Edge {
    /// creates a new edge. Designed to not bother with the Box as an
    /// implementation detail. Might be changed in the future if needed.
    pub fn new<F: Fn(&mut [Scalar]) + 'static, G: Mapping + 'static>(
        clustering: F,
        mapping: G,
    ) -> Self {
        Edge {
            clustering: Box::new(clustering),
            mapping: Box::new(mapping),
        }
    }
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
    pub coords: Array2d<Vec2d>,
    pub edge_i_min: Edge,
    pub edge_i_max: Edge,
    pub edge_j_min: Edge,
    pub edge_j_max: Edge,
}

// #[derive(Copy, Clone)]
// pub enum Edge {
//     IMin,
//     IMax,
//     JMin,
//     JMax,
// }

// pub struct EdgeSegmentProps {
//     pub edge: Edge,
//     pub clustering: Box<dyn Fn(&mut [Scalar])>,
//     pub mapping: Box<dyn Fn(&[Scalar], &mut [Vec2d])>,
// }

impl Block2d {
    /// new block given a name and the number of cells as the shape of the block.
    /// The number of cells in every direction needs to be greater than 0
    pub fn new(
        name: String,
        n_cells: (Index, Index),
        edge_i_min: Edge,
        edge_i_max: Edge,
        edge_j_min: Edge,
        edge_j_max: Edge,
    ) -> Self {
        // check input
        assert!(
            n_cells.0 > 0 && n_cells.1 > 0,
            "size in every direction needs to be greater than 0."
        );

        let n_points = (n_cells.0 + 1, n_cells.1 + 1);

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

        // apply input to the clustering and mappings
        (edge_i_min.clustering)(&mut s1);
        edge_i_min
            .mapping
            .computational_to_physical(&s1, &mut x_i_min);
        (edge_i_max.clustering)(&mut s2);
        edge_i_max
            .mapping
            .computational_to_physical(&s1, &mut x_i_max);
        (edge_j_min.clustering)(&mut t1);
        edge_j_min
            .mapping
            .computational_to_physical(&s1, &mut x_j_min);
        (edge_j_max.clustering)(&mut t2);
        edge_j_max
            .mapping
            .computational_to_physical(&s1, &mut x_j_max);

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

    // /// new block given a name and the number of cells as the shape of the block.
    // /// The number of cells in every direction needs to be greater than 0
    // pub fn new(name: String, n_cells: (Index, Index), edges: [EdgeSegmentProps; 4]) -> Self {
    //     // check input
    //     assert!(
    //         n_cells.0 > 0 && n_cells.1 > 0,
    //         "size in every direction needs to be greater than 0."
    //     );

    //     let n_points = (n_cells.0 + 1, n_cells.1 + 1);

    //     // allocate 1d arrays representing the clustering in the intermediate
    //     // control domain with values 0<=u/v<=1
    //     let mut s1 = vec![Scalar::NAN; n_points.0]; // i_min edge
    //     let mut s2 = vec![Scalar::NAN; n_points.0]; // i_max edge
    //     let mut t1 = vec![Scalar::NAN; n_points.1]; // j_min edge
    //     let mut t2 = vec![Scalar::NAN; n_points.1]; // j_max edge

    //     // allocate 1d arrays representing the physical coordinates on the edges
    //     let mut x_i_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.0];
    //     let mut x_i_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.0];
    //     let mut x_j_min = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.1];
    //     let mut x_j_max = vec![Vec2d(Scalar::NAN, Scalar::NAN); n_points.1];

    //     // apply input to the clustering and mappings
    //     for edge in edges.iter() {
    //         match edge.edge {
    //             Edge::IMin => {
    //                 (edge.clustering)(&mut s1);
    //                 (edge.mapping)(&s1, &mut x_i_min);
    //             }
    //             Edge::IMax => {
    //                 (edge.clustering)(&mut s2);
    //                 (edge.mapping)(&s2, &mut x_i_max);
    //             }
    //             Edge::JMin => {
    //                 (edge.clustering)(&mut t1);
    //                 (edge.mapping)(&t1, &mut x_j_min);
    //             }
    //             Edge::JMax => {
    //                 (edge.clustering)(&mut t2);
    //                 (edge.mapping)(&t2, &mut x_j_max);
    //             }
    //         }
    //     }

    //     let mut x = Array2d::<Vec2d>::new([n_points.0, n_points.1]);

    //     tfi_linear_2d(
    //         &s1, &s2, &t1, &t2, &x_i_min, &x_i_max, &x_j_min, &x_j_max, &mut x,
    //     );

    //     // allocate & initialize memory
    //     // initialize in intermediate control domain to give a sain default
    //     Block2d {
    //         name,
    //         coords: x,
    //         // TODO save clusterings and mappings to struct
    //     }
    // }

    /// returns the number of cells in every direction
    pub fn cells(&self) -> [Index; 2] {
        let [i_points, j_points] = self.points();
        [i_points - 1, j_points - 1]
    }

    pub fn points(&self) -> [Index; 2] {
        self.coords.shape
    }
}
