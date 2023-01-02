// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::block2d::{Edge, Mapping};
use crate::types::{Block2d, Scalar, Vec2d};

fn equidistant_clustering(u: &mut [Scalar]) {
    let n = u.len();
    for (i, u) in u.iter_mut().enumerate() {
        *u = i as Scalar / (n - 1) as Scalar
    }
}

/// representing a straigt line from x_min to x_max
struct Line2d {
    x_min: Vec2d,
    x_max: Vec2d,
}

impl Line2d {
    /// returns a Line2d from x0 to x1
    fn new(x0: Vec2d, x1: Vec2d) -> Self {
        Line2d {
            x_min: x0,
            x_max: x1,
        }
    }
}

impl Mapping for Line2d {
    fn computational_to_physical(&self, u: &[Scalar], x: &mut [Vec2d]) {
        let dx = self.x_max - self.x_min;
        for (&u, x) in u.iter().zip(x.iter_mut()) {
            *x = self.x_min + u * dx;
        }
    }
}

struct Rectangle {
    origin: Vec2d,
    size: Vec2d,
}

impl Rectangle {
    fn new(origin: Vec2d, size: Vec2d) -> Self {
        assert!(
            size.0 > 0.0 && size.1 > 0.0,
            "Size must be greater than zero"
        );
        Rectangle { origin, size }
    }

    fn i_min(&self) -> Line2d {
        Line2d::new(self.origin, self.origin + Vec2d(self.size.0, 0.0))
    }

    fn i_max(&self) -> Line2d {
        let x0 = self.origin + Vec2d(0.0, self.size.1);
        Line2d::new(x0, x0 + Vec2d(self.size.0, 0.0))
    }

    fn j_min(&self) -> Line2d {
        Line2d::new(self.origin, self.origin + Vec2d(0.0, self.size.1))
    }

    fn j_max(&self) -> Line2d {
        let x0 = self.origin + Vec2d(self.size.0, 0.0);
        Line2d::new(x0, x0 + Vec2d(0.0, self.size.1))
    }
}

/// test the block creation with a simple rectangle
#[test]
fn test_default_block_creation() {
    let rect = Rectangle::new(Vec2d(0.0, 0.0), Vec2d(1.0, 1.0));

    let block = Block2d::new(
        String::from("block"),
        (2, 2),
        Edge::new(equidistant_clustering, rect.i_min()),
        Edge::new(equidistant_clustering, rect.i_max()),
        Edge::new(equidistant_clustering, rect.j_min()),
        Edge::new(equidistant_clustering, rect.j_max()),
    );

    println!("{:?}", block.coords);

    // test block name
    assert_eq!(block.name, String::from("block"));

    // test contained number of points in i and j directions
    assert_eq!(block.coords.shape, [3, 3]);

    // test total number of contained points
    assert_eq!(block.coords.size(), 9);

    // test defaultly set block coordinates
    assert_eq!(block.coords[[0, 0]], Vec2d(0., 0.));
    assert_eq!(block.coords[[1, 0]], Vec2d(0.5, 0.));
    assert_eq!(block.coords[[2, 0]], Vec2d(1.0, 0.));

    assert_eq!(block.coords[[0, 1]], Vec2d(0., 0.5));
    assert_eq!(block.coords[[1, 1]], Vec2d(0.5, 0.5));
    assert_eq!(block.coords[[2, 1]], Vec2d(1.0, 0.5));

    assert_eq!(block.coords[[0, 2]], Vec2d(0., 1.0));
    assert_eq!(block.coords[[1, 2]], Vec2d(0.5, 1.0));
    assert_eq!(block.coords[[2, 2]], Vec2d(1.0, 1.0));
}

// /// test the block creation with corners (and linear, equidistant edges)
// #[test]
// fn test_block_creation_with_corners() {
//     let block = Block2d::new(
//         String::from("block"),
//         (2, 2),
//         vec![
//             Box::new(|u| {
//                 Vec::from_iter(
//                     u.iter()
//                         .map(|&u| Vec2d(-1.0, -1.0) + u * (Vec2d(1.0, -1.0) - Vec2d(-1.0, -1.0))),
//                 )
//                 .into_boxed_slice()
//             }),
//             Box::new(|u| {
//                 Vec::from_iter(
//                     u.iter()
//                         .map(|&u| Vec2d(1.0, -1.0) + u * (Vec2d(1.0, 1.0) - Vec2d(1.0, -1.0))),
//                 )
//                 .into_boxed_slice()
//             }),
//             Box::new(|u| {
//                 Vec::from_iter(
//                     u.iter()
//                         .map(|&u| Vec2d(1.0, 1.0) + u * (Vec2d(-1.0, 1.0) - Vec2d(1.0, 1.0))),
//                 )
//                 .into_boxed_slice()
//             }),
//             Box::new(|u| {
//                 Vec::from_iter(
//                     u.iter()
//                         .map(|&u| Vec2d(-1.0, 1.0) + u * (Vec2d(-1.0, -1.0) - Vec2d(-1.0, 1.0))),
//                 )
//                 .into_boxed_slice()
//             }),
//         ],
//     );
//     // .corners(
//     //     Vec2d(-1.0, -1.0),
//     //     Vec2d(1.0, -1.0),
//     //     Vec2d(1.0, 1.0),
//     //     Vec2d(-1.0, 1.0),
//     // );

//     // test block coordinates
//     assert_eq!(block.coords[[0, 0]], Vec2d(-1.0, -1.0));
//     assert_eq!(block.coords[[1, 0]], Vec2d(0.0, -1.0));
//     assert_eq!(block.coords[[2, 0]], Vec2d(1.0, -1.0));

//     assert_eq!(block.coords[[0, 1]], Vec2d(-1.0, 0.0));
//     assert_eq!(block.coords[[1, 1]], Vec2d(0.0, 0.0));
//     assert_eq!(block.coords[[2, 1]], Vec2d(1.0, 0.0));

//     assert_eq!(block.coords[[0, 2]], Vec2d(-1.0, 1.0));
//     assert_eq!(block.coords[[1, 2]], Vec2d(0.0, 1.0));
//     assert_eq!(block.coords[[2, 2]], Vec2d(1.0, 1.0));
// }

// /// test linear block creation
// ///           x 2
// ///         /   \
// ///       x 1    x 3
// ///         \   /
// ///           x 0
// #[test]
// fn test_block_creation_linear() {
//     let block = Block2d::new(String::from("block"), (2, 2)).corners(
//         Vec2d(0.0, -1.0),
//         Vec2d(-1.0, 0.0),
//         Vec2d(0.0, 1.0),
//         Vec2d(1.0, 0.0),
//     );

//     // test block coordinates
//     assert_eq!(block.coords[[0, 0]], Vec2d(0.0, -1.0));
//     assert_eq!(block.coords[[1, 0]], Vec2d(-0.5, -0.5));
//     assert_eq!(block.coords[[2, 0]], Vec2d(-1.0, 0.0));

//     assert_eq!(block.coords[[0, 1]], Vec2d(0.5, -0.5));
//     assert_eq!(block.coords[[1, 1]], Vec2d(0.0, 0.0));
//     assert_eq!(block.coords[[2, 1]], Vec2d(-0.5, 0.5));

//     assert_eq!(block.coords[[0, 2]], Vec2d(1.0, 0.0));
//     assert_eq!(block.coords[[1, 2]], Vec2d(0.5, 0.5));
//     assert_eq!(block.coords[[2, 2]], Vec2d(0.0, 1.0));
// }

// /// test creation of block with curved edges
// ///           x 2
// ///         /   \
// ///       x 1    x 3
// ///         \   /
// ///           x 0
// #[test]
// fn test_block_creation_curved() {
//     let linear_map = ||

//     let block = Block2d::new(String::from("block"), (2, 2), ).

//     // .corners(
//     //     Vec2d(0.0, -1.0),
//     //     Vec2d(-1.0, 0.0),
//     //     Vec2d(0.0, 1.0),
//     //     Vec2d(1.0, 0.0),
//     // );

//     // // test block coordinates
//     // assert_eq!(block.coords[[0, 0]], Vec2d(0.0, -1.0));
//     // assert_eq!(block.coords[[1, 0]], Vec2d(-0.5, -0.5));
//     // assert_eq!(block.coords[[2, 0]], Vec2d(-1.0, 0.0));

//     // assert_eq!(block.coords[[0, 1]], Vec2d(0.5, -0.5));
//     // assert_eq!(block.coords[[1, 1]], Vec2d(0.0, 0.0));
//     // assert_eq!(block.coords[[2, 1]], Vec2d(-0.5, 0.5));

//     // assert_eq!(block.coords[[0, 2]], Vec2d(1.0, 0.0));
//     // assert_eq!(block.coords[[1, 2]], Vec2d(0.5, 0.5));
//     // assert_eq!(block.coords[[2, 2]], Vec2d(0.0, 1.0));
// }
