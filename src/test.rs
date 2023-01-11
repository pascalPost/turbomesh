// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::clustering::uniform_clustering;
use crate::geometry::Rectangle;
use crate::types::{Block2d, Edge, Vec2d};

// TODO add block creation given a geometry and the edge clusterings.
// The geometry must implement a Geometry trait specifying the mapping functions
// to the edges

/// test the block creation with a simple rectangle
#[test]
fn test_default_block_creation() {
    let rect = Rectangle::new(Vec2d(0.0, 0.0), Vec2d(1.0, 1.0));

    let block = Block2d::new(
        String::from("block"),
        (2, 2),
        Edge::new(uniform_clustering, rect.i_min()),
        Edge::new(uniform_clustering, rect.i_max()),
        Edge::new(uniform_clustering, rect.j_min()),
        Edge::new(uniform_clustering, rect.j_max()),
    );

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
