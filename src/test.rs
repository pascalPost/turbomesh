// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::block2d::{Edge, EdgeSegmentProps};
use crate::types::{Block2d, Scalar, Vec2d};

fn equidistant_clustering(u: &mut [Scalar]) {
    let n = u.len();
    for (i, u) in u.iter_mut().enumerate() {
        *u = i as Scalar / (n - 1) as Scalar
    }
    // u.iter_mut()
    //     .enumerate()
    //     .for_each(|(i, u)| );
}

// fn linear_edge(start: Vec2d, end: Vec2d) -> Box<dyn Fn(&[Scalar], &[Vec2d])> {
//     Box::new(|u, v| for (u, v) in u.iter().zip(v.iter_mut()) {

//     })
// }

/// test the default block creation (done in computational space equidistant
/// with values 0<= x,y <=1)
#[test]
fn test_default_block_creation() {
    let block = Block2d::new(
        String::from("block"),
        (2, 2),
        vec![
            EdgeSegmentProps {
                edge: Edge::IMin,
                start: 0,
                end: 3,
                clustering: Box::new(equidistant_clustering),
                mapping: Box::new(|u, x| {
                    for (u, x) in u.iter().zip(x.iter_mut()) {
                        *x = Vec2d(*u, 0.0);
                    }
                }),
            },
            EdgeSegmentProps {
                edge: Edge::IMax,
                start: 0,
                end: 3,
                clustering: Box::new(equidistant_clustering),
                mapping: Box::new(|u, x| {
                    for (u, x) in u.iter().zip(x.iter_mut()) {
                        *x = Vec2d(*u, 1.0);
                    }
                }),
            },
            EdgeSegmentProps {
                edge: Edge::JMin,
                start: 0,
                end: 3,
                clustering: Box::new(equidistant_clustering),
                mapping: Box::new(|u, x| {
                    for (u, x) in u.iter().zip(x.iter_mut()) {
                        *x = Vec2d(0.0, *u);
                    }
                }),
            },
            EdgeSegmentProps {
                edge: Edge::JMax,
                start: 0,
                end: 3,
                clustering: Box::new(equidistant_clustering),
                mapping: Box::new(|u, x| {
                    for (u, x) in u.iter().zip(x.iter_mut()) {
                        *x = Vec2d(1.0, *u);
                    }
                }),
            },
        ],
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
