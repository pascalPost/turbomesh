use crate::types::{Block2d, Vec2d};

/// test the default block creation (done in computational space equidistant
/// with values 0<= x,y <=1)
#[test]
fn test_default_block_creation() {
    let block = Block2d::new(String::from("block"), (2, 2));

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

/// test the block creation with corners (and linear, equidistant edges)
#[test]
fn test_block_creation_with_corners() {
    let block = Block2d::new(String::from("block"), (2, 2)).corners(
        Vec2d(-1.0, -1.0),
        Vec2d(1.0, -1.0),
        Vec2d(1.0, 1.0),
        Vec2d(-1.0, 1.0),
    );

    // test block coordinates
    assert_eq!(block.coords[[0, 0]], Vec2d(-1.0, -1.0));
    assert_eq!(block.coords[[1, 0]], Vec2d(0.0, -1.0));
    assert_eq!(block.coords[[2, 0]], Vec2d(1.0, -1.0));

    assert_eq!(block.coords[[0, 1]], Vec2d(-1.0, 0.0));
    assert_eq!(block.coords[[1, 1]], Vec2d(0.0, 0.0));
    assert_eq!(block.coords[[2, 1]], Vec2d(1.0, 0.0));

    assert_eq!(block.coords[[0, 2]], Vec2d(-1.0, 1.0));
    assert_eq!(block.coords[[1, 2]], Vec2d(0.0, 1.0));
    assert_eq!(block.coords[[2, 2]], Vec2d(1.0, 1.0));
}

/// test linear block creation
///           x 2
///         /   \
///       x 1    x 3
///         \   /
///           x 0
#[test]
fn test_block_creation_linear() {
    let block = Block2d::new(String::from("block"), (2, 2)).corners(
        Vec2d(0.0, -1.0),
        Vec2d(-1.0, 0.0),
        Vec2d(0.0, 1.0),
        Vec2d(1.0, 0.0),
    );

    // test block coordinates
    assert_eq!(block.coords[[0, 0]], Vec2d(0.0, -1.0));
    assert_eq!(block.coords[[1, 0]], Vec2d(-0.5, -0.5));
    assert_eq!(block.coords[[2, 0]], Vec2d(-1.0, 0.0));

    assert_eq!(block.coords[[0, 1]], Vec2d(0.5, -0.5));
    assert_eq!(block.coords[[1, 1]], Vec2d(0.0, 0.0));
    assert_eq!(block.coords[[2, 1]], Vec2d(-0.5, 0.5));

    assert_eq!(block.coords[[0, 2]], Vec2d(1.0, 0.0));
    assert_eq!(block.coords[[1, 2]], Vec2d(0.5, 0.5));
    assert_eq!(block.coords[[2, 2]], Vec2d(0.0, 1.0));
}
