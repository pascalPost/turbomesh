// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

mod block2d;
mod boundary;
mod edge;
mod segment;
mod vec2d;

pub type Scalar = f64;
pub type Index = usize;

pub use crate::types::block2d::{Block2d, EdgeIndex};
pub use crate::types::boundary::{
    edge_view_mut, BlockBoundary, BlockBoundaryRange, BlockBoundaryRangeNew, BlockConnection,
};
pub use crate::types::edge::{BlockEdgeData, Edge, EdgeView};
pub use crate::types::segment::ClusteringFunction;
pub use crate::types::segment::DiscretizableCurve;
pub use crate::types::segment::Segment;
pub use crate::types::segment::SegmentFunction;
pub use crate::types::vec2d::Vec2d;
