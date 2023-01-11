// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

mod array2d;
mod block2d;
mod vec2d;

pub type Scalar = f64;
pub type Index = usize;

pub use super::types::array2d::Array2d;
pub use crate::types::block2d::Block2d;
pub use crate::types::block2d::Edge;
pub use crate::types::block2d::Mapping;
pub use crate::types::vec2d::Vec2d;
