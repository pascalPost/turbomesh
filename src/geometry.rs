// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Scalar, Vec2d};
use crate::Mapping;

/// representing a straigt line from x_min to x_max
pub struct Line2d {
    x_min: Vec2d,
    x_max: Vec2d,
}

impl Line2d {
    /// returns a Line2d from x0 to x1
    pub fn new(x0: Vec2d, x1: Vec2d) -> Self {
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

pub struct Rectangle {
    origin: Vec2d,
    size: Vec2d,
}

impl Rectangle {
    pub fn new(origin: Vec2d, size: Vec2d) -> Self {
        assert!(
            size.0 > 0.0 && size.1 > 0.0,
            "Size must be greater than zero"
        );
        Rectangle { origin, size }
    }

    pub fn i_min(&self) -> Line2d {
        Line2d::new(self.origin, self.origin + Vec2d(self.size.0, 0.0))
    }

    pub fn i_max(&self) -> Line2d {
        let x0 = self.origin + Vec2d(0.0, self.size.1);
        Line2d::new(x0, x0 + Vec2d(self.size.0, 0.0))
    }

    pub fn j_min(&self) -> Line2d {
        Line2d::new(self.origin, self.origin + Vec2d(0.0, self.size.1))
    }

    pub fn j_max(&self) -> Line2d {
        let x0 = self.origin + Vec2d(self.size.0, 0.0);
        Line2d::new(x0, x0 + Vec2d(0.0, self.size.1))
    }
}
