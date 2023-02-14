// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::Scalar;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2d(pub Scalar, pub Scalar);

// impl [Vec2d] {}

// impl std::slice {
//     pub fn as_slice(&self) -> &[[Scalar; 2]] {
//         unsafe {
//             std::slice::from_raw_parts(
//                 self as *const Self as *const [Scalar; 2],
//                 std::mem::size_of::<Scalar>() * 2,
//             )
//         }
//     }
// }

impl Vec2d {
    pub fn abs(&self) -> Scalar {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }

    // /// reinterpret as a slice of Scalar
    // pub fn as_slice(&self) -> &[[Scalar; 2]] {
    //     // assert!(std::mem::size_of::Vec2d() <= isize::MAX as _);
    //     unsafe {
    //         std::slice::from_raw_parts(
    //             self as *const Self as *const [Scalar; 2],
    //             std::mem::size_of::<Scalar>() * 2,
    //         )
    //     }
    // }

    // /// reinterpret as a mut slice of Scalar
    // pub fn as_mut_slice(&mut self) -> &mut [[Scalar; 2]] {
    //     // assert!(std::mem::size_of::Vec2d() <= isize::MAX as _);
    //     unsafe {
    //         std::slice::from_raw_parts_mut(
    //             self as *mut Self as *mut [Scalar; 2],
    //             std::mem::size_of::<Scalar>() * 2,
    //         )
    //     }
    // }
}

impl std::convert::Into<[Scalar; 2]> for Vec2d {
    fn into(self) -> [Scalar; 2] {
        [self.0, self.1]
    }
}

// impl std::convert::Into<[[Scalar; 2]]> for [Vec2d] {
//     fn into(self) -> [Scalar; 2] {
//         [self.0, self.1]
//     }
// }

impl std::ops::Add for Vec2d {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl std::ops::Sub for Vec2d {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl std::ops::Mul<Vec2d> for Scalar {
    type Output = Vec2d;

    fn mul(self, rhs: Vec2d) -> Self::Output {
        Vec2d(self * rhs.0, self * rhs.1)
    }
}

pub trait From<T>: Sized {
    fn from(value: T) -> Self;
}

impl From<(Scalar, Scalar)> for Vec2d {
    fn from(value: (Scalar, Scalar)) -> Self {
        Vec2d(value.0, value.1)
    }
}

impl std::fmt::Display for Vec2d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

use float_cmp::{ApproxEq, F64Margin};
// use float_cmp::Margin;

impl ApproxEq for &Vec2d {
    // TODO implement Margin for f32
    type Margin = F64Margin;

    fn approx_eq<T: Into<Self::Margin>>(self, other: Self, margin: T) -> bool {
        let margin = margin.into();
        self.0.approx_eq(other.0, margin) && self.1.approx_eq(other.1, margin)
    }
}
