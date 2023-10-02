// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::interpolation::FittingSpline;
use crate::types::{DiscretizableCurve, Scalar, Vec2d};
// use std::rc::Rc;

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

    // assumption of u : u in [0, 1], i.e. in computational space
    fn discretize(&self, u: &[Scalar], x: &mut [Vec2d]) {
        let dx = self.x_max - self.x_min;
        for (&u, x) in u.iter().zip(x.iter_mut()) {
            *x = self.x_min + u * dx;
        }
    }
}

fn add_to_given_clustering(u: &mut [Scalar], u_new: &Vec<Scalar>) {
    // add to given clustering to possibly existing clustering
    assert!(u_new.starts_with(&[0.0]));
    let mut start = 0.0;
    if u[0].is_nan() {
        u[0] = 0.0;
    } else {
        start = u[0];
    }
    u[1..]
        .iter_mut()
        .zip(u_new[1..].iter())
        .for_each(|(u, u_ref)| *u = start + *u_ref);
}

impl DiscretizableCurve for Line2d {
    fn arclength(&self) -> Scalar {
        (self.x_max - self.x_min).abs()
    }

    fn discretize_curve(
        &self,
        clustering: &impl crate::types::ClusteringFunction,
        u: &mut [Scalar],
        x: &mut [Vec2d],
    ) {
        // clustering in computational space
        let mut u_new = clustering.get_clustering(u.len());

        // compute discretized coordinates on line
        self.discretize(&u_new[..], x);

        // convert clustering to physical space
        let scale = self.arclength();
        u_new.iter_mut().for_each(|u| *u *= scale);

        add_to_given_clustering(u, &u_new);
    }
}

#[derive(Clone)]
pub struct Spline {
    spline: FittingSpline<2>,
}

impl Spline {
    pub fn new(spline: FittingSpline<2>) -> Self {
        Self { spline }
    }

    pub fn interpolate_val(&self, u: Scalar) -> Vec2d {
        let x = self.spline.interpolate(&[u])[0];
        Vec2d(x[0], x[1])
    }

    pub fn interpolate(&self, u: &[f64]) -> Box<[[f64; 2]]> {
        self.spline.interpolate(u)
    }

    fn interpolate_into(&self, u: &[Scalar], x: &mut [Vec2d]) {
        unsafe {
            let (_, x_conv, _) = x.align_to_mut::<[Scalar; 2]>();
            self.spline.interpolate_into(u, x_conv);
        }
    }
}

impl DiscretizableCurve for Spline {
    fn discretize_curve(
        &self,
        clustering: &impl crate::types::ClusteringFunction,
        u: &mut [Scalar],
        x: &mut [Vec2d],
    ) {
        // clustering in computational space
        let mut u_new = clustering.get_clustering(u.len());

        // (x,y)=f(u) where from 0<=u<=1
        self.interpolate_into(&u_new[..], x);

        // convert clustering to physical space
        let scale = self.arclength();
        u_new.iter_mut().for_each(|u| *u /= scale);

        add_to_given_clustering(u, &u_new);
    }

    fn arclength(&self) -> Scalar {
        self.spline.integrate(0.0, 1.0)
    }
}

// pub struct SplineRef {
//     spline: Rc<FittingSpline<2>>,
//     v_start: Scalar,
//     v_end: Scalar,
// }

// impl SplineRef {
//     pub fn new(spline: &Rc<FittingSpline<2>>, v_start: Scalar, v_end: Scalar) -> Self {
//         Self {
//             spline: spline.clone(),
//             v_start,
//             v_end,
//         }
//     }

//     pub fn split_at(&self, v: Scalar) -> (Self, Self) {
//         assert!(
//             self.v_start < v && self.v_end > v,
//             "Spline must be split in between its range [{}, {}]",
//             self.v_start,
//             self.v_end
//         );
//         (
//             SplineRef::new(&self.spline, self.v_start, v),
//             SplineRef::new(&self.spline, v, self.v_end),
//         )
//     }
// }

// impl MappingFunction for SplineRef {
//     /// (x,y)=f(u) where from 0<=u<=1
//     fn computational_to_physical(&self, u: &[Scalar], x: &mut [Vec2d]) {
//         // map 0<=u<=1 to the arclength of the spline v
//         let mut v = vec![Scalar::NAN; u.len()];

//         let v_scale = self.v_end - self.v_start;
//         for (&u, v) in u.iter().zip(v.iter_mut()) {
//             *v = self.v_start + u * v_scale;
//         }

//         unsafe {
//             let (_, x_conv, _) = x.align_to_mut::<[Scalar; 2]>();
//             self.spline.interpolate_into(v.as_slice(), x_conv);
//         }
//     }
// }

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
