// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{ClusteringFunction, Scalar};
use serde::Deserialize;

pub struct UniformClustering {}

impl UniformClustering {
    pub fn new() -> Self {
        Self {}
    }
}

impl ClusteringFunction for UniformClustering {
    fn get_clustering(&self, points: usize) -> Vec<Scalar> {
        let n = points;
        let mut u = vec![0.0; n];
        u.iter_mut()
            .enumerate()
            .for_each(|(i, u)| *u = i as Scalar / (n - 1) as Scalar);
        u
    }
}

/// SingleHyperbolicTangentClustering implements the hyperbolic tangent clustering function
/// matching the specified spacing of the first cell approximately.
///
/// Details can be find in:
///
/// Vinokur, Marcel. “On One-Dimensional Stretching Functions for Finite-Difference Calculations.”
/// Journal of Computational Physics 50, no. 2 (May 1983): 215–34.
/// https://doi.org/10.1016/0021-9991(83)90065-7.
///
/// Thompson, Joe F. “A General Three-Dimensional Elliptic Grid Generation System on a Composite Block Structure.”
/// Computer Methods in Applied Mechanics and Engineering 64, no. 1–3 (October 1987): 377–411.
/// https://doi.org/10.1016/0045-7825(87)90047-8.
#[derive(Clone, Copy, Debug, Deserialize)]
pub struct SingleHyperbolicTangentClustering {
    delta_s: Scalar,
}

impl SingleHyperbolicTangentClustering {
    pub fn new(delta_s: Scalar) -> Self {
        Self { delta_s }
    }
}

impl ClusteringFunction for SingleHyperbolicTangentClustering {
    fn get_clustering(&self, points: usize) -> Vec<Scalar> {
        // TODO enhance implementation: it seems that only the derivative at the wall is specified not
        // the spacing itself, see https://www.cfd-online.com/Wiki/Structured_mesh_generation
        // or in https://www.osti.gov/etdeweb/servlets/purl/632740 page 9:
        // ds/deta is approximately equal to the normalized height of the first cell delta_1 / delta_s
        // For this, the mor general implementation of the Vinokur clustering might be needed, where
        // the spacing can be specified at any point

        // compute B
        let I = points - 1;
        let B = I as Scalar * self.delta_s;

        let y = 1.0 / B;

        // eq. 63 to 67 in Vinokur 1983
        let delta = if y < 2.7829681 {
            let y_bar = y - 1.0;
            (6.0 * y_bar).sqrt()
                * (1.0
                    + y_bar
                        * (-0.15
                            + y_bar
                                * (0.057321429
                                    + y_bar
                                        * (-0.024907295
                                            + y_bar * (0.0077424461 - 0.0010794123 * y_bar)))))
        } else {
            let w = 1.0 / y - 0.028527431;
            let v = y.ln();
            v + (1.0 + 1.0 / v) * (2.0 * v).ln() - 0.02041793
                + w * (0.24902722 + w * (1.9496443 + w * (-2.6294547 + 8.56795911 * w)))
        };

        let mut xi: Vec<Scalar> = (0..points)
            .map(|x| x as Scalar / (points - 1) as Scalar)
            .collect();

        for xi in xi[1..].iter_mut() {
            let s = 1.0 + (0.5 * delta * (*xi - 1.0)).tanh() / (0.5 * delta).tanh();
            *xi = s;
        }

        assert_eq!(xi[0], 0.0);
        assert_eq!(xi[points - 1], 1.0);

        xi
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
pub struct SingleExponentialClustering {
    alpha: Scalar,
}

impl SingleExponentialClustering {
    pub fn new(alpha: Scalar) -> Self {
        Self { alpha }
    }
}

impl ClusteringFunction for SingleExponentialClustering {
    fn get_clustering(&self, points: usize) -> Vec<Scalar> {
        let alpha = self.alpha;

        let n = points;

        let mut x: Vec<Scalar> = (0..n).map(|x| x as Scalar / (n - 1) as Scalar).collect();
        for v in x.iter_mut() {
            *v = ((alpha * *v).exp() - 1.0) / (alpha.exp() - 1.0);
        }

        x
    }
}

/// Roberts cluster function, see
/// https://github.com/luohancfd/CFCFD-NG/blob/dev/lib/nm/source/fobject.cxx
/// alpha = 0.5 cluster at both ends
/// alpha = 0.0 cluster toward t=1.0
/// stretching factor 1.0 < beta < +inf, closer to 1.0 gives stronger clustering.
/// for free function implementation, see roberts_clustering
#[derive(Clone, Copy, Debug, Deserialize)]
pub struct RobertsClustering {
    alpha: Scalar,
    beta: Scalar,
}

/// Roberts cluster function, see
/// https://github.com/luohancfd/CFCFD-NG/blob/dev/lib/nm/source/fobject.cxx
/// alpha = 0.5 cluster at both ends
/// alpha = 0.0 cluster toward t=1.0
/// stretching factor 1.0 < beta < +inf, closer to 1.0 gives stronger clustering
impl RobertsClustering {
    pub fn new(alpha: Scalar, beta: Scalar) -> Self {
        Self { alpha, beta }
    }
}

impl ClusteringFunction for RobertsClustering {
    fn get_clustering(&self, points: usize) -> Vec<Scalar> {
        let alpha = self.alpha;
        let beta = self.beta;

        let n = points;

        // init w/ values between 0..1 by u = i / n
        let mut u: Vec<Scalar> = (0..n).map(|x| x as Scalar / (n - 1) as Scalar).collect();

        u.iter_mut().for_each(|u| {
            let tmp = ((beta + 1.0) / (beta - 1.0)).powf((*u - alpha) / (1.0 - alpha));
            let tbar = (beta + 2.0 * alpha) * tmp - beta + 2.0 * alpha;
            *u = tbar / ((2.0 * alpha + 1.0) * (1.0 + tmp));
        });

        u
    }
}
