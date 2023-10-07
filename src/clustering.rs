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

// pub fn single_exponential_clustering(n: Index, a: Scalar) -> Box<[Scalar]> {
//     let mut x: Vec<Scalar> = (0..n).map(|x| x as Scalar / (n - 1) as Scalar).collect();
//     for v in x.iter_mut() {
//         *v = ((a * *v).exp() - 1.0) / (a.exp() - 1.0);
//     }
//     x.into_boxed_slice()
// }

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
