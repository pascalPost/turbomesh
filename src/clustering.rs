// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Index, Scalar};

pub fn uniform_clustering(u: &mut [Scalar]) {
    let n = u.len();
    for (i, u) in u.iter_mut().enumerate() {
        *u = i as Scalar / (n - 1) as Scalar
    }
}

pub fn single_exponential_clustering(n: Index, a: Scalar) -> Box<[Scalar]> {
    let mut x: Vec<Scalar> = (0..n).map(|x| x as Scalar / (n - 1) as Scalar).collect();
    for v in x.iter_mut() {
        *v = ((a * *v).exp() - 1.0) / (a.exp() - 1.0);
    }
    x.into_boxed_slice()
}

/// Roberts cluster function, see
/// https://github.com/luohancfd/CFCFD-NG/blob/dev/lib/nm/source/fobject.cxx
/// alpha = 0.5 cluster at both ends
/// alpha = 0.0 cluster toward t=1.0
/// stretching factor 1.0 < beta < +inf, closer to 1.0 gives stronger clustering
pub fn roberts_clustering(n: Index, alpha: Scalar, beta: Scalar) -> Box<[Scalar]> {
    let mut x: Vec<Scalar> = (0..n).map(|x| x as Scalar / (n - 1) as Scalar).collect();

    for v in x.iter_mut() {
        let tmp = ((beta + 1.0) / (beta - 1.0)).powf((*v - alpha) / (1.0 - alpha));
        let tbar = (beta + 2.0 * alpha) * tmp - beta + 2.0 * alpha;
        *v = tbar / ((2.0 * alpha + 1.0) * (1.0 + tmp))
    }

    x.into_boxed_slice()
}
