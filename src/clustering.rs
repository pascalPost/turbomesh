// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::{Index, Scalar};

pub fn uniform_clustering(n: Index) -> ndarray::Array1<Scalar> {
    ndarray::Array::linspace(0.0, 1.0, n)
}

pub fn single_exponential_clustering(n: Index, a: Scalar) -> ndarray::Array1<Scalar> {
    let x = ndarray::Array::linspace(0.0, 1.0, n);
    x.mapv_into(|v| ((a * v).exp() - 1.0) / (a.exp() - 1.0))
}

/// Roberts cluster function, see
/// https://github.com/luohancfd/CFCFD-NG/blob/dev/lib/nm/source/fobject.cxx
/// alpha = 0.5 cluster at both ends
/// alpha = 0.0 cluster toward t=1.0
/// stretching factor 1.0 < beta < +inf, closer to 1.0 gives stronger clustering
pub fn roberts_clustering(n: Index, alpha: Scalar, beta: Scalar) -> ndarray::Array1<Scalar> {
    let x = ndarray::Array::linspace(0.0, 1.0, n);
    x.mapv_into(|v| {
        let tmp = ((beta + 1.0) / (beta - 1.0)).powf((v - alpha) / (1.0 - alpha));
        let tbar = (beta + 2.0 * alpha) * tmp - beta + 2.0 * alpha;
        tbar / ((2.0 * alpha + 1.0) * (1.0 + tmp))
    })
}

pub fn uniform_distribution(start: Scalar, end: Scalar, n: Index) -> ndarray::Array1<Scalar> {
    ndarray::Array::linspace(start, end, n)
}

// fn uniform_distribution
