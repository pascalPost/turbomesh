// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

pub type Scalar = f64;
pub type Index = usize;

#[derive(Debug, Clone, Copy)]
pub struct Vec2d(pub Scalar, pub Scalar);

impl Vec2d {
    pub fn abs(&self) -> Scalar {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }
}

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

#[derive(Debug)]
pub struct Block2d {
    name: String,
    pub coords: ndarray::Array2<Vec2d>,
}

pub enum Edge {
    IMin,
    IMax,
    JMin,
    JMax,
}

impl Block2d {
    pub fn new(name: String, size: (Index, Index)) -> Self {
        Block2d {
            name,
            coords: ndarray::Array2::from_shape_simple_fn((size.0, size.1), || {
                Vec2d(Scalar::NAN, Scalar::NAN)
            }),
        }
    }

    fn size(&self) -> (Index, Index) {
        self.coords.dim()
    }

    pub fn i_len(&self) -> Index {
        self.coords.len_of(ndarray::Axis(0))
    }

    pub fn j_len(&self) -> Index {
        self.coords.len_of(ndarray::Axis(1))
    }

    pub fn edge_mut(&mut self, edge: Edge) -> ndarray::ArrayViewMut1<'_, Vec2d> {
        match edge {
            Edge::IMin => self.coords.column_mut(0),
            Edge::IMax => self.coords.column_mut(self.coords.ncols() - 1),
            Edge::JMin => self.coords.row_mut(0),
            Edge::JMax => self.coords.row_mut(self.coords.nrows() - 1),
        }
    }

    fn mod_edge<F>(mut self, edge: Edge, f: F) -> Self
    where
        F: FnMut(Vec2d) -> Vec2d,
    {
        self.edge_mut(edge).mapv_inplace(f);
        self
    }

    // fn mod_edge_test<F>(mut self, edge: Edge, f: F) -> Self
    // {
    //     self.edge_mut(edge).indexed_iter_mut();
    //     self
    // }
}
