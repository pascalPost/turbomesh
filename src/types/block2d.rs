use crate::types::{Array2d, Index, Vec2d};

#[derive(Debug)]
pub struct Block2d {
    pub name: String,
    pub coords: Array2d<Vec2d>,
}

// pub enum Edge {
//     IMin,
//     IMax,
//     JMin,
//     JMax,
// }

impl Block2d {
    pub fn new(name: String, shape: (Index, Index)) -> Self {
        Block2d {
            name,
            coords: Array2d::new((shape.0, shape.1)),
        }
    }

    pub fn shape(&self) -> (Index, Index) {
        self.coords.shape
    }

    pub fn _name(&self) -> &str {
        &self.name
    }

    // pub fn edge_mut(&mut self, edge: Edge) -> ndarray::ArrayViewMut1<'_, Vec2d> {
    //     match edge {
    //         Edge::IMin => self.coords.column_mut(0),
    //         Edge::IMax => self.coords.column_mut(self.coords.ncols() - 1),
    //         Edge::JMin => self.coords.row_mut(0),
    //         Edge::JMax => self.coords.row_mut(self.coords.nrows() - 1),
    //     }
    // }

    // fn mod_edge<F>(mut self, edge: Edge, f: F) -> Self
    // where
    //     F: FnMut(Vec2d) -> Vec2d,
    // {
    //     self.edge_mut(edge).mapv_inplace(f);
    //     self
    // }

    // fn mod_edge_test<F>(mut self, edge: Edge, f: F) -> Self
    // {
    //     self.edge_mut(edge).indexed_iter_mut();
    //     self
    // }
}
