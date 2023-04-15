// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

/// represents a point index in a block. INCLUDES_GHOST_LAYER indicates if the
/// ghost point layer is included or not.
pub struct PointIndex<
    IndexT: num_traits::int::PrimInt,
    const INCLUDES_GHOST_LAYER: bool,
    const DIMS: usize = 2,
> {
    pub index: [IndexT; DIMS],
}

impl<IndexT: num_traits::int::PrimInt, const INCLUDES_GHOST_LAYER: bool, const DIMS: usize>
    PointIndex<IndexT, INCLUDES_GHOST_LAYER, DIMS>
{
    fn new(index: [IndexT; DIMS]) -> Self {
        Self { index }
    }
}

impl<IndexT: num_traits::int::PrimInt + std::ops::SubAssign, const DIMS: usize>
    PointIndex<IndexT, true, DIMS>
{
    fn exclude_ghost_layer(&self) -> PointIndex<IndexT, false, DIMS> {
        let mut index = self.index.clone();
        index.iter_mut().for_each(|x| *x -= IndexT::one());
        PointIndex::<IndexT, false, DIMS>::new(index)
    }
}

impl<IndexT: num_traits::int::PrimInt + std::ops::AddAssign, const DIMS: usize>
    PointIndex<IndexT, false, DIMS>
{
    fn include_ghost_layer(&self) -> PointIndex<IndexT, true, DIMS> {
        let mut index = self.index.clone();
        index.iter_mut().for_each(|x| *x += IndexT::one());
        PointIndex::<IndexT, true, DIMS>::new(index)
    }
}
