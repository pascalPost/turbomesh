// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::types::Index;
use std::alloc;
use std::mem;
use std::ptr::NonNull;

/// contigous 2d array with column major memory layout
pub struct Array2d<T> {
    pub shape: [Index; 2],
    ptr: NonNull<T>,
}

impl<T> Array2d<T> {
    pub fn size(&self) -> Index {
        self.shape[0] * self.shape[1]
    }

    pub fn new(shape: [Index; 2]) -> Self {
        assert!(mem::size_of::<T>() != 0, "We're not ready to handle ZSTs");
        assert!(shape[0] * shape[1] > 0, "Size must be greater than zero.");

        Self {
            shape,
            ptr: {
                let layout = alloc::Layout::array::<T>(shape[0] * shape[1]).unwrap();

                unsafe {
                    let ptr = alloc::alloc_zeroed(layout);
                    if ptr.is_null() {
                        alloc::handle_alloc_error(layout);
                    }

                    // If allocation fails, `new_ptr` will be null, in which case we abort.
                    match NonNull::new(ptr as *mut T) {
                        Some(p) => p,
                        None => alloc::handle_alloc_error(layout),
                    }
                }
            },
        }
    }

    fn index(&self, idx: [usize; 2]) -> Index {
        let i = idx[0] + idx[1] * self.shape[0];
        assert!(i < self.size());
        i
    }
}

impl<T> Drop for Array2d<T> {
    fn drop(&mut self) {
        let layout = alloc::Layout::array::<T>(self.size()).unwrap();
        unsafe {
            alloc::dealloc(self.ptr.as_ptr() as *mut u8, layout);
        }
    }
}

impl<T> std::ops::Index<[usize; 2]> for Array2d<T> {
    type Output = T;

    fn index(&self, idx: [usize; 2]) -> &T {
        unsafe { &*self.ptr.as_ptr().add(self.index(idx)) }
    }
}

impl<T> std::ops::IndexMut<[usize; 2]> for Array2d<T> {
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut T {
        unsafe { &mut *self.ptr.as_ptr().add(self.index(idx)) }
    }
}

impl<T> std::fmt::Debug for Array2d<T>
where
    T: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for j in 0..self.shape[1] {
            for i in 0..self.shape[0] {
                writeln!(f, "({} {}) [{:?}]", i, j, self[[i, j]])?;
            }
        }
        Ok(())
    }
}
