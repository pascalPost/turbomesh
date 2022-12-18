// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

mod output;
mod test;
mod tfi;
mod types;

use types::Index;

use crate::output::write_block_to_legacy_vtk;
use crate::tfi::tfi_linear_2d;
use crate::types::Vec2d;
use crate::types::{Block2d, Scalar};

fn main() {
    let block_0 = Block2d::new(String::from("block_0"), (10, 20));

    // .corners(
    //     Vec2d(-1.0, -1.0),
    //     Vec2d(1.0, -1.0),
    //     Vec2d(1.0, 1.0),
    //     Vec2d(-1.0, 1.0),
    // );

    // .clustering(
    //     (equidistant(), constant()),
    //     (constant(), equidistant()),
    //     (equidistant(), constant()),
    //     (constant(), equidistant()),
    // );

    // let mut block_0 = Block2d::new(String::from("block_0"), (10, 20));

    // let x_min = -1.0;
    // let x_max = 1.0;
    // let y_min = -1.0;
    // let y_max = 1.0;

    // // i min edge
    // for i in 0..block_0.shape().0 {
    //     block_0.coords[[i, 0]] = Vec2d(
    //         x_min + (x_max - x_min) * i as Scalar / (block_0.shape().0 - 1) as Scalar,
    //         y_min,
    //     );
    // }

    // // i max edge
    // for i in 0..block_0.shape().0 {
    //     let j_len = block_0.shape().1;
    //     block_0.coords[[i, j_len - 1]] = Vec2d(
    //         x_min + (x_max - x_min) * i as Scalar / (block_0.shape().0 - 1) as Scalar,
    //         y_max,
    //     );
    // }

    // // j min edge
    // for j in 0..block_0.shape().1 {
    //     block_0.coords[[0, j]] = Vec2d(
    //         x_min,
    //         y_min + (y_max - y_min) * j as Scalar / (block_0.shape().1 - 1) as Scalar,
    //     );
    // }

    // // j max edge
    // for j in 0..block_0.shape().1 {
    //     let i_len = block_0.shape().0;
    //     block_0.coords[[i_len - 1, j]] = Vec2d(
    //         x_max,
    //         y_min + (y_max - y_min) * j as Scalar / (block_0.shape().1 - 1) as Scalar,
    //     );
    // }

    // tfi_linear_2d(&mut block_0);

    let file_name = "turbomesh.vtk";

    write_block_to_legacy_vtk(file_name, &block_0)
        .unwrap_or_else(|err| println!("Error encountered writing to file: {:?}", err));
}
