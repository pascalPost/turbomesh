// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

mod tfi;
mod types;

use crate::tfi::tfi_linear_2d;
use crate::types::{Block2d, Edge, Scalar};
use std::fs::File;
use std::io::Write;

fn write_block_to_legacy_vtk(file_name: &str, block: &Block2d) -> Result<(), std::io::Error> {
    println!("Writing outlut to {file_name}");

    let i_len = block.i_len();
    let j_len = block.j_len();

    let mut output = File::create(file_name)?;
    writeln!(output, "# vtk DataFile Version 3.0")?;
    writeln!(output, "turbomesh")?;
    writeln!(output, "ASCII")?;
    writeln!(output, "DATASET STRUCTURED_GRID")?;
    writeln!(output, "DIMENSIONS {i_len} {j_len} 1")?;
    writeln!(output, "POINTS {} double", i_len * j_len)?;

    for j in 0..j_len {
        for i in 0..i_len {
            writeln!(
                output,
                "{} {} 0.0",
                block.coords[[i, j]].0,
                block.coords[[i, j]].1
            )?;
        }
    }

    Ok(())
}

fn main() {
    let i_len = 10;
    let j_len = 20;

    // let block_0 = Block2d::new(String::from("block_0"), (I, J))
    //     .mod_edge(Edge::IMin, |v| {})
    //     .mod_edge(Edge::IMax)
    //     .mod_edge(Edge::JMin)
    //     .mod_edge(Edge::JMax);

    let mut block_0 = Block2d::new(String::from("block_0"), (i_len, j_len));

    let x_min = -1.0;
    let x_max = 1.0;
    let y_min = -1.0;
    let y_max = 1.0;

    for (i, v) in block_0.edge_mut(Edge::IMin).indexed_iter_mut() {
        v.0 = x_min + (x_max - x_min) * i as Scalar / (i_len - 1) as Scalar;
        v.1 = y_min;
    }

    for (i, v) in block_0.edge_mut(Edge::IMax).indexed_iter_mut() {
        v.0 = x_min + (x_max - x_min) * i as Scalar / (i_len - 1) as Scalar;
        v.1 = y_max;
    }

    for (j, v) in block_0.edge_mut(Edge::JMin).indexed_iter_mut() {
        v.0 = x_min;
        v.1 = y_min + (y_max - y_min) * j as Scalar / (j_len - 1) as Scalar;
    }

    for (j, v) in block_0.edge_mut(Edge::JMax).indexed_iter_mut() {
        v.0 = x_max;
        v.1 = y_min + (y_max - y_min) * j as Scalar / (j_len - 1) as Scalar;
    }

    tfi_linear_2d(&mut block_0);

    let file_name = "turbomesh.vtk";

    write_block_to_legacy_vtk(file_name, &block_0)
        .unwrap_or_else(|err| println!("Error encountered writing to file: {:?}", err));
}
