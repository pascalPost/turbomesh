// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

// use crate::types::Block2d;
// use crate::Mesh;
// use std::fs::File;
// use std::io::Write;

// pub fn write_block_to_legacy_vtk(file_name: &str, block: &Block2d) -> Result<(), std::io::Error> {
//     println!("Writing outlut to {file_name}");

//     let [i_len, j_len] = block.points();

//     let mut output = File::create(file_name)?;
//     writeln!(output, "# vtk DataFile Version 3.0")?;
//     writeln!(output, "turbomesh")?;
//     writeln!(output, "ASCII")?;
//     writeln!(output, "DATASET STRUCTURED_GRID")?;
//     writeln!(output, "DIMENSIONS {i_len} {j_len} 1")?;
//     writeln!(output, "POINTS {} double", i_len * j_len)?;

//     for j in 0..j_len {
//         for i in 0..i_len {
//             writeln!(
//                 output,
//                 "{} {} 0.0",
//                 block.coords[[i, j]].0,
//                 block.coords[[i, j]].1
//             )?;
//         }
//     }

//     Ok(())
// }

// pub fn write_mesh_to_legacy_vtk(file_name: &str, mesh: &Mesh) -> Result<(), std::io::Error> {
//     println!("Writing outlut to {file_name}");

//     let [i_len, j_len] = block.points();

//     let mut output = File::create(file_name)?;
//     writeln!(output, "# vtk DataFile Version 3.0")?;
//     writeln!(output, "turbomesh")?;
//     writeln!(output, "ASCII")?;
//     writeln!(output, "DATASET STRUCTURED_GRID")?;
//     writeln!(output, "DIMENSIONS {i_len} {j_len} 1")?;
//     writeln!(output, "POINTS {} double", i_len * j_len)?;

//     for j in 0..j_len {
//         for i in 0..i_len {
//             writeln!(
//                 output,
//                 "{} {} 0.0",
//                 block.coords[[i, j]].0,
//                 block.coords[[i, j]].1
//             )?;
//         }
//     }

//     Ok(())
// }
