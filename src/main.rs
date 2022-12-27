// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use turbomesh::{Block2d, Mesh, Vec2d};

fn main() {
    let mut mesh = Mesh::new();

    mesh.blocks = vec![
        Block2d::new(String::from("block_0"), (5, 5)).corners(
            Vec2d(0.0, 0.0),
            Vec2d(1.0, 0.0),
            Vec2d(1.0, 1.0),
            Vec2d(0.0, 1.0),
        ),
        Block2d::new(String::from("block_1"), (5, 5)).corners(
            Vec2d(1.0, 0.0),
            Vec2d(2.0, 0.0),
            Vec2d(2.0, 1.0),
            Vec2d(0.0, 1.0),
        ),
    ];

    // let file_name = "turbomesh.vtk";

    mesh.save("turbomesh.cgns").expect("cgns writing error.");

    // write_block_to_legacy_vtk(file_name, &block_0)
    //     .unwrap_or_else(|err| println!("Error encountered writing to file: {:?}", err));
}
