// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::cgns::interface::*;
use crate::types::Block2d;

/// mesh data structure
pub struct Mesh {
    pub blocks: Vec<Block2d>,
}

impl Mesh {
    /// creates a new mesh instance
    pub fn new() -> Self {
        Self { blocks: vec![] }
    }

    pub fn add_block(&mut self, block: Block2d) {
        self.blocks.push(block);
    }

    /// writes the mesh to a structured cgns
    pub fn save(&self, file_name: &str) -> Result<(), Box<dyn Error>> {
        let mut i_file = open(file_name, CgMode::Write)?;
        let i_base = base_write(i_file, "Base", 2, 2)?;

        for block in &self.blocks {
            let size: Vec<CgSizeT> = vec![
                block.points()[0] as CgSizeT,
                block.points()[1] as CgSizeT,
                block.points()[0] as CgSizeT - 1,
                block.points()[1] as CgSizeT - 1,
                0,
                0,
            ];

            let i_zone = zone_write(
                i_file,
                i_base,
                block.name.as_str(),
                &size,
                &ZoneType::Structured,
            )?;

            {
                let mut coords_x = vec![0.0; block.coords.size()];
                let mut coords_y = vec![0.0; block.coords.size()];

                for i in 0..block.coords.size() {
                    coords_x[i] = block.coords.as_slice()[i].0;
                    coords_y[i] = block.coords.as_slice()[i].1;
                }

                coord_write(i_file, i_base, i_zone, "CoordinateX", &coords_x.as_slice())?;
                coord_write(i_file, i_base, i_zone, "CoordinateY", &coords_y.as_slice())?;
            }
        }

        close(&mut i_file)?;
        Ok(())
    }
}
