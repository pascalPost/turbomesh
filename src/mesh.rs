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

    /// writes the mesh to a structured cgns
    pub fn save(&self, file_name: &str) -> Result<(), Box<dyn Error>> {
        let mut i_file = open(file_name, CgMode::Write)?;
        let i_base = base_write(i_file, "Base", 2, 2)?;

        for block in &self.blocks {
            let size: Vec<CgSizeT> = vec![
                block.points()[0] as CgSizeT,
                block.points()[0] as CgSizeT - 1,
                0,
                block.points()[1] as CgSizeT,
                block.points()[1] as CgSizeT - 1,
                0,
            ];

            let _i_zone = zone_write(
                i_file,
                i_base,
                block.name.as_str(),
                &size,
                ZoneType::Structured,
            )?;
        }

        close(&mut i_file)?;
        Ok(())
    }
}
