// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use crate::cgns::interface::*;
use crate::smoothing::{smooth_block, smooth_mesh, SmoothingMethod};
use crate::types::{Block2d, BlockBoundary};

/// mesh data structure
pub struct Mesh {
    pub blocks: Vec<Block2d>,
    pub edges: Vec<BlockBoundary>,
}

impl Mesh {
    /// creates a new mesh instance
    pub fn new() -> Self {
        Self {
            blocks: vec![],
            edges: vec![],
        }
    }

    /// add the given block to the mesh
    pub fn add_block(&mut self, block: Block2d) {
        self.blocks.push(block);
    }

    /// returns the index of the block w/ given name if found
    pub fn block_id(&self, name: &str) -> Option<usize> {
        self.blocks.iter().position(|block| block.name == name)
    }

    /// smooth the given mesh with the specified method
    pub fn smooth(&mut self, method: &SmoothingMethod) {
        match method {
            SmoothingMethod::None => {}
            SmoothingMethod::Global { iterations } => smooth_mesh(self, *iterations).unwrap(),
            SmoothingMethod::BlockInternal { iterations } => {
                self.blocks.iter_mut().for_each(|block| {
                    println!("block: {}", block.name);
                    smooth_block(block, *iterations).unwrap()
                });
            }
        }
    }

    /// returns the total number of points of the mesh
    pub fn points(&self) -> usize {
        self.blocks.iter().map(|b| b.coords.len()).sum()
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
                // assumption coordinates in row-major (C) ordering -> needs
                // transformation for cgns
                // TODO add detection

                let coords_x: Vec<f64> = block.coords.t().iter().map(|x| x.0).collect();
                let coords_y: Vec<f64> = block.coords.t().iter().map(|x| x.1).collect();

                coord_write(i_file, i_base, i_zone, "CoordinateX", &coords_x.as_slice())?;
                coord_write(i_file, i_base, i_zone, "CoordinateY", &coords_y.as_slice())?;
            }
        }

        close(&mut i_file)?;
        Ok(())
    }
}
