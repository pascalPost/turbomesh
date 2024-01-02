// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::clustering::{RobertsClustering, SingleHyperbolicTangentClustering, UniformClustering};
use crate::geometry::{Line2d, Spline};
use crate::interpolation::FittingSpline;
use crate::smoothing::SmoothingMethod;
use crate::types::{
    BlockBoundary, BlockBoundaryRange, BlockConnection, Edge, EdgeIndex, EdgeView,
    PeriodicBlockConnection,
};
use crate::{Block2d, Mesh, Segment, Vec2d};
use libc::initgroups;
use ndarray::Array;
use plotters::prelude::*;
use serde::Deserialize;
use std::error::Error;

/// read the profile coordinates from given file and returns them as a vec
pub fn read_profile<T>(file_name: &str) -> Result<Vec<[T; 2]>, Box<dyn Error>>
where
    T: num_traits::float::Float + std::str::FromStr,
    <T as std::str::FromStr>::Err: std::error::Error,
    <T as std::str::FromStr>::Err: 'static,
{
    let mut v: Vec<[T; 2]> = vec![];

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .delimiter(b' ')
        .from_path(file_name)?;

    // does not work as iter jumps to end of file
    // TODO insert reserve for better efficiency
    // TODO implement own CSV reader as this seems to be overly complicated
    // v.reserve(iter.count());

    for result in rdr.records() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here.
        let record = result?;

        assert!(record.len() == 2);

        let x = record[0].parse::<T>()?;
        let y = record[1].parse::<T>()?;

        v.push([x, y]);
    }

    assert!(v.len() > 0, "parsed profile data is empty.");

    Ok(v)
}

pub fn blade_profile(
    ps_csv_path: &str,
    ss_csv_path: &str,
) -> Result<(Spline, Spline), Box<dyn Error>> {
    // read suction and pressure side coordinates from csv and save into vec

    let mut ps = read_profile(ps_csv_path).expect("error in reading pressure side coordinates");

    if ps.first().unwrap()[0] > ps.last().unwrap()[0] {
        ps.reverse();
    }

    let mut ss = read_profile(ss_csv_path).expect("Error in reading suction side coordinates");

    if ss.first().unwrap()[0] > ss.last().unwrap()[0] {
        ss.reverse();
    }

    assert!(
        ps[0] == ss[0] || ps[0] == ss[ss.len() - 1],
        "Leading edge of suction and pressure side must be equal."
    );
    assert!(
        ps[ps.len() - 1] == ss[ss.len() - 1] || ps[ps.len() - 1] == ss[0],
        "Trailing edge of suction and pressure side must be equal."
    );

    // interpolate
    let ps_spline =
        FittingSpline::new(ps.as_slice(), None, None, 3).expect("Error in spline interpolation");
    let ss_spline =
        FittingSpline::new(ss.as_slice(), None, None, 3).expect("Error in spline interpolation");

    // TODO enhance interpolation by enforcing matching derivatives at the boundaries

    Ok((Spline::new(ps_spline), Spline::new(ss_spline)))
}

#[derive(Deserialize, Debug)]
#[serde(tag = "type")]
enum BladeClusteringFunction {
    RobertsClustering(RobertsClustering),
}

#[derive(Deserialize, Debug)]
pub struct TurbineTemplate {
    ps_csv_path: std::path::PathBuf,
    ss_csv_path: std::path::PathBuf,
    pitch: f64,
    blade_clustering_function: BladeClusteringFunction,
    num_cells_blade_half: usize,
    num_cells_middle_blocks_on_blade_half: usize,
    num_cells_next_to_middle_blocks_on_blade: usize,
    num_cells_away_from_blade: usize,
    num_cells_inlet: usize,
    num_cells_outlet: usize,
    smoothing: SmoothingMethod,
}

impl TurbineTemplate {
    /// appends the given root path to the file paths of the template
    pub fn append_root_path(&mut self, root_path: &std::path::Path) {
        self.ps_csv_path = root_path.join(&self.ps_csv_path);
        self.ss_csv_path = root_path.join(&self.ss_csv_path);
    }

    /// runs the template and returns the resulting geometry and mesh
    pub fn run(&self) {
        let ps_csv_path = self.ps_csv_path.to_str().unwrap();
        let ss_csv_path = self.ss_csv_path.to_str().unwrap();

        let pitch = self.pitch;
        let num_cells_blade_half = self.num_cells_blade_half;
        let num_cells_middle_blocks_on_blade_half = self.num_cells_middle_blocks_on_blade_half;
        let num_cells_next_to_middle_blocks_on_blade =
            self.num_cells_next_to_middle_blocks_on_blade;
        let num_cells_away_from_blade = self.num_cells_away_from_blade;

        // TODO can be computed automatically based on the average size
        let num_cells_inlet = self.num_cells_inlet;
        let num_cells_outlet = self.num_cells_outlet;

        let (ps_spline, ss_spline) = blade_profile(ps_csv_path, ss_csv_path)
            .expect("Error in profile input and interpolation");

        // blade discretization on which the 2d blocking is based
        let blade_clustering_function = match self.blade_clustering_function {
            BladeClusteringFunction::RobertsClustering(roberts_clustering) => roberts_clustering,
        };

        let ps_edge = EdgeView::new(Edge::new(
            "Pressure_Side_Edge".to_string(),
            &vec![Box::new(Segment::new(
                num_cells_blade_half + 1,
                blade_clustering_function,
                ps_spline.clone(),
            ))],
        ));

        let ss_edge = EdgeView::new(Edge::new(
            "Suction_Side_Edge".to_string(),
            &vec![Box::new(Segment::new(
                num_cells_blade_half + 1,
                blade_clustering_function,
                ss_spline.clone(),
            ))],
        ));

        // o-grid around blade for viscous computations

        // TODO make this runtime dependent

        // TODO replace with a percentage value of the chord length
        let d = 0.01;

        // compute o-grid target by projecting the blade normal outward

        let mut ss_outer = vec![Vec2d(0.0, 0.0); num_cells_blade_half + 1];

        // handle all points except the first and last
        for i in 1..num_cells_blade_half {
            let x_i = ss_edge.point_coord(i);
            let x_im1 = ss_edge.point_coord(i - 1);
            let x_ip1 = ss_edge.point_coord(i + 1);

            let x_xi = 0.5 * (x_ip1 - x_im1);

            let n = (1.0 / x_xi.abs()) * Vec2d(-x_xi.1, x_xi.0);

            ss_outer[i] = x_i + d * n;
        }

        // handle first point
        {
            let i = 0;
            let x_i = ss_edge.point_coord(i);
            let x_ip1 = ss_edge.point_coord(i + 1);

            let x_xi = -x_i + x_ip1;

            let n = (1.0 / x_xi.abs()) * Vec2d(-x_xi.1, x_xi.0);

            ss_outer[i] = x_i + d * n;
        }

        // handle last point
        {
            let i = num_cells_blade_half;
            let x_i = ss_edge.point_coord(i);
            let x_im1 = ss_edge.point_coord(i - 1);

            let x_xi = -x_im1 + x_i;

            let n = (1.0 / x_xi.abs()) * Vec2d(-x_xi.1, x_xi.0);

            ss_outer[i] = x_i + d * n;
        }

        let ss_outer_edge = EdgeView::new(Edge::new_fixed(
            "Suction_Side_Outer_Edge".to_string(),
            ss_outer,
        ));

        let mut ps_outer = vec![Vec2d(0.0, 0.0); num_cells_blade_half + 1];

        ps_outer[0] = ss_outer_edge.point_coord(ss_outer_edge.start);
        ps_outer[num_cells_blade_half] = ss_outer_edge.point_coord(ss_outer_edge.end);

        // handle all points except the first and last
        for i in 1..num_cells_blade_half {
            let x_i = ps_edge.point_coord(i);
            let x_im1 = ps_edge.point_coord(i - 1);
            let x_ip1 = ps_edge.point_coord(i + 1);

            let x_xi = 0.5 * (x_ip1 - x_im1);

            let n = (1.0 / x_xi.abs()) * Vec2d(x_xi.1, -x_xi.0);

            ps_outer[i] = x_i + d * n;
        }

        let ps_outer_edge = EdgeView::new(Edge::new_fixed(
            "Pressure_Side_Outer_Edge".to_string(),
            ps_outer,
        ));

        // split blade distribution

        let (ps_edge_in_middle, ps_edge_rest) =
            ps_outer_edge.split_at(num_cells_middle_blocks_on_blade_half);
        let (ps_edge_in_lower, ps_edge_rest) =
            ps_edge_rest.split_at(num_cells_next_to_middle_blocks_on_blade);
        let (ps_edge_pll1, ps_edge_ex_middle) =
            ps_edge_rest.split_at(ps_edge_rest.len() - 1 - num_cells_middle_blocks_on_blade_half);

        let (ss_edge_in_middle, ss_edge_rest) =
            ss_outer_edge.split_at(num_cells_middle_blocks_on_blade_half);
        let (ss_edge_in_ss, ss_edge_rest) = ss_edge_rest.split_at(
            ss_edge_rest.len()
                - 1
                - num_cells_middle_blocks_on_blade_half
                - num_cells_next_to_middle_blocks_on_blade,
        );
        let (ss_edge_ex_ss, ss_edge_ex_middle) =
            ss_edge_rest.split_at(ss_edge_rest.len() - 1 - num_cells_middle_blocks_on_blade_half);

        // TODO add approximation of the leading and trailing edge by approximating
        // the chamber line and computing its intersection with the profile.
        // For now, the values of the given profiles are used.
        let leading_edge = ps_spline.interpolate_val(0.0);
        let trailing_edge = ps_spline.interpolate_val(1.0);

        let mut mesh = Mesh::new();

        let row_prefix = "row_01_";

        let viscous = true;
        let o_grid_blocks: usize = 2;

        if viscous {
            {
                // |          /
                // |         |          ^
                // |        |           |
                // |-------< x_00      i_min
                //         LE

                let block_name = "blade_ss";

                let x_00 = ss_edge.point_coord(ss_edge.start);
                let x_10 = ss_edge.point_coord(ss_edge.end);
                let x_01 = ss_outer_edge.point_coord(ss_outer_edge.start);
                let x_11 = ss_outer_edge.point_coord(ss_outer_edge.end);

                let block = Block2d::new(
                    row_prefix.to_owned() + block_name,
                    vec![Box::new(ss_edge)],
                    vec![
                        Box::new(ss_edge_in_middle),
                        Box::new(ss_edge_in_ss),
                        Box::new(ss_edge_ex_ss),
                        Box::new(ss_edge_ex_middle),
                    ],
                    vec![Box::new(Segment::new(
                        num_cells_away_from_blade + 1,
                        SingleHyperbolicTangentClustering::new(0.01),
                        Line2d::new(x_00, x_01),
                    ))],
                    vec![Box::new(Segment::new(
                        num_cells_away_from_blade + 1,
                        SingleHyperbolicTangentClustering::new(0.01),
                        Line2d::new(x_10, x_11),
                    ))],
                );

                mesh.add_block(block);

                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    0,
                    EdgeIndex::IMin,
                    0..1,
                )));

                // TODO remove
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    0,
                    EdgeIndex::IMax,
                    1..3,
                )));
            }

            {
                //         LE
                // |-------< x_00      i_min
                // |        |            |
                // |         |           v
                // |          /

                let block_name = "blade_ps";

                let x_00 = ps_edge.point_coord(ps_edge.start);
                let x_10 = ps_edge.point_coord(ps_edge.end);
                let x_01 = ps_outer_edge.point_coord(ps_outer_edge.start);
                let x_11 = ps_outer_edge.point_coord(ps_outer_edge.end);

                let block = Block2d::new(
                    row_prefix.to_owned() + block_name,
                    vec![Box::new(ps_edge)],
                    vec![
                        Box::new(ps_edge_in_middle),
                        Box::new(ps_edge_in_lower),
                        Box::new(ps_edge_pll1),
                        Box::new(ps_edge_ex_middle),
                    ],
                    vec![Box::new(Segment::new(
                        num_cells_away_from_blade + 1,
                        SingleHyperbolicTangentClustering::new(0.01),
                        Line2d::new(x_00, x_01),
                    ))],
                    vec![Box::new(Segment::new(
                        num_cells_away_from_blade + 1,
                        SingleHyperbolicTangentClustering::new(0.01),
                        Line2d::new(x_10, x_11),
                    ))],
                );

                mesh.add_block(block);

                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    1,
                    EdgeIndex::IMin,
                    0..1,
                )));

                // // TODO remove
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     1,
                //     EdgeIndex::IMax,
                //     3..,
                // )));

                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMax, 0..1),
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMax, 0..1),
                        ),
                    )));

                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMin, 0..1),
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMin, 0..1),
                        ),
                    )));
            }
        }

        {
            // x_01
            // |---------- x_00
            // |         |
            // |        |           |
            // |       < LE       i_min
            // |        |           |
            // |         |          V
            // |---------- x_10
            // x_11

            let block_name = "inlet_middle";

            let edge_i_min_0 = mesh.blocks[0].edge_segment(EdgeIndex::IMax, 0);
            let edge_i_min_1 = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 0);

            let x_00 = *edge_i_min_0.x.last().unwrap();
            let x_10 = *edge_i_min_1.x.last().unwrap();

            // TODO remove hard coded position for block
            let x_01 = x_00 - Vec2d(0.02, -0.001);

            // TODO remove hard coded position for block
            let x_11 = x_10 - Vec2d(0.02, 0.02);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(edge_i_min_0.reverse()), Box::new(edge_i_min_1)],
                vec![Box::new(Segment::new(
                    num_cells_middle_blocks_on_blade_half * 2 + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    num_cells_away_from_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    num_cells_away_from_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            if !viscous {
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    0,
                    EdgeIndex::IMin,
                    0..=1,
                )));
            } else {
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 0..1),
                            BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMin, 0..1).reverse(),
                        ),
                    )));
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 0..1),
                            BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMin, 1..2),
                        ),
                    )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    2,
                    EdgeIndex::IMax,
                    0..1,
                )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    2,
                    EdgeIndex::JMin,
                    0..1,
                )));
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     2,
                //     EdgeIndex::JMax,
                //     0..1,
                // )));
            }
        }

        {
            let block_name = "inlet_ps";

            // copy j max edge of in_middle block (id: 0 or 2)
            let edge_j_min = mesh.blocks[0 + o_grid_blocks].edge_data(EdgeIndex::JMax);

            let edge_i_min = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 1);

            let x_01 = *edge_j_min.x.last().unwrap();
            let x_10 = *edge_i_min.x.last().unwrap();

            // on periodic bc
            let x_11 = leading_edge + Vec2d(0.0, -0.5 * pitch);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(edge_i_min)],
                vec![Box::new(Segment::new(
                    num_cells_next_to_middle_blocks_on_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(edge_j_min)],
                vec![Box::new(Segment::new(
                    num_cells_away_from_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 0 + o_grid_blocks, EdgeIndex::JMax, 0..1),
                        BlockBoundaryRange::new(&mesh, 1 + o_grid_blocks, EdgeIndex::JMin, 0..1),
                    ),
                )));

            if !viscous {
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    1,
                    EdgeIndex::IMin,
                    0..1,
                )));
            } else {
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 1..2),
                            BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMin, ..),
                        ),
                    )));

                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     3,
                //     EdgeIndex::IMin,
                //     ..,
                // )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    3,
                    EdgeIndex::IMax,
                    ..,
                )));
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     3,
                //     EdgeIndex::JMin,
                //     ..,
                // )));
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     3,
                //     EdgeIndex::JMax,
                //     ..,
                // )));
            }
        }

        {
            let block_name = "exit_ps";

            // copy j max edge of in_lower block (id: 1)
            let edge_j_min = mesh.blocks[1 + o_grid_blocks].edge_data(EdgeIndex::JMax);

            let edge_i_min_00 = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 2);

            let x_01 = *edge_j_min.x.last().unwrap();
            let x_blade_end = *edge_i_min_00.x.last().unwrap();

            // TODO remove hard coded coordinate
            let x_10 = x_blade_end + Vec2d(0.02, -0.025);

            // on periodic bc
            let x_11 = trailing_edge + Vec2d(0.0, -0.5 * pitch);

            let len_i_max = edge_i_min_00.len() + num_cells_away_from_blade;

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![
                    Box::new(edge_i_min_00),
                    Box::new(Segment::new(
                        num_cells_away_from_blade + 1,
                        UniformClustering::new(),
                        Line2d::new(x_blade_end, x_10),
                    )),
                ],
                vec![Box::new(Segment::new(
                    len_i_max,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(edge_j_min)],
                vec![Box::new(Segment::new(
                    num_cells_away_from_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 1 + o_grid_blocks, EdgeIndex::JMax, 0..1),
                        BlockBoundaryRange::new(&mesh, 2 + o_grid_blocks, EdgeIndex::JMin, 0..1),
                    ),
                )));

            if !viscous {
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    2,
                    EdgeIndex::IMin,
                    0..1,
                )));
            } else {
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 2..3),
                            BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMin, 0..1),
                        ),
                    )));

                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     4,
                //     EdgeIndex::IMin,
                //     1..,
                // )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    4,
                    EdgeIndex::IMax,
                    ..,
                )));
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     4,
                //     EdgeIndex::JMin,
                //     ..,
                // )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    4,
                    EdgeIndex::JMax,
                    ..,
                )));
            }
        }

        {
            let block_name = "exit_middle";

            // copy last segment of i min edge of pll1 block (id: 2)
            let edge_j_min = mesh.blocks[2 + o_grid_blocks].edge_segment(EdgeIndex::IMin, 1);

            // ps_edge_ex_middle
            let edge_i_min_00 = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 3);
            // ss_edge_ex_middle
            let edge_i_min_01 = mesh.blocks[0].edge_segment(EdgeIndex::IMax, 3);

            let x_01 = *edge_j_min.x.last().unwrap();
            let x_10 = *edge_i_min_01.x.first().unwrap();

            // TODO remove hard coded coordinate
            let x_11 = x_10 + Vec2d(0.007, 0.001);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(edge_i_min_00), Box::new(edge_i_min_01.reverse())],
                vec![Box::new(Segment::new(
                    num_cells_middle_blocks_on_blade_half * 2 + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(edge_j_min)],
                vec![Box::new(Segment::new(
                    num_cells_away_from_blade + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 2 + o_grid_blocks, EdgeIndex::IMin, 1..2),
                        BlockBoundaryRange::new(&mesh, 3 + o_grid_blocks, EdgeIndex::JMin, 0..1),
                    ),
                )));

            if !viscous {
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    3,
                    EdgeIndex::IMin,
                    0..2,
                )));
            } else {
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 3..=3),
                            BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 0..=0),
                        ),
                    )));
                mesh.edges
                    .push(BlockBoundary::Connection(BlockConnection::new(
                        &mesh,
                        (
                            BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 3..=3).reverse(),
                            BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 1..=1),
                        ),
                    )));

                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     5,
                //     EdgeIndex::IMin,
                //     ..,
                // )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    5,
                    EdgeIndex::IMax,
                    ..,
                )));
                // mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                //     &mesh,
                //     5,
                //     EdgeIndex::JMin,
                //     ..,
                // )));
                mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
                    &mesh,
                    5,
                    EdgeIndex::JMax,
                    ..,
                )));
            }
        }

        // {
        //     let block_name = "exit_ss";
        //
        //     let edge_j_min = mesh.blocks.last().unwrap().edge_data(EdgeIndex::JMax);
        //
        //     let x_01 = *edge_j_min.x.last().unwrap();
        //     let x_10 = ss_edge_ex_ss.point_coord(ss_edge_ex_ss.start);
        //     let x_11 = trailing_edge + Vec2d(0.0, 0.5 * pitch);
        //
        //     let block = Block2d::new(
        //         row_prefix.to_owned() + block_name,
        //         vec![Box::new(ss_edge_ex_ss.rev())],
        //         vec![Box::new(Segment::new(
        //             ss_edge_ex_ss.len(),
        //             UniformClustering::new(),
        //             Line2d::new(x_01, x_11),
        //         ))],
        //         vec![Box::new(edge_j_min)],
        //         vec![Box::new(Segment::new(
        //             num_cells_away_from_blade + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_10, x_11),
        //         ))],
        //     );
        //
        //     mesh.add_block(block);
        //
        //     mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        //         &mesh,
        //         4,
        //         EdgeIndex::IMin,
        //         0..1,
        //     )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 3, EdgeIndex::JMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMin, 0..1),
        //             ),
        //         )));
        // }
        //
        // {
        //     let block_name = "inlet_ss";
        //
        //     let edge_j_min = mesh.blocks.last().unwrap().edge_data(EdgeIndex::JMax);
        //     let edge_i_min_last_segment = mesh.blocks[0].edge_data(EdgeIndex::JMin);
        //
        //     let x_11 = leading_edge + Vec2d(0.0, 0.5 * pitch);
        //     let x_10 = *edge_i_min_last_segment.x.last().unwrap();
        //     let x_01 = *edge_j_min.x.last().unwrap();
        //
        //     let num_points_i_max = edge_i_min_last_segment.len() + ss_edge_in_ss.len() - 1;
        //
        //     let block = Block2d::new(
        //         row_prefix.to_owned() + block_name,
        //         vec![
        //             Box::new(ss_edge_in_ss.rev()),
        //             Box::new(edge_i_min_last_segment),
        //         ],
        //         vec![Box::new(Segment::new(
        //             num_points_i_max,
        //             UniformClustering::new(),
        //             Line2d::new(x_01, x_11),
        //         ))],
        //         vec![Box::new(edge_j_min)],
        //         vec![Box::new(Segment::new(
        //             num_cells_away_from_blade + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_10, x_11),
        //         ))],
        //     );
        //
        //     mesh.add_block(block);
        //
        //     mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        //         &mesh,
        //         5,
        //         EdgeIndex::IMin,
        //         0..1,
        //     )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMin, 0..1),
        //             ),
        //         )));
        //     mesh.edges.push(BlockBoundary::PeriodicConnection(
        //         PeriodicBlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMax, 0..1).reverse(),
        //             ),
        //             Vec2d(0.0, pitch),
        //         ),
        //     ));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMin, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 1..2),
        //             ),
        //         )));
        // }
        //
        // {
        //     let block_name = "inlet";
        //
        //     // copy from inlet_lower
        //     let edge_j_max_seg_0 = mesh.blocks[1].edge_data(EdgeIndex::IMax);
        //     // copy from inlet_middle
        //     let edge_j_max_seg_1 = mesh.blocks[0].edge_data(EdgeIndex::IMax);
        //     // copy from inlet_ss
        //     let edge_j_max_seg_2 = mesh.blocks[5].edge_data(EdgeIndex::JMax);
        //
        //     let num_points_j =
        //         edge_j_max_seg_0.len() + edge_j_max_seg_1.len() + edge_j_max_seg_2.len() - 2;
        //
        //     let x_10 = *edge_j_max_seg_0.x.last().unwrap();
        //     let x_11 = *edge_j_max_seg_2.x.last().unwrap();
        //     let x_00 = x_10 + Vec2d(-0.05, 0.0);
        //     let x_01 = x_00 + Vec2d(0.0, pitch);
        //
        //     let block = Block2d::new(
        //         row_prefix.to_owned() + block_name,
        //         vec![Box::new(Segment::new(
        //             num_cells_inlet + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_00, x_10),
        //         ))],
        //         vec![Box::new(Segment::new(
        //             num_cells_inlet + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_01, x_11),
        //         ))],
        //         vec![Box::new(Segment::new(
        //             num_points_j,
        //             UniformClustering::new(),
        //             Line2d::new(x_00, x_01),
        //         ))],
        //         vec![
        //             Box::new(edge_j_max_seg_0.reverse()),
        //             Box::new(edge_j_max_seg_1.reverse()),
        //             Box::new(edge_j_max_seg_2),
        //         ],
        //     );
        //
        //     mesh.add_block(block);
        //
        //     mesh.edges
        //         .push(BlockBoundary::Inlet(BlockBoundaryRange::new(
        //             &mesh,
        //             6,
        //             EdgeIndex::JMin,
        //             0..1,
        //         )));
        //     mesh.edges.push(BlockBoundary::PeriodicConnection(
        //         PeriodicBlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMin, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMax, 0..1),
        //             ),
        //             Vec2d(0.0, pitch),
        //         ),
        //     ));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 0..1).reverse(),
        //                 BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 0..1),
        //             ),
        //         )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 0..1).reverse(),
        //                 BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 1..2),
        //             ),
        //         )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 2..3),
        //             ),
        //         )));
        // }
        //
        // {
        //     let block_name = "exit";
        //
        //     // copy from exit_lower
        //     let edge_j_min_seg_0 = mesh.blocks[2].edge_data(EdgeIndex::JMax);
        //     // copy from exit_middle
        //     let edge_j_min_seg_1 = mesh.blocks[3].edge_data(EdgeIndex::IMax);
        //     // copy from exit_ss
        //     let edge_j_min_seg_2 = mesh.blocks[4].edge_data(EdgeIndex::IMax);
        //
        //     let num_points_j =
        //         edge_j_min_seg_0.len() + edge_j_min_seg_1.len() + edge_j_min_seg_2.len() - 2;
        //
        //     let x_00 = *edge_j_min_seg_0.x.last().unwrap();
        //     let x_01 = *edge_j_min_seg_2.x.last().unwrap();
        //     let x_10 = x_00 + Vec2d(0.05, 0.0);
        //     let x_11 = x_10 + Vec2d(0.0, pitch);
        //
        //     let block = Block2d::new(
        //         row_prefix.to_owned() + block_name,
        //         vec![Box::new(Segment::new(
        //             num_cells_outlet + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_00, x_10),
        //         ))],
        //         vec![Box::new(Segment::new(
        //             num_cells_outlet + 1,
        //             UniformClustering::new(),
        //             Line2d::new(x_01, x_11),
        //         ))],
        //         vec![
        //             Box::new(edge_j_min_seg_0.reverse()),
        //             Box::new(edge_j_min_seg_1),
        //             Box::new(edge_j_min_seg_2),
        //         ],
        //         vec![Box::new(Segment::new(
        //             num_points_j,
        //             UniformClustering::new(),
        //             Line2d::new(x_10, x_11),
        //         ))],
        //     );
        //
        //     mesh.add_block(block);
        //
        //     mesh.edges
        //         .push(BlockBoundary::Outlet(BlockBoundaryRange::new(
        //             &mesh,
        //             7,
        //             EdgeIndex::JMax,
        //             0..1,
        //         )));
        //     mesh.edges.push(BlockBoundary::PeriodicConnection(
        //         PeriodicBlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMin, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMax, 0..1),
        //             ),
        //             Vec2d(0.0, pitch),
        //         ),
        //     ));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 2, EdgeIndex::JMax, 0..1).reverse(),
        //                 BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 0..1),
        //             ),
        //         )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 1..2),
        //             ),
        //         )));
        //     mesh.edges
        //         .push(BlockBoundary::Connection(BlockConnection::new(
        //             &mesh,
        //             (
        //                 BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMax, 0..1),
        //                 BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 2..3),
        //             ),
        //         )));
        // }

        mesh.save("turbomesh_linear.cgns").unwrap();

        mesh.smooth(&self.smoothing);

        mesh.save("turbomesh.cgns").unwrap();

        // // plot blocking
        // let ps_spline_int =
        //     ps_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
        // let ss_spline_int =
        //     ss_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
        //
        // let ps_clustering_int = &ps_edge.edge().coords;
        // let ss_clustering_int = &ss_edge.edge().coords;
        //
        // let root = BitMapBackend::new("blocking.png", (1900, 1200)).into_drawing_area();
        // root.fill(&WHITE).unwrap();
        // let mut chart = ChartBuilder::on(&root)
        //     // .caption("y=x^2", ("sans-serif", 50).into_font())
        //     .margin(5)
        //     .x_label_area_size(30)
        //     .y_label_area_size(30)
        //     .build_cartesian_2d(0.98f32..1.2f32, -0.1f32..0.1f32)
        //     .unwrap();
        //
        // chart.configure_mesh().draw().unwrap();
        //
        // chart
        //     .draw_series(LineSeries::new(
        //         ps_spline_int.iter().map(|v| (v[0] as f32, v[1] as f32)),
        //         &RED,
        //     ))
        //     .unwrap()
        //     .label("PS");
        // chart
        //     .draw_series(LineSeries::new(
        //         ss_spline_int.iter().map(|v| (v[0] as f32, v[1] as f32)),
        //         &BLUE,
        //     ))
        //     .unwrap()
        //     .label("SS");
        //
        // chart
        //     .draw_series(
        //         ps_clustering_int
        //             .iter()
        //             .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLACK.filled())),
        //     )
        //     .unwrap();
        // chart
        //     .draw_series(
        //         ss_clustering_int
        //             .iter()
        //             .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLACK.filled())),
        //     )
        //     .unwrap();
        //
        // for block in mesh.blocks.iter() {
        //     chart
        //         .draw_series(
        //             block
        //                 .coords
        //                 .iter()
        //                 .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLUE.filled())),
        //         )
        //         .unwrap();
        // }
        //
        // chart
        //     .configure_series_labels()
        //     .background_style(&WHITE.mix(0.8))
        //     .border_style(&BLACK)
        //     .draw()
        //     .unwrap();
        //
        // root.present().unwrap();
    }
}
