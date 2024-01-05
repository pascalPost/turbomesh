// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::clustering::{RobertsClustering, SingleHyperbolicTangentClustering, UniformClustering};
use crate::geometry::Line2d;
use crate::smoothing::SmoothingMethod;
use crate::templates::helper::blade_profile;
use crate::types::{
    BlockBoundary, BlockBoundaryRange, BlockConnection, Edge, EdgeIndex, EdgeView,
    PeriodicBlockConnection,
};
use crate::{Block2d, Mesh, Segment, Vec2d};
use plotters::prelude::*;
use serde::Deserialize;

#[derive(Deserialize, Debug)]
#[serde(tag = "type")]
enum BladeClusteringFunction {
    RobertsClustering(RobertsClustering),
}

#[derive(Deserialize, Debug)]
pub struct O4HTemplate {
    ps_csv_path: std::path::PathBuf,
    ss_csv_path: std::path::PathBuf,
    pitch: f64,
    blade_clustering_function: BladeClusteringFunction,

    num_cells_ogrid: usize,
    num_cells_in_half_i: usize,
    num_cells_in_j: usize,
    num_cells_out_half_i: usize,
    num_cells_out_j: usize,
    num_cells_down_i: usize,
    num_cells_up_i: usize,
    num_cells_middle_j: usize,

    // TODO can be made optional if computed automatically based on the average size
    num_cells_upstream_i: usize,
    num_cells_downstream_i: usize,

    smoothing: SmoothingMethod,
}

impl O4HTemplate {
    /// appends the given root path to the file paths of the template
    pub fn append_root_path(&mut self, root_path: &std::path::Path) {
        self.ps_csv_path = root_path.join(&self.ps_csv_path);
        self.ss_csv_path = root_path.join(&self.ss_csv_path);
    }

    /// runs the template and returns the resulting geometry and mesh
    pub fn run(&self) {
        let ps_csv_path = self.ps_csv_path.to_str().unwrap();
        let ss_csv_path = self.ss_csv_path.to_str().unwrap();

        let (ps_spline, ss_spline) = blade_profile(ps_csv_path, ss_csv_path)
            .expect("Error in profile input and interpolation");

        // blade discretization on which the 2d blocking is based
        let blade_clustering_function = match self.blade_clustering_function {
            BladeClusteringFunction::RobertsClustering(roberts_clustering) => roberts_clustering,
        };

        // number of cells on the ss and ps of the blade
        let num_cells_blade_half =
            self.num_cells_in_half_i + self.num_cells_middle_j + self.num_cells_out_half_i;

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

        let (ps_edge_in, ps_edge_rest) = ps_outer_edge.split_at(self.num_cells_in_half_i);
        let (ps_edge_middle, ps_edge_out) = ps_edge_rest.split_at(self.num_cells_middle_j);

        let (ss_edge_in, ss_edge_rest) = ss_outer_edge.split_at(self.num_cells_in_half_i);
        let (ss_edge_middle, ss_edge_out) = ss_edge_rest.split_at(self.num_cells_middle_j);

        // TODO add approximation of the leading and trailing edge by approximating
        // the chamber line and computing its intersection with the profile.
        // For now, the values of the given profiles are used.
        let leading_edge = ps_spline.interpolate_val(0.0);
        let trailing_edge = ps_spline.interpolate_val(1.0);

        let mut mesh = Mesh::new();

        let row_prefix = "row_01_";

        {
            // |          /
            // |         |          ^
            // |        |           |
            // |-------< x_00      i_min
            //         LE

            let block_name = "ss";

            let x_00 = ss_edge.point_coord(ss_edge.start);
            let x_10 = ss_edge.point_coord(ss_edge.end);
            let x_01 = ss_outer_edge.point_coord(ss_outer_edge.start);
            let x_11 = ss_outer_edge.point_coord(ss_outer_edge.end);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(ss_edge)],
                vec![
                    Box::new(ss_edge_in),
                    Box::new(ss_edge_middle),
                    Box::new(ss_edge_out),
                ],
                vec![Box::new(Segment::new(
                    self.num_cells_ogrid + 1,
                    SingleHyperbolicTangentClustering::new(0.01),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_ogrid + 1,
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
        }

        {
            //         LE
            // |-------< x_00      i_min
            // |        |            |
            // |         |           v
            // |          /

            let block_name = "ps";

            let x_00 = ps_edge.point_coord(ps_edge.start);
            let x_10 = ps_edge.point_coord(ps_edge.end);
            let x_01 = ps_outer_edge.point_coord(ps_outer_edge.start);
            let x_11 = ps_outer_edge.point_coord(ps_outer_edge.end);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(ps_edge)],
                vec![
                    Box::new(ps_edge_in),
                    Box::new(ps_edge_middle),
                    Box::new(ps_edge_out),
                ],
                vec![Box::new(Segment::new(
                    self.num_cells_ogrid + 1,
                    SingleHyperbolicTangentClustering::new(0.01),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_ogrid + 1,
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

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMax, ..),
                        BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMax, ..),
                    ),
                )));

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMin, ..),
                        BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMin, ..),
                    ),
                )));
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

            let block_name = "in";

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
                    self.num_cells_in_half_i * 2 + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_in_j + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_in_j + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 0..=0),
                        BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMin, 0..=0).reverse(),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 0..=0),
                        BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMin, 1..=1),
                    ),
                )));
        }

        {
            let block_name = "out";

            let edge_i_min_0 = mesh.blocks[0].edge_segment(EdgeIndex::IMax, 2);
            let edge_i_min_1 = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 2);

            let x_00 = *edge_i_min_0.x.first().unwrap();
            let x_10 = *edge_i_min_1.x.first().unwrap();

            // TODO remove hard coded coordinate
            let x_01 = x_00 + Vec2d(0.02, 0.001);

            // TODO remove hard coded coordinate
            let x_11 = x_10 + Vec2d(0.015, -0.04);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(edge_i_min_0), Box::new(edge_i_min_1.reverse())],
                vec![Box::new(Segment::new(
                    self.num_cells_out_half_i * 2 + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_out_j + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_out_j + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 2..=2),
                        BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMin, 0..=0),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 2..=2),
                        BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMin, 1..=1).reverse(),
                    ),
                )));
        }

        {
            let block_name = "down";

            let edge_i_min_0 = mesh.blocks[2].edge_segment(EdgeIndex::JMax, 0);
            let edge_i_min_1 = mesh.blocks[1].edge_segment(EdgeIndex::IMax, 1);
            let edge_i_min_2 = mesh.blocks[3].edge_segment(EdgeIndex::JMax, 0);

            let len_edge_i = edge_i_min_0.len() + edge_i_min_1.len() + edge_i_min_2.len() - 3;

            let x_00 = *edge_i_min_0.x.last().unwrap();
            let x_10 = *edge_i_min_2.x.last().unwrap();

            // on periodic bc
            let x_01 = leading_edge - Vec2d(0.0, 0.5 * self.pitch);
            let x_11 = trailing_edge - Vec2d(0.0, 0.5 * self.pitch);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![
                    Box::new(edge_i_min_0.reverse()),
                    Box::new(edge_i_min_1),
                    Box::new(edge_i_min_2),
                ],
                vec![Box::new(Segment::new(
                    len_edge_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_down_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_down_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMin, 0..=0),
                        BlockBoundaryRange::new(&mesh, 2, EdgeIndex::JMax, ..).reverse(),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMin, 1..=1),
                        BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 1..=1),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMin, 2..=2),
                        BlockBoundaryRange::new(&mesh, 3, EdgeIndex::JMax, ..),
                    ),
                )));
        }

        {
            let block_name = "up";

            let edge_i_min_0 = mesh.blocks[2].edge_segment(EdgeIndex::JMin, 0);
            let edge_i_min_1 = mesh.blocks[0].edge_segment(EdgeIndex::IMax, 1);
            let edge_i_min_2 = mesh.blocks[3].edge_segment(EdgeIndex::JMin, 0);

            let len_edge_i = edge_i_min_0.len() + edge_i_min_1.len() + edge_i_min_2.len() - 3;

            let x_00 = *edge_i_min_0.x.last().unwrap();
            let x_10 = *edge_i_min_2.x.last().unwrap();

            // on periodic bc
            let x_01 = leading_edge + Vec2d(0.0, 0.5 * self.pitch);
            let x_11 = trailing_edge + Vec2d(0.0, 0.5 * self.pitch);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![
                    Box::new(edge_i_min_0.reverse()),
                    Box::new(edge_i_min_1),
                    Box::new(edge_i_min_2),
                ],
                vec![Box::new(Segment::new(
                    len_edge_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_down_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_down_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 0..=0),
                        BlockBoundaryRange::new(&mesh, 2, EdgeIndex::JMin, ..).reverse(),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 1..=1),
                        BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 1..=1),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 2..=2),
                        BlockBoundaryRange::new(&mesh, 3, EdgeIndex::JMin, ..),
                    ),
                )));
            mesh.edges.push(BlockBoundary::PeriodicConnection(
                PeriodicBlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMax, ..),
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMax, ..),
                    ),
                    Vec2d(0.0, self.pitch),
                ),
            ));
        }

        {
            let block_name = "upstream";

            let edge_j_max_seg_0 = mesh.blocks[4].edge_data(EdgeIndex::JMin);
            let edge_j_max_seg_1 = mesh.blocks[2].edge_data(EdgeIndex::IMax);
            let edge_j_max_seg_2 = mesh.blocks[5].edge_data(EdgeIndex::JMin);

            let num_points_j =
                edge_j_max_seg_0.len() + edge_j_max_seg_1.len() + edge_j_max_seg_2.len() - 2;

            let x_10 = *edge_j_max_seg_0.x.last().unwrap();
            let x_11 = *edge_j_max_seg_2.x.last().unwrap();
            let x_00 = x_10 + Vec2d(-0.05, 0.0);
            let x_01 = x_00 + Vec2d(0.0, self.pitch);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(Segment::new(
                    self.num_cells_upstream_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_10),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_upstream_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![Box::new(Segment::new(
                    num_points_j,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_01),
                ))],
                vec![
                    Box::new(edge_j_max_seg_0.reverse()),
                    Box::new(edge_j_max_seg_1.reverse()),
                    Box::new(edge_j_max_seg_2),
                ],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Inlet(BlockBoundaryRange::new(
                    &mesh,
                    6,
                    EdgeIndex::JMin,
                    ..,
                )));
            mesh.edges.push(BlockBoundary::PeriodicConnection(
                PeriodicBlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMin, ..),
                        BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMax, ..),
                    ),
                    Vec2d(0.0, self.pitch),
                ),
            ));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMin, ..).reverse(),
                        BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 0..=0),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMax, ..).reverse(),
                        BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 1..=1),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMin, ..),
                        BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 2..=2),
                    ),
                )));
        }

        {
            let block_name = "downstream";

            let edge_j_min_seg_0 = mesh.blocks[4].edge_data(EdgeIndex::JMax);
            let edge_j_min_seg_1 = mesh.blocks[3].edge_data(EdgeIndex::IMax);
            let edge_j_min_seg_2 = mesh.blocks[5].edge_data(EdgeIndex::JMax);

            let num_points_j =
                edge_j_min_seg_0.len() + edge_j_min_seg_1.len() + edge_j_min_seg_2.len() - 2;

            let x_00 = *edge_j_min_seg_0.x.last().unwrap();
            let x_01 = *edge_j_min_seg_2.x.last().unwrap();
            let x_10 = x_00 + Vec2d(0.05, 0.0);
            let x_11 = x_10 + Vec2d(0.0, self.pitch);

            let block = Block2d::new(
                row_prefix.to_owned() + block_name,
                vec![Box::new(Segment::new(
                    self.num_cells_downstream_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_00, x_10),
                ))],
                vec![Box::new(Segment::new(
                    self.num_cells_downstream_i + 1,
                    UniformClustering::new(),
                    Line2d::new(x_01, x_11),
                ))],
                vec![
                    Box::new(edge_j_min_seg_0.reverse()),
                    Box::new(edge_j_min_seg_1.reverse()),
                    Box::new(edge_j_min_seg_2),
                ],
                vec![Box::new(Segment::new(
                    num_points_j,
                    UniformClustering::new(),
                    Line2d::new(x_10, x_11),
                ))],
            );

            mesh.add_block(block);

            mesh.edges
                .push(BlockBoundary::Outlet(BlockBoundaryRange::new(
                    &mesh,
                    7,
                    EdgeIndex::JMax,
                    ..,
                )));
            mesh.edges.push(BlockBoundary::PeriodicConnection(
                PeriodicBlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMin, ..),
                        BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMax, ..),
                    ),
                    Vec2d(0.0, self.pitch),
                ),
            ));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMax, ..).reverse(),
                        BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 0..=0),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMax, ..).reverse(),
                        BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 1..=1),
                    ),
                )));
            mesh.edges
                .push(BlockBoundary::Connection(BlockConnection::new(
                    &mesh,
                    (
                        BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMax, ..),
                        BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 2..=2),
                    ),
                )));
        }

        mesh.save("turbomesh_linear.cgns").unwrap();

        mesh.smooth(&self.smoothing);

        mesh.save("turbomesh.cgns").unwrap();
    }
}
