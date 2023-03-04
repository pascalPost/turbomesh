// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::clustering::{RobertsClustering, UniformClustering};
use crate::geometry::{Line2d, Spline};
use crate::interpolation::FittingSpline;
use crate::types::{BlockBoundary, BlockBoundaryRange, BlockConnection, Edge, EdgeIndex, EdgeView};
use crate::{Block2d, Geometry, Mesh, Segment, Vec2d};
use ndarray::Array;
use plotters::prelude::*;
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

pub fn run_turbine_template(ps_csv_path: &str, ss_csv_path: &str) -> (Geometry, Mesh) {
    // let mut geometry = Geometry {};

    // TODO move all of these parameters to function input
    let pitch = 88.36 * 1e-3; //m
    let num_cells_blade_half: usize = 60;
    let num_cells_middle_blocks_on_blade_half: usize = 5;
    let num_cells_next_to_middle_blocks_on_blade: usize = 20;
    let num_cells_away_from_blade: usize = 10;

    // TODO can be computed automatically based on the average size
    let num_cells_inlet = 20;
    let num_cells_outlet = 20;

    let (ps_spline, ss_spline) =
        blade_profile(ps_csv_path, ss_csv_path).expect("Error in profile input and interpolation");

    // blade discretization on which the 2d blocking is based
    let blade_clustering_function = RobertsClustering::new(0.5, 1.01);

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

    // split blade distribution

    let (ps_edge_in_middle, ps_edge_rest) = ps_edge.split_at(num_cells_middle_blocks_on_blade_half);
    let (ps_edge_in_lower, ps_edge_rest) =
        ps_edge_rest.split_at(num_cells_next_to_middle_blocks_on_blade);
    let (ps_edge_pll1, ps_edge_ex_middle) =
        ps_edge_rest.split_at(ps_edge_rest.len() - 1 - num_cells_middle_blocks_on_blade_half);

    let (ss_edge_in_middle, ss_edge_rest) = ss_edge.split_at(num_cells_middle_blocks_on_blade_half);
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

        let x_00 = ss_edge_in_middle.point_coord(ss_edge_in_middle.end);
        let x_10 = ps_edge_in_middle.point_coord(ps_edge_in_middle.end);

        // TODO remove hard coded position for block
        let x_01 = x_00 - Vec2d(0.007, -0.025);

        // TODO remove hard coded position for block
        let x_11 = x_10 - Vec2d(0.007, 0.001);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![
                Box::new(ss_edge_in_middle.rev()),
                Box::new(ps_edge_in_middle),
            ],
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                0,
                EdgeIndex::IMin,
                0..=1,
            )));
    }

    {
        let block_name = "inlet_ps";

        // copy j max edge of in_middle block (id: 0)
        let edge_j_min = mesh.blocks[0].edge_data(EdgeIndex::JMax);

        let x_01 = *edge_j_min.x.last().unwrap();
        let x_10 = ps_edge_in_lower.point_coord(ps_edge_in_lower.end);

        // on periodic bc
        let x_11 = leading_edge + Vec2d(0.0, -0.5 * pitch);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![Box::new(ps_edge_in_lower)],
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                1,
                EdgeIndex::IMin,
                0..1,
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMin, 0..1),
                ),
            )));
    }

    {
        let block_name = "exit_ss";

        // copy j max edge of in_lower block (id: 1)
        let edge_j_min = mesh.blocks[1].edge_data(EdgeIndex::JMax);

        let x_01 = *edge_j_min.x.last().unwrap();
        let x_blade_end = ps_edge_pll1.point_coord(ps_edge_pll1.end);

        // TODO remove hard coded coordinate
        let x_10 = x_blade_end + Vec2d(0.007, -0.025);

        // on periodic bc
        let x_11 = trailing_edge + Vec2d(0.0, -0.5 * pitch);

        let len_i_max = ps_edge_pll1.len() + num_cells_away_from_blade;

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![
                Box::new(ps_edge_pll1),
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                2,
                EdgeIndex::IMin,
                0..1,
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 1, EdgeIndex::JMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 2, EdgeIndex::JMin, 0..1),
                ),
            )));
    }

    {
        let block_name = "exit_middle";

        // copy last segment of i min edge of pll1 block (id: 2)
        let edge_j_min = mesh.blocks[2].edge_segment(EdgeIndex::IMin, 1);

        let x_01 = *edge_j_min.x.last().unwrap();
        let x_10 = ss_edge_ex_middle.point_coord(ss_edge_ex_middle.start);

        // TODO remove hard coded coordinate
        let x_11 = x_10 + Vec2d(0.007, 0.001);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![
                Box::new(ps_edge_ex_middle),
                Box::new(ss_edge_ex_middle.rev()),
            ],
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                3,
                EdgeIndex::IMin,
                0..2,
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMin, 1..2),
                    BlockBoundaryRange::new(&mesh, 3, EdgeIndex::JMin, 0..1),
                ),
            )));
    }

    {
        let block_name = "exit_ss";

        let edge_j_min = mesh.blocks.last().unwrap().edge_data(EdgeIndex::JMax);

        let x_01 = *edge_j_min.x.last().unwrap();
        let x_10 = ss_edge_ex_ss.point_coord(ss_edge_ex_ss.start);
        let x_11 = trailing_edge + Vec2d(0.0, 0.5 * pitch);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![Box::new(ss_edge_ex_ss.rev())],
            vec![Box::new(Segment::new(
                ss_edge_ex_ss.len(),
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                4,
                EdgeIndex::IMin,
                0..1,
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 3, EdgeIndex::JMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMin, 0..1),
                ),
            )));
    }

    {
        let block_name = "inlet_ss";

        let edge_j_min = mesh.blocks.last().unwrap().edge_data(EdgeIndex::JMax);
        let edge_i_min_last_segment = mesh.blocks[0].edge_data(EdgeIndex::JMin);

        let x_11 = leading_edge + Vec2d(0.0, 0.5 * pitch);
        let x_10 = *edge_i_min_last_segment.x.last().unwrap();
        let x_01 = *edge_j_min.x.last().unwrap();

        let num_points_i_max = edge_i_min_last_segment.len() + ss_edge_in_ss.len() - 1;

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![
                Box::new(ss_edge_in_ss.rev()),
                Box::new(edge_i_min_last_segment),
            ],
            vec![Box::new(Segment::new(
                num_points_i_max,
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

        mesh.boundaries
            .push(BlockBoundary::Wall(BlockBoundaryRange::new(
                &mesh,
                5,
                EdgeIndex::IMin,
                0..1,
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 4, EdgeIndex::JMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMin, 0..1),
                ),
            )));
        mesh.boundaries.push(BlockBoundary::PeriodicConnection {
            connection: (
                BlockBoundaryRange::new(&mesh, 2, EdgeIndex::IMax, 0..1),
                BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMax, 0..1),
            ),
            translation: pitch,
        });
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 0, EdgeIndex::JMin, 0..1),
                    BlockBoundaryRange::new(&mesh, 5, EdgeIndex::IMin, 1..2),
                ),
            )));
    }

    {
        let block_name = "inlet";

        // copy from inlet_lower
        let edge_j_max_seg_0 = mesh.blocks[1].edge_data(EdgeIndex::IMax);
        // copy from inlet_middle
        let edge_j_max_seg_1 = mesh.blocks[0].edge_data(EdgeIndex::IMax);
        // copy from inlet_ss
        let edge_j_max_seg_2 = mesh.blocks[5].edge_data(EdgeIndex::JMax);

        let num_points_j =
            edge_j_max_seg_0.len() + edge_j_max_seg_1.len() + edge_j_max_seg_2.len() - 2;

        let x_10 = *edge_j_max_seg_0.x.last().unwrap();
        let x_11 = *edge_j_max_seg_2.x.last().unwrap();
        let x_00 = x_10 + Vec2d(-0.05, 0.0);
        let x_01 = x_00 + Vec2d(0.0, pitch);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![Box::new(Segment::new(
                num_cells_inlet + 1,
                UniformClustering::new(),
                Line2d::new(x_00, x_10),
            ))],
            vec![Box::new(Segment::new(
                num_cells_inlet + 1,
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

        mesh.boundaries
            .push(BlockBoundary::Inlet(BlockBoundaryRange::new(
                &mesh,
                6,
                EdgeIndex::JMin,
                0..1,
            )));
        mesh.boundaries.push(BlockBoundary::PeriodicConnection {
            connection: (
                BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMin, 0..1),
                BlockBoundaryRange::new(&mesh, 6, EdgeIndex::IMax, 0..1),
            ),
            translation: pitch,
        });
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 1, EdgeIndex::IMax, 0..1).reverse(),
                    BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 0..1),
                ),
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 0, EdgeIndex::IMax, 0..1).reverse(),
                    BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 1..2),
                ),
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 5, EdgeIndex::JMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 6, EdgeIndex::JMax, 2..3),
                ),
            )));
    }

    {
        let block_name = "exit";

        // copy from exit_lower
        let edge_j_min_seg_0 = mesh.blocks[2].edge_data(EdgeIndex::JMax);
        // copy from exit_middle
        let edge_j_min_seg_1 = mesh.blocks[3].edge_data(EdgeIndex::IMax);
        // copy from exit_ss
        let edge_j_min_seg_2 = mesh.blocks[4].edge_data(EdgeIndex::IMax);

        let num_points_j =
            edge_j_min_seg_0.len() + edge_j_min_seg_1.len() + edge_j_min_seg_2.len() - 2;

        let x_00 = *edge_j_min_seg_0.x.last().unwrap();
        let x_01 = *edge_j_min_seg_2.x.last().unwrap();
        let x_10 = x_00 + Vec2d(0.05, 0.0);
        let x_11 = x_10 + Vec2d(0.0, pitch);

        let block = Block2d::new(
            row_prefix.to_owned() + block_name,
            vec![Box::new(Segment::new(
                num_cells_outlet + 1,
                UniformClustering::new(),
                Line2d::new(x_00, x_10),
            ))],
            vec![Box::new(Segment::new(
                num_cells_outlet + 1,
                UniformClustering::new(),
                Line2d::new(x_01, x_11),
            ))],
            vec![
                Box::new(edge_j_min_seg_0.reverse()),
                Box::new(edge_j_min_seg_1),
                Box::new(edge_j_min_seg_2),
            ],
            vec![Box::new(Segment::new(
                num_points_j,
                UniformClustering::new(),
                Line2d::new(x_10, x_11),
            ))],
        );

        mesh.add_block(block);

        mesh.boundaries
            .push(BlockBoundary::Outlet(BlockBoundaryRange::new(
                &mesh,
                7,
                EdgeIndex::JMax,
                0..1,
            )));
        mesh.boundaries.push(BlockBoundary::PeriodicConnection {
            connection: (
                BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMin, 0..1),
                BlockBoundaryRange::new(&mesh, 7, EdgeIndex::IMax, 0..1),
            ),
            translation: pitch,
        });
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 2, EdgeIndex::JMax, 0..1).reverse(),
                    BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 0..1),
                ),
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 3, EdgeIndex::IMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 1..2),
                ),
            )));
        mesh.boundaries
            .push(BlockBoundary::Connection(BlockConnection::new(
                &mesh,
                (
                    BlockBoundaryRange::new(&mesh, 4, EdgeIndex::IMax, 0..1),
                    BlockBoundaryRange::new(&mesh, 7, EdgeIndex::JMin, 2..3),
                ),
            )));
    }

    // mesh.blocks
    //     .iter()
    //     .for_each(|block| println!("{:#?}", block.coords.dim()));

    // println!("{:#?}", mesh.boundaries);

    mesh.smooth();

    // smooth_mesh(&mut mesh);

    // compute_derivatives(&mesh);

    // mesh.blocks.iter_mut().for_each(|block| {
    //     println!("block: {}", block.name);
    //     smooth_block(block, 100).unwrap()
    // });

    // plot blocking
    let ps_spline_int = ps_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
    let ss_spline_int = ss_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());

    let ps_clustering_int = &ps_edge.edge().coords;
    let ss_clustering_int = &ss_edge.edge().coords;

    let root = BitMapBackend::new("blocking.png", (1900, 1200)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        // .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.98f32..1.2f32, -0.1f32..0.1f32)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            ps_spline_int.iter().map(|v| (v[0] as f32, v[1] as f32)),
            &RED,
        ))
        .unwrap()
        .label("PS");
    chart
        .draw_series(LineSeries::new(
            ss_spline_int.iter().map(|v| (v[0] as f32, v[1] as f32)),
            &BLUE,
        ))
        .unwrap()
        .label("SS");

    chart
        .draw_series(
            ps_clustering_int
                .iter()
                .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLACK.filled())),
        )
        .unwrap();
    chart
        .draw_series(
            ss_clustering_int
                .iter()
                .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLACK.filled())),
        )
        .unwrap();

    for block in mesh.blocks.iter() {
        chart
            .draw_series(
                block
                    .coords
                    .iter()
                    .map(|v| Circle::new((v.0 as f32, v.1 as f32), 3, BLUE.filled())),
            )
            .unwrap();
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();

    // // mesh creation
    // let mut mesh = Mesh::new();

    // mesh.add_block(Block2d::new(
    //     String::from(row_prefix.to_owned() + "in_m"),
    //     (ninm + 1, ncut),
    // ));

    // mesh.save("turbomesh.cgns").unwrap();

    (Geometry {}, mesh)
}
