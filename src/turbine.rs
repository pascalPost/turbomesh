// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::clustering::{RobertsClustering, UniformClustering};
use crate::geometry::{Line2d, Spline};
use crate::interpolation::FittingSpline;
use crate::types::{ClusteringFunction, Edge, EdgeView, MappingFunction};
use crate::{Block2d, Geometry, Mesh, Scalar, Segment, Vec2d};
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
    let num_cells_blade: usize = 60;
    let num_cells_in_middle_half: usize = 5;
    let num_cells_scut: usize = 20;
    let num_cells_cut: usize = 10;

    let (ps_spline, ss_spline) =
        blade_profile(ps_csv_path, ss_csv_path).expect("Error in profile input and interpolation");

    // blade discretization on which the 2d blocking is based
    let blade_clustering_function = RobertsClustering::new(0.5, 1.01);

    let ps_edge = EdgeView::new(Edge::new(
        "Pressure_Side_Edge".to_string(),
        &vec![Segment::new(
            num_cells_blade + 1,
            blade_clustering_function,
            ps_spline.clone(),
        )],
    ));

    let ss_edge = EdgeView::new(Edge::new(
        "Suction_Side_Edge".to_string(),
        &vec![Segment::new(
            num_cells_blade + 1,
            blade_clustering_function,
            ss_spline.clone(),
        )],
    ));

    // split blade distribution

    let (ps_edge_in_middle, ps_edge_rest) = ps_edge.split_at(num_cells_in_middle_half);
    let (ss_edge_in_middle, ss_edge_rest) = ss_edge.split_at(num_cells_in_middle_half);

    // TODO add approximation of the leading and trailing edge by approximating
    // the chamber line and computing its intersection with the profile.
    // For now, the values of the given profiles are used.
    // let leading_edge = ps_spline.interpolate_val(0.0);

    let mut mesh = Mesh::new();

    let row_prefix = "row_01_";

    // in_middle block

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

        let x_00 = ss_edge_in_middle.point_coord(ss_edge_in_middle.end);
        let x_10 = ps_edge_in_middle.point_coord(ps_edge_in_middle.end);

        // TODO remove hard coded position for block
        let x_01 = x_00 - Vec2d(0.007, -0.025);

        // TODO remove hard coded position for block
        let x_11 = x_10 - Vec2d(0.007, 0.001);

        let block_in_middle = Block2d::new(
            row_prefix.to_owned() + "in_middle",
            // (num_cells_in_middle_half * 2, num_cells_cut),
            vec![
                Box::new(ss_edge_in_middle.rev()),
                Box::new(ps_edge_in_middle),
            ],
            vec![Box::new(Segment::new(
                num_cells_in_middle_half * 2 + 1,
                UniformClustering::new(),
                Line2d::new(x_01, x_11),
            ))],
            vec![Box::new(Segment::new(
                num_cells_cut + 1,
                UniformClustering::new(),
                Line2d::new(x_00, x_01),
            ))],
            vec![Box::new(Segment::new(
                num_cells_cut + 1,
                UniformClustering::new(),
                Line2d::new(x_10, x_11),
            ))],
        );

        mesh.add_block(block_in_middle);
    }

    // // in_lower block

    // {
    //     // on blade start
    //     let x_00 = mesh.blocks[0].edge_j_max[0]
    //         .mapping
    //         .computational_to_physical_val(0.0);

    //     let x_01 = mesh.blocks[0].edge_j_max[0]
    //         .mapping
    //         .computational_to_physical_val(1.0);

    //     // on blade end
    //     let x_10 = ps_spline
    //         .interpolate(&[blade_clustering[num_cells_in_middle_half + num_cells_scut]])[0];
    //     let x_10 = Vec2d(x_10[0], x_10[1]);

    //     // on periodic bc
    //     let x_leading_edge = ss_spline.interpolate(&[0.0]);
    //     let x_leading_edge = Vec2d(x_leading_edge[0][0], x_leading_edge[0][1]);
    //     let x_11 = x_leading_edge + Vec2d(0.0, -0.5 * pitch);

    //     let i_min_clustering = {
    //         let blade_clustering = blade_clustering.clone();
    //         move |u: &mut [Scalar]| {
    //             u[0] = 0.0;

    //             let idx_start = num_cells_in_middle_half + 1;
    //             let idx_end = idx_start + num_cells_scut + 1;

    //             blade_clustering[idx_start..idx_end]
    //                 .windows(2)
    //                 .zip(u[1..].iter_mut())
    //                 .for_each(|(x, u)| *u = x[1] - x[0]);

    //             // cum sum
    //             let u_max = u.iter_mut().fold(0.0, |acc, x| {
    //                 *x += acc;
    //                 *x
    //             });

    //             // renormalization to 0<=u<=1
    //             u.iter_mut().for_each(|u| *u /= u_max);
    //         }
    //     };

    //     let idx_start = num_cells_in_middle_half;
    //     let idx_end = idx_start + num_cells_scut;

    //     let v_start = blade_clustering[idx_start];
    //     let v_end = blade_clustering[idx_end];

    //     let i_min_curve = SplineSegment::new(&ps_spline, v_start, v_end);

    //     println!("{}", x_00);
    //     println!(
    //         "{:?}",
    //         ps_spline.interpolate(&[blade_clustering[num_cells_in_middle_half]])[0]
    //     );
    //     println!("{}", i_min_curve.computational_to_physical_val(0.0));

    //     let block_in_lower = Block2d::new(
    //         row_prefix.to_owned() + "in_lower",
    //         (num_cells_scut, num_cells_cut),
    //         vec![EdgeSegment::new(i_min_clustering, i_min_curve).into()],
    //         vec![EdgeSegment::new(uniform_clustering, Line2d::new(x_01, x_11)).into()],
    //         vec![EdgeSegment::new(uniform_clustering, Line2d::new(x_00, x_01)).into()],
    //         vec![EdgeSegment::new(uniform_clustering, Line2d::new(x_10, x_11)).into()],
    //     );

    //     mesh.add_block(block_in_lower);
    // }

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
        .build_cartesian_2d(1.00f32..1.13f32, -0.05f32..0.02f32)
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
                    .as_slice()
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
