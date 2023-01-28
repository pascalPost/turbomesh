// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::clustering::{roberts_clustering, uniform_clustering};
use crate::geometry::Line2d;
use crate::interpolation::FittingSpline;
use crate::{Block2d, Edge, Geometry, Mesh, Scalar, Vec2d};
use ndarray::Array;
use plotters::prelude::*;
use std::error::Error;
use std::rc::Rc;

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
) -> Result<(FittingSpline<2>, FittingSpline<2>), Box<dyn Error>> {
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

    Ok((ps_spline, ss_spline))
}

pub fn run_turbine_template(ps_csv_path: &str, ss_csv_path: &str) -> (Geometry, Mesh) {
    let mut geometry = Geometry {};
    // use plotters::style::full_palette::ORANGE;

    let num_cells_blade: usize = 60;
    let num_cells_in_middle_half: usize = 5;
    let scut: usize = 20;
    let num_cells_cut: usize = 10;

    let (ps_spline, ss_spline) =
        blade_profile(ps_csv_path, ss_csv_path).expect("Error in profile input and interpolation");

    // clustering around the blade
    let blade_clustering = Rc::<[Scalar]>::from(roberts_clustering(num_cells_blade, 0.5, 1.01));

    // TODO add approximation of the leading and trailing edge by approximating
    // the chamber line and computing its intersection with the profile.
    // For now, the values of the given profiles are used.

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

        let u_on_blade = blade_clustering[num_cells_in_middle_half];

        // TODO insert function to convert [Scalar; 2] to Vec2d

        let x_00 = ss_spline.interpolate(&[u_on_blade])[0];
        let x_00 = Vec2d(x_00[0], x_00[1]);

        // TODO remove hard coded position for block
        let x_01 = x_00 - Vec2d(0.007, -0.025);

        let x_10 = ps_spline.interpolate(&[u_on_blade])[0];
        let x_10 = Vec2d(x_10[0], x_10[1]);

        // TODO remove hard coded position for block
        let x_11 = x_10 - Vec2d(0.007, 0.001);

        println!("{x_01} {x_11}");

        let i_min_clustering = {
            let blade_clustering = blade_clustering.clone();
            move |u: &mut [Scalar]| {
                // add 0 at beginning for cum sum at the end
                u[0] = 0.0;

                // set the initial clustering for the first half as the
                // difference of the blade clustering
                blade_clustering[0..num_cells_in_middle_half + 1]
                    .windows(2)
                    .rev()
                    .zip(u[1..].iter_mut())
                    .for_each(|(x, u)| *u = x[1] - x[0]);

                // set the second half from the already computed half as the
                // clustering is symmetric around the leading edge
                let (u_lower, u_upper) = u.split_at_mut(num_cells_in_middle_half + 1);
                u_upper
                    .iter_mut()
                    .zip(u_lower.iter().rev())
                    .for_each(|(u, val)| *u = *val);

                // cum sum
                let u_max = u.iter_mut().fold(0.0, |acc, x| {
                    *x += acc;
                    *x
                });

                // renormalization to 0<=u<=1
                u.iter_mut().for_each(|u| *u /= u_max);
            }
        };

        let block_in_middle = Block2d::new(
            row_prefix.to_owned() + "in_middle",
            (num_cells_in_middle_half * 2, num_cells_cut),
            Edge::new(i_min_clustering, Line2d::new(x_00, x_10)),
            Edge::new(uniform_clustering, Line2d::new(x_01, x_11)),
            Edge::new(uniform_clustering, Line2d::new(x_00, x_01)),
            Edge::new(uniform_clustering, Line2d::new(x_10, x_11)),
        );

        mesh.add_block(block_in_middle);
    }

    // plot blocking
    let ps_spline_int = ps_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
    let ss_spline_int = ss_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());

    let ps_clustering_int = ps_spline.interpolate(blade_clustering.as_ref());
    let ss_clustering_int = ss_spline.interpolate(blade_clustering.as_ref());

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
                .map(|v| Circle::new((v[0] as f32, v[1] as f32), 3, BLACK.filled())),
        )
        .unwrap();
    chart
        .draw_series(
            ss_clustering_int
                .iter()
                .map(|v| Circle::new((v[0] as f32, v[1] as f32), 3, BLACK.filled())),
        )
        .unwrap();

    for block in mesh.blocks.iter() {
        println!("{:#?}", block.coords);
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

    (geometry, mesh)
}
