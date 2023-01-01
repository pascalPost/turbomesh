// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

#[cfg(test)]
use crate::clustering::roberts_clustering;

#[cfg(test)]
use crate::{Block2d, Mesh, Vec2d};

#[cfg(test)]
use ndarray::Array;

#[cfg(test)]
use plotters::prelude::*;

use crate::interpolation::FittingSpline;
use std::error::Error;

/// read the profile coordinates from given file and returns them as a vec
fn read_profile(file_name: &str) -> Result<Vec<[f64; 2]>, Box<dyn Error>> {
    let mut v: Vec<[f64; 2]> = vec![];

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

        let x: f64 = record[0].parse()?;
        let y: f64 = record[1].parse()?;

        v.push([x, y]);
    }

    assert!(v.len() > 0, "parsed profile data is empty.");

    Ok(v)
}

pub fn blade_profile() -> Result<(FittingSpline<2>, FittingSpline<2>), Box<dyn Error>> {
    // read suction and pressure side coordinates from csv and save into vec

    let ps = read_profile("examples/T106/T106_ps.dat")
        .expect("error in reading pressure side coordinates");

    let ss = read_profile("examples/T106/T106_ss.dat")
        .expect("Error in reading suction side coordinates");

    // println!("ps[0]: {:?}", ps[0]);
    // println!("ps[-1]: {:?}", ps[ps.len() - 1]);
    // println!("ss[0]: {:?}", ss[0]);
    // println!("ss[-1]: {:?}", ss[ss.len() - 1]);

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

    // // plot data for debugging
    // let ps_spline_int = ps_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
    // let ss_spline_int = ss_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());

    // let mut profile = ps.clone();
    // profile.extend_from_slice(&ss[1..ss.len() - 2]);
    // let profile = profile;

    // let profile_spline = FittingSpline::new(profile.as_slice(), None, None, 3)
    //     .expect("Error in spline interpolation");
    // let profile_int = ss_spline.interpolate(Array::linspace(0., 1., 2000).as_slice().unwrap());

    // let root = BitMapBackend::new("profile.png", (1900, 1200)).into_drawing_area();
    // root.fill(&WHITE).unwrap();
    // let mut chart = ChartBuilder::on(&root)
    //     // .caption("y=x^2", ("sans-serif", 50).into_font())
    //     .margin(5)
    //     .x_label_area_size(30)
    //     .y_label_area_size(30)
    //     .build_cartesian_2d(1.04f32..1.13f32, -0.05f32..0.02f32)
    //     .unwrap();

    // chart.configure_mesh().draw().unwrap();

    // chart
    //     .draw_series(LineSeries::new(
    //         profile_int.iter().map(|v| (v[0] as f32, v[1] as f32)),
    //         &BLACK,
    //     ))
    //     .unwrap()
    //     .label("profile");
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
    //         &ORANGE_A100,
    //     ))
    //     .unwrap()
    //     .label("SS");

    // chart
    //     .configure_series_labels()
    //     .background_style(&WHITE.mix(0.8))
    //     .border_style(&BLACK)
    //     .draw()
    //     .unwrap();

    // root.present().unwrap();

    Ok((ps_spline, ss_spline))
}

#[cfg(test)]
#[test]
fn blocking() {
    use plotters::style::full_palette::ORANGE;

    let num_cells_blade: usize = 60;
    let num_cells_in_middle_half: usize = 5;
    let scut: usize = 20;
    let num_cells_cut: usize = 10;

    let (ps_spline, ss_spline) = blade_profile().expect("Error in profile input and interpolation");

    // clustering around the blade
    let blade_clustering = roberts_clustering(num_cells_blade, 0.5, 1.01);

    // TODO add approximation of the leading and trailing edge by approximating
    // the chamber line and computing its intersection with the profile.
    // For now, the values of the given profiles are used.

    let mut mesh = Mesh::new();

    let row_prefix = "row_01_";

    // in_middle block

    let block_in_middle = Block2d::new(
        String::from(row_prefix.to_owned() + "block_in_middle"),
        (num_cells_in_middle_half * 2, num_cells_cut),
        vec![],
    )
    .corners(
        Vec2d(0.0, 0.0),
        Vec2d(1.0, 0.0),
        Vec2d(1.0, 1.0),
        Vec2d(0.0, 1.0),
    );

    mesh.add_block(block_in_middle);

    // let x_in_middle_lower_start = x_blade[-(num_cells_in_middle / 2) - 1];

    // plt.plot(x_inm_lower_start, y_inm_lower_start, '.r')

    // x_inm_upper_start = x_blade[ninm // 2]
    // y_inm_upper_start = y_blade[ninm // 2]

    // plt.plot(x_inm_upper_start, y_inm_upper_start, '.r')

    // x_inm_lower_end = x_blade[-(ninm // 2) - 1] - 0.01
    // y_inm_lower_end = y_blade[-(ninm // 2) - 1] - 0.01

    // x_inm_upper_end = x_blade[ninm // 2] - 0.01
    // y_inm_upper_end = y_blade[ninm // 2] + 0.01

    // plt.plot([x_inm_lower_start, x_inm_lower_end], [
    //          y_inm_lower_start, y_inm_lower_end], 'r')
    // plt.plot([x_inm_upper_start, x_inm_upper_end], [
    //          y_inm_upper_start, y_inm_upper_end], 'r')
    // plt.plot([x_inm_lower_end, x_inm_upper_end], [
    //          y_inm_lower_end, y_inm_upper_end], 'r')

    // plot blocking
    let ps_spline_int = ps_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());
    let ss_spline_int = ss_spline.interpolate(Array::linspace(0., 1., 1000).as_slice().unwrap());

    let ps_clustering_int = ps_spline.interpolate(blade_clustering.as_slice().unwrap());
    let ss_clustering_int = ss_spline.interpolate(blade_clustering.as_slice().unwrap());

    let root = BitMapBackend::new("blocking.png", (1900, 1200)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        // .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(1.04f32..1.13f32, -0.05f32..0.02f32)
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
            &ORANGE,
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
}
