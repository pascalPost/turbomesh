// Copyright (c) Pascal Post. All Rights Reserved.
// Licensed under AGPLv3 license (see LICENSE.txt for details)

use crate::geometry::Spline;
use crate::interpolation::FittingSpline;
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
