// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

fn main() {
    let (_, _) = turbomesh::turbine::run_turbine_template(
        "examples/T106/T106_ps.dat",
        "examples/T106/T106_ss.dat",
    );
}
