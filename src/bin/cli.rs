// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use clap::Parser;
use figment::{
    providers::{Format, Yaml},
    Figment,
};
use turbomesh::turbine::TurbineTemplate;

#[derive(Parser)]
#[command(author = "Pascal Post", about = "turbomesh CLI", long_about = None)]
struct Cli {
    /// The path to the config file to read
    path: std::path::PathBuf,
}

fn main() {
    env_logger::init();
    let args = Cli::parse();

    let config_file_path = args.path;

    let mut turbine_template: TurbineTemplate = Figment::new()
        .merge(Yaml::file(config_file_path.to_str().unwrap()))
        .extract()
        .unwrap();

    // update file pathes to be relative to config file
    turbine_template.append_root_path(config_file_path.parent().unwrap());

    println!("{:?}", turbine_template);

    turbine_template.run();
}
