// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use clap::Parser;
use figment::{
    providers::{Format, Yaml},
    Figment,
};
use log::debug;
use serde::Deserialize;
use turbomesh::turbine::TurbineTemplate;

/// struct representing the CLI arguments
#[derive(Parser)]
#[command(author = "Pascal Post", about = "turbomesh CLI", long_about = None)]
struct Cli {
    /// The path to the config file to read
    path: std::path::PathBuf,
}

/// enum listing all available meshing templates
#[derive(Deserialize, Debug)]
#[serde(tag = "template")]
enum Template {
    Turbine(TurbineTemplate),
}

fn main() {
    env_logger::init();
    let args = Cli::parse();

    let config_file_path = args.path;

    let template: Template = Figment::new()
        .merge(Yaml::file(config_file_path.to_str().unwrap()))
        .extract()
        .unwrap();

    debug!("cli arguments: {:?}", template);

    match template {
        Template::Turbine(mut turbine_template) => {
            // update file pathes to be relative to config file
            turbine_template.append_root_path(config_file_path.parent().unwrap());

            turbine_template.run();
        }
    }
}
