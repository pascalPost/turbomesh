use std::error::Error;
use turbomesh::types::boundary::{BlockBoundary, BlockBoundaryRange, PeriodicBlockConnection};
use turbomesh::types::EdgeIndex;
use turbomesh::{
    clustering::{RobertsClustering, UniformClustering},
    geometry::Rectangle,
    smoothing::smooth_mesh,
    Block2d, Mesh, Segment, Vec2d,
};

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();

    let mut mesh = Mesh::new();

    let pitch = 1.0;

    let rect = Rectangle::new(Vec2d(0.0, 0.0), Vec2d(pitch, 1.0));

    let num_points = (25, 20);

    let alpha = 0.5;
    let beta = 1.03;

    let block = Block2d::new(
        String::from("block"),
        vec![Box::new(Segment::new(
            num_points.0,
            UniformClustering::new(),
            rect.i_min(),
        ))],
        vec![Box::new(Segment::new(
            num_points.0,
            UniformClustering::new(),
            rect.i_max(),
        ))],
        vec![Box::new(Segment::new(
            num_points.1,
            RobertsClustering::new(alpha, beta),
            rect.j_min(),
        ))],
        vec![Box::new(Segment::new(
            num_points.1,
            RobertsClustering::new(alpha, beta),
            rect.j_max(),
        ))],
    );

    mesh.add_block(block);

    mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        &mesh,
        0,
        EdgeIndex::IMin,
        0..1,
    )));

    mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        &mesh,
        0,
        EdgeIndex::IMax,
        0..1,
    )));

    mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        &mesh,
        0,
        EdgeIndex::JMin,
        0..1,
    )));

    mesh.edges.push(BlockBoundary::Wall(BlockBoundaryRange::new(
        &mesh,
        0,
        EdgeIndex::JMax,
        0..1,
    )));

    mesh.save("test.cgns")?;

    smooth_mesh(&mut mesh, 1)?;

    mesh.save("test_smoothed.cgns")?;

    Ok(())
}
