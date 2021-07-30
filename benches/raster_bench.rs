extern crate rust_gfx;

use criterion::{criterion_group, criterion_main, black_box, Criterion};

use rust_gfx::camera::default_camera;
use rust_gfx::matrix::translation;
use rust_gfx::rasterize::*;
use rust_gfx::shape::colour_cube;
use rust_gfx::triangle::{Triangle, TriangleComputations};
use rust_gfx::vector::Vector3;

fn triangle_benchmark(c: &mut Criterion) {
    let t = Triangle::new([
        Vector3::new(150.0, 450.0, 0.0),
        Vector3::new(150.0, 150.0, 0.0),
        Vector3::new(450.0, 450.0, 0.0),
    ]);

    let comps = TriangleComputations::new(&t, 314.0 + 0.5, 372.0 + 0.5);

    c.bench_function("triangle_computations", |b| {
        b.iter(|| black_box(TriangleComputations::new(&t, 314.0 + 0.5, 372.0 + 0.5)))
    });

    c.bench_function("get_overlap", |b| b.iter(|| black_box(comps.get_overlap())));
}

fn raster_benchmark(c: &mut Criterion) {
    let width = 100;
    let height = 100;

    let canvas = Canvas::new(width, height);
    let camera = default_camera(width as f32, height as f32);

    let mut z_buffer = vec![0.0f32; (width * height) as usize];
    let mut frame = vec![0u8; (4 * width * height) as usize];

    let mut shapes = vec![];
    shapes.push(colour_cube());

    let translation = translation(0.0, 0.0, -5.0);
    shapes.iter_mut().for_each(|s| s.transform = translation);

    c.bench_function("to_raster_space", |b| {
        b.iter(|| black_box(to_raster_space(-0.5, -0.5, width, height)))
    });

    c.bench_function("set_pixel", |b| {
        b.iter(|| canvas.set_pixel(50, 50, RED, frame.as_mut_slice()))
    });
    
    c.bench_function("draw_cube", |b| {
        b.iter(|| {
            canvas.draw(
                &camera,
                &shapes,
                frame.as_mut_slice(),
                z_buffer.as_mut_slice(),
            )
        })
    });
}

criterion_group!(benches, triangle_benchmark, raster_benchmark);
criterion_main!(benches);
