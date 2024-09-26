#[macro_use]
extern crate approx;

pub mod camera;
pub mod matrix;
pub mod rasterize;
pub mod shape;
pub mod triangle;
pub mod vector;

// Why not f32::EPSILON? This is defined as 1.1920929e-07, which is too small.
// In graphics applications, due to larger scales and rounding errors from multiple operations (like transformations, shading, etc.), 
// we want a larger tolerance threshold. A value like 0.0001 or 0.001 is often used.
pub const EPSILON: f32 = 0.001;
