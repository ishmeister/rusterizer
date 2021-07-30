#[macro_use]
extern crate approx;

pub mod camera;
pub mod matrix;
pub mod rasterize;
pub mod shape;
pub mod triangle;
pub mod vector;

pub const EPSILON: f32 = 0.0001;
