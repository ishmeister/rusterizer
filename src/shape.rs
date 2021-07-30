use super::matrix::{identity, Matrix4x4};
use super::rasterize::*;
use super::triangle::Triangle;
use super::vector::Vector3;

pub struct Shape {
    pub transform: Matrix4x4<f32>,
    pub triangles: Vec<Triangle<f32>>,
}

pub fn colour_cube() -> Shape {
    let vertices = [
        [-1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0]
    ];

    let triangle_idxs = [
        [0, 1, 2], [2, 3, 0], 
        [7, 6, 5], [5, 4, 7],
        [4, 5, 1], [1 ,0, 4],
        [3, 2, 6], [6 ,7, 3],
        [4, 0, 3], [3, 7, 4],
        [1, 5, 6], [6, 2, 1]
    ];

    let colours = [
        BLUE, RED, GREEN, GREEN, RED, BLUE, 
        BLUE, RED, GREEN, GREEN, RED, BLUE, 
        RED, GREEN, RED, RED, BLUE, RED, 
        RED, GREEN, RED, RED, BLUE, RED, 
        RED, BLUE, RED, RED, BLUE, RED, 
        RED, GREEN, RED, RED, GREEN, RED,
    ];

    let mut triangles = vec![];

    for (i, tri_idx) in triangle_idxs.iter().enumerate() {
        let mut triangle = Triangle::new([
            Vector3::from_array(vertices[tri_idx[0]]),
            Vector3::from_array(vertices[tri_idx[1]]),
            Vector3::from_array(vertices[tri_idx[2]]),
        ]);

        let col_idx = i * 3;
        
        triangle.colours = [
            colours[col_idx], 
            colours[col_idx + 1], 
            colours[col_idx + 2]];

        triangle.double_side = false;
        triangles.push(triangle);
    }

    Shape {
        transform: identity(),
        triangles,
    }
}
