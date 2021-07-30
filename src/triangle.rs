use crate::rasterize::Colour;
use crate::vector::Vector3;

#[derive(Debug, Clone)]
pub struct Triangle<T> {
    pub points: [Vector3<T>; 3],
    pub colours: [Colour; 3],
    pub normal: Vector3<f32>,
    pub double_side: bool,
}

impl Triangle<f32> {
    pub fn new(points: [Vector3<f32>; 3]) -> Self {
        let colours = [crate::rasterize::GREY; 3];
        let edge0 = points[1] - points[0];
        let edge1 = points[2] - points[0];
        let normal = edge0.cross(edge1).normalize();

        Triangle::<f32> {
            points,
            colours,
            normal,
            double_side: false,
        }
    }
}

#[derive(Debug)]
pub struct TriangleComputations<T> {
    pub edges: [Vector3<T>; 3],
    z_inv: [f32; 3],
    w0_params: [f32; 5],
    w1_params: [f32; 5],
    w2_params: [f32; 5],
}

impl TriangleComputations<f32> {
    pub fn new(t: &Triangle<f32>, x_start: f32, y_start: f32) -> Self {
        let edges = [
            t.points[1] - t.points[2],
            t.points[2] - t.points[0],
            t.points[0] - t.points[1],
        ];

        // area of parallelogram (2x triangle area)
        let area = edge_function(t.points[0], t.points[1], t.points[2]);

        // pre-compute linear edge function variables
        let mut w0_params = edge_function_params(t.points[1], t.points[2], x_start, y_start);
        let mut w1_params = edge_function_params(t.points[2], t.points[0], x_start, y_start);
        let mut w2_params = edge_function_params(t.points[0], t.points[1], x_start, y_start);

        // scale each param by area to avoid division in the main loop
        for i in 0..5 {
            w0_params[i] /= area;
            w1_params[i] /= area;
            w2_params[i] /= area;
        }

        // pre-computed reciprocals
        let z_inv = [
            1.0 / t.points[0].z,
            1.0 / t.points[1].z,
            1.0 / t.points[2].z,
        ];

        TriangleComputations::<f32> {
            edges: edges,
            z_inv: z_inv,
            w0_params,
            w1_params,
            w2_params,
        }
    }

    pub fn step_x(&mut self) {
        self.w0_params[3] += self.w0_params[2];
        self.w1_params[3] += self.w1_params[2];
        self.w2_params[3] += self.w2_params[2];
    }

    pub fn step_y(&mut self) {
        self.w0_params[1] += self.w0_params[0];
        self.w0_params[3] = self.w0_params[4];
        self.w1_params[1] += self.w1_params[0];
        self.w1_params[3] = self.w1_params[4];
        self.w2_params[1] += self.w2_params[0];
        self.w2_params[3] = self.w2_params[4];
    }

    pub fn get_overlap(&self) -> (bool, f32, f32, f32, f32) {
        let (w0, w1, w2) = self.get_barycentric_coords();

        let overlaps = (w0 > 0.0 || (w0 == 0.0 && is_top_or_left_edge(self.edges[0])))
            && (w1 > 0.0 || (w1 == 0.0 && is_top_or_left_edge(self.edges[1])))
            && (w2 > 0.0 || (w2 == 0.0 && is_top_or_left_edge(self.edges[2])));

        let mut point_z = f32::MAX;
        if overlaps {
            // find the correct z of point for the depth buffer by linear interpolation
            // the reciprocal of z should already be pre-computed
            point_z = 1.0 / (self.z_inv[0] * w0 + self.z_inv[1] * w1 + self.z_inv[2] * w2);
        }

        (overlaps, w0, w1, w2, point_z)
    }

    // These barycentric coordinates are in the range [0,1] with sum equal to 1
    // Any point on the triangle is a weighted average of the the triangle vertices and the barycentric coords
    // Use to determine the relative weighting of triangle attributes from each vertex
    pub fn get_barycentric_coords(&self) -> (f32, f32, f32) {
        let w0 = self.w0_params[1] + self.w0_params[3];
        let w1 = self.w1_params[1] + self.w1_params[3];
        let w2 = self.w2_params[1] + self.w2_params[3];
        (w0, w1, w2)
    }
}

fn edge_function_params(
    p1: Vector3<f32>,
    p2: Vector3<f32>,
    x_start: f32,
    y_start: f32,
) -> [f32; 5] {
    let edge = p1 - p2;
    let mut w0_params = [0.0f32; 5];
    w0_params[0] = edge.x;
    w0_params[1] = w0_params[0] * (y_start - p1.y);
    w0_params[2] = -edge.y;
    w0_params[3] = w0_params[2] * (x_start - p1.x);
    w0_params[4] = w0_params[3];
    w0_params
}

// top edges are flat and left edges have a rising y
fn is_top_or_left_edge(edge: Vector3<f32>) -> bool {
    edge.y == 0.0 && edge.x > 0.0 || edge.y > 0.0
}

fn edge_function(a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>) -> f32 {
    (a.x - b.x) * (c.y - a.y) - (a.y - b.y) * (c.x - a.x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EPSILON;

    fn test_triangle1() -> Triangle<f32> {
        Triangle::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ])
    }

    fn test_triangle2() -> Triangle<f32> {
        Triangle::new([
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ])
    }

    #[test]
    fn get_overlap_inside() {
        let t1 = TriangleComputations::new(&test_triangle1(), 0.4, 0.4);
        let t2 = TriangleComputations::new(&test_triangle2(), 0.8, 0.8);
        assert_eq!(true, t1.get_overlap().0);
        assert_eq!(true, t2.get_overlap().0);
    }

    #[test]
    fn get_overlap_on_top_edge() {
        let t = TriangleComputations::new(&test_triangle2(), 0.5, 1.0);
        assert_eq!(true, t.get_overlap().0);
    }

    #[test]
    fn get_overlap_on_left_edge() {
        let t1 = TriangleComputations::new(&test_triangle1(), 0.0, 0.5);
        let t2 = TriangleComputations::new(&test_triangle2(), 0.5, 0.5);
        assert_eq!(true, t1.get_overlap().0);
        assert_eq!(true, t2.get_overlap().0);
    }

    #[test]
    fn does_not_contain_point_on_right_edge() {
        let t1 = TriangleComputations::new(&test_triangle1(), 0.5, 0.5);
        let t2 = TriangleComputations::new(&test_triangle2(), 1.0, 0.5);
        assert_eq!(false, t1.get_overlap().0);
        assert_eq!(false, t2.get_overlap().0);
    }

    #[test]
    fn does_not_contain_point_outside() {
        let t1 = TriangleComputations::new(&test_triangle1(), 1.0, 1.0);
        let t2 = TriangleComputations::new(&test_triangle2(), 0.0, 0.0);
        assert_eq!(false, t1.get_overlap().0);
        assert_eq!(false, t2.get_overlap().0);
    }

    #[test]
    fn barycentric_coordinates_1() {
        let t = TriangleComputations::new(&test_triangle1(), 0.0, 0.0);
        let (w0, w1, w2) = t.get_barycentric_coords();
        assert_eq!(1.0, w0 + w1 + w2);
        assert_eq!(w0, 1.0);
        assert_eq!(w1, 0.0);
        assert_eq!(w2, 0.0);
    }

    #[test]
    fn barycentric_coordinates_2() {
        let t = TriangleComputations::new(&test_triangle1(), 0.0, 1.0);
        let (w0, w1, w2) = t.get_barycentric_coords();
        assert_eq!(1.0, w0 + w1 + w2);
        assert_eq!(w0, 0.0);
        assert_eq!(w1, 0.0);
        assert_eq!(w2, 1.0);
    }

    #[test]
    fn barycentric_coordinates_3() {
        let t = TriangleComputations::new(&test_triangle1(), 1.0, 0.0);
        let (w0, w1, w2) = t.get_barycentric_coords();
        assert_eq!(1.0, w0 + w1 + w2);
        assert_eq!(w0, 0.0);
        assert_eq!(w1, 1.0);
        assert_eq!(w2, 0.0);
    }

    #[test]
    fn barycentric_coordinates_4() {
        let t = TriangleComputations::new(&test_triangle1(), 0.5, 0.5);
        let (w0, w1, w2) = t.get_barycentric_coords();
        assert_eq!(1.0, w0 + w1 + w2);
        assert_eq!(w0, 0.0);
        assert_eq!(w1, 0.5);
        assert_eq!(w2, 0.5);
    }

    #[test]
    fn barycentric_coordinates_5() {
        let t = TriangleComputations::new(&test_triangle1(), 0.25, 0.25);
        let (w0, w1, w2) = t.get_barycentric_coords();
        assert_eq!(1.0, w0 + w1 + w2);
        assert_eq!(w0, 0.5);
        assert_eq!(w1, 0.25);
        assert_eq!(w2, 0.25);
    }

    #[test]
    fn get_overlap_true() {
        let screen_z = 0.9509509509509508;
        let t = Triangle::new([
            Vector3::new(150.0, 450.0, screen_z),
            Vector3::new(150.0, 150.0, screen_z),
            Vector3::new(450.0, 450.0, screen_z),
        ]);

        let r = TriangleComputations::new(&t, 314.0 + 0.5, 372.0 + 0.5);
        let (overlaps, w0, w1, w2, z) = r.get_overlap();

        assert_eq!(overlaps, true);
        assert_abs_diff_eq!(w0, 0.193333, epsilon = EPSILON);
        assert_abs_diff_eq!(w1, 0.258333, epsilon = EPSILON);
        assert_abs_diff_eq!(w2, 0.548333, epsilon = EPSILON);
        assert_abs_diff_eq!(z, 0.9509509, epsilon = EPSILON);
    }
}
