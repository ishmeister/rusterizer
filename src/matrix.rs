use crate::vector::Vector3;
use crate::EPSILON;
use std::f32;
use std::ops;
use std::ops::{Index, IndexMut};

#[derive(Clone, Copy, Debug)]
pub struct Matrix4x4<T> {
    elems: [[T; 4]; 4],
}

impl Default for Matrix4x4<f32> {
    fn default() -> Self {
        Self::new()
    }
}

impl Matrix4x4<f32> {
    pub fn new() -> Self {
        Matrix4x4::<f32> {
            elems: [[0.0f32; 4]; 4],
        }
    }

    pub fn of(elems: [[f32; 4]; 4]) -> Self {
        Matrix4x4::<f32> { elems }
    }

    fn mul_point(&self, p: Vector3<f32>) -> Vector3<f32> {
        let mut x = p.x * self[0][0] + p.y * self[1][0] + p.z * self[2][0] + self[3][0];
        let mut y = p.x * self[0][1] + p.y * self[1][1] + p.z * self[2][1] + self[3][1];
        let mut z = p.x * self[0][2] + p.y * self[1][2] + p.z * self[2][2] + self[3][2];
        let w = p.x * self[0][3] + p.y * self[1][3] + p.z * self[2][3] + self[3][3];

        // homogenous coordinate
        // the last column is almost always (0, 0, 0, 1) except when projection matrices are used
        // so the point needs to be normalised so that w = 1 again
        if abs_diff_ne!(w, 1.0, epsilon = EPSILON) && abs_diff_ne!(w, 0.0, epsilon = EPSILON) {
            x /= w;
            y /= w;
            z /= w;
        }

        Vector3::new(x, y, z)
    }

    /// optimised 4x4 matrix inverse function based on: http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
    /// code from: https://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform
    pub fn inverse(&self) -> Matrix4x4<f32> {
        let a = self;

        let s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1];
        let s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
        let s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3];
        let s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
        let s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3];
        let s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];

        let c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3];
        let c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
        let c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2];
        let c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
        let c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2];
        let c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];

        let det = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
        assert_abs_diff_ne!(det, 0.0, epsilon = EPSILON);

        let invdet = 1.0 / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

        let mut b = Matrix4x4::default();

        b[0][0] = (a[1][1] * c5 - a[1][2] * c4 + a[1][3] * c3) * invdet;
        b[0][1] = (-a[0][1] * c5 + a[0][2] * c4 - a[0][3] * c3) * invdet;
        b[0][2] = (a[3][1] * s5 - a[3][2] * s4 + a[3][3] * s3) * invdet;
        b[0][3] = (-a[2][1] * s5 + a[2][2] * s4 - a[2][3] * s3) * invdet;

        b[1][0] = (-a[1][0] * c5 + a[1][2] * c2 - a[1][3] * c1) * invdet;
        b[1][1] = (a[0][0] * c5 - a[0][2] * c2 + a[0][3] * c1) * invdet;
        b[1][2] = (-a[3][0] * s5 + a[3][2] * s2 - a[3][3] * s1) * invdet;
        b[1][3] = (a[2][0] * s5 - a[2][2] * s2 + a[2][3] * s1) * invdet;

        b[2][0] = (a[1][0] * c4 - a[1][1] * c2 + a[1][3] * c0) * invdet;
        b[2][1] = (-a[0][0] * c4 + a[0][1] * c2 - a[0][3] * c0) * invdet;
        b[2][2] = (a[3][0] * s4 - a[3][1] * s2 + a[3][3] * s0) * invdet;
        b[2][3] = (-a[2][0] * s4 + a[2][1] * s2 - a[2][3] * s0) * invdet;

        b[3][0] = (-a[1][0] * c3 + a[1][1] * c1 - a[1][2] * c0) * invdet;
        b[3][1] = (a[0][0] * c3 - a[0][1] * c1 + a[0][2] * c0) * invdet;
        b[3][2] = (-a[3][0] * s3 + a[3][1] * s1 - a[3][2] * s0) * invdet;
        b[3][3] = (a[2][0] * s3 - a[2][1] * s1 + a[2][2] * s0) * invdet;

        b
    }
}

impl Index<usize> for Matrix4x4<f32> {
    type Output = [f32; 4];
    fn index(&self, index: usize) -> &Self::Output {
        &self.elems[index]
    }
}

impl IndexMut<usize> for Matrix4x4<f32> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elems[index]
    }
}

impl ops::Mul<Matrix4x4<f32>> for Matrix4x4<f32> {
    type Output = Matrix4x4<f32>;
    fn mul(self, rhs: Matrix4x4<f32>) -> Matrix4x4<f32> {
        let mut m = Matrix4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] = self[i][0] * rhs[0][j]
                    + self[i][1] * rhs[1][j]
                    + self[i][2] * rhs[2][j]
                    + self[i][3] * rhs[3][j];
            }
        }
        m
    }
}

impl ops::Mul<Vector3<f32>> for Matrix4x4<f32> {
    type Output = Vector3<f32>;
    fn mul(self, rhs: Vector3<f32>) -> Vector3<f32> {
        self.mul_point(rhs)
    }
}

// for row-major multiplications with a point
impl ops::Mul<Matrix4x4<f32>> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn mul(self, rhs: Matrix4x4<f32>) -> Vector3<f32> {
        rhs.mul_point(self)
    }
}

impl PartialEq for Matrix4x4<f32> {
    fn eq(&self, rhs: &Self) -> bool {
        for i in 0..4 {
            for j in 0..4 {
                if !abs_diff_eq!(self[i][j], rhs[i][j], epsilon = EPSILON) {
                    return false;
                }
            }
        }
        true
    }
}

pub const fn identity() -> Matrix4x4<f32> {
    Matrix4x4::<f32> {
        elems: [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    }
}

pub fn translation(x: f32, y: f32, z: f32) -> Matrix4x4<f32> {
    Matrix4x4::of([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [x, y, z, 1.0],
    ])
}

pub fn scaling(x: f32, y: f32, z: f32) -> Matrix4x4<f32> {
    Matrix4x4::of([
        [x, 0.0, 0.0, 0.0],
        [0.0, y, 0.0, 0.0],
        [0.0, 0.0, z, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
}

pub fn rotation_x(r: f32) -> Matrix4x4<f32> {
    let cos_r = r.cos();
    let sin_r = r.sin();
    Matrix4x4::of([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, cos_r, sin_r, 0.0],
        [0.0, -sin_r, cos_r, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
}

pub fn rotation_y(r: f32) -> Matrix4x4<f32> {
    let cos_r = r.cos();
    let sin_r = r.sin();
    Matrix4x4::of([
        [cos_r, 0.0, -sin_r, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [sin_r, 0.0, cos_r, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
}

pub fn rotation_z(r: f32) -> Matrix4x4<f32> {
    let cos_r = r.cos();
    let sin_r = r.sin();
    Matrix4x4::of([
        [cos_r, sin_r, 0.0, 0.0],
        [-sin_r, cos_r, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix4x4_identity() {
        let m = identity();
        assert_eq!(m[0][0], 1.0);
        assert_eq!(m[0][1], 0.0);
        assert_eq!(m[0][2], 0.0);
        assert_eq!(m[0][3], 0.0);
        assert_eq!(m[1][0], 0.0);
        assert_eq!(m[1][1], 1.0);
        assert_eq!(m[1][2], 0.0);
        assert_eq!(m[1][3], 0.0);
        assert_eq!(m[2][0], 0.0);
        assert_eq!(m[2][1], 0.0);
        assert_eq!(m[2][2], 1.0);
        assert_eq!(m[2][3], 0.0);
        assert_eq!(m[3][0], 0.0);
        assert_eq!(m[3][1], 0.0);
        assert_eq!(m[3][2], 0.0);
        assert_eq!(m[3][3], 1.0);
    }

    #[test]
    fn matrix4x4_identity_multiply() {
        let a = Matrix4x4::of([
            [0.0, 1.0, 2.0, 4.0],
            [1.0, 2.0, 4.0, 8.0],
            [2.0, 4.0, 8.0, 16.0],
            [4.0, 8.0, 16.0, 32.0],
        ]);

        let b = a * identity();

        assert_eq!(a, b);
    }

    #[test]
    fn matrix4x4_multiply() {
        let a = Matrix4x4::of([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);

        let b = Matrix4x4::of([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);

        let c = a * b;

        let expected = Matrix4x4::of([
            [20.0, 22.0, 50.0, 48.0],
            [44.0, 54.0, 114.0, 108.0],
            [40.0, 58.0, 110.0, 102.0],
            [16.0, 26.0, 46.0, 42.0],
        ]);

        assert_eq!(c, expected);
    }

    #[test]
    fn multiply_point_by_translation_matrix() {
        let m = translation(5.0, -3.0, 2.0);
        let p1 = Vector3::new(-3.0, 4.0, 5.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(2.0, 1.0, 7.0), p2);
    }

    #[test]
    fn multiply_point_by_scaling_matrix() {
        let m = scaling(2.0, 3.0, 4.0);
        let p1 = Vector3::new(-4.0, 6.0, 8.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(-8.0, 18.0, 32.0), p2);
    }

    #[test]
    fn multiply_point_by_rotation_x_matrix_90_degrees() {
        let m = rotation_x((90.0f32).to_radians());
        let p1 = Vector3::new(0.0, 1.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(0.0, 0.0, 1.0), p2);
    }

    #[test]
    fn multiply_point_by_rotation_x_matrix_45_degrees() {
        let m = rotation_x((45.0f32).to_radians());
        let p1 = Vector3::new(0.0, 1.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(0.0, 0.7071, 0.7071), p2);
    }

    #[test]
    fn multiply_point_by_rotation_y_matrix_90_degrees() {
        let m = rotation_y((90.0f32).to_radians());
        let p1 = Vector3::new(1.0, 0.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(0.0, 0.0, -1.0), p2);
    }

    #[test]
    fn multiply_point_by_rotation_y_matrix_45_degrees() {
        let m = rotation_y((45.0f32).to_radians());
        let p1 = Vector3::new(1.0, 0.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(0.7071, 0.0, -0.7071), p2);
    }

    #[test]
    fn multiply_point_by_rotation_z_matrix_90_degrees() {
        let m = rotation_z((90.0f32).to_radians());
        let p1 = Vector3::new(0.0, 1.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(-1.0, 0.0, 0.0), p2);
    }

    #[test]
    fn multiply_point_by_rotation_z_matrix_45_degrees() {
        let m = rotation_z((45.0f32).to_radians());
        let p1 = Vector3::new(0.0, 1.0, 0.0);
        let p2 = p1 * m;
        assert_eq!(Vector3::new(-0.7071, 0.7071, 0.0), p2);
    }

    #[test]
    fn multiply_point_by_combined_rotation_xyz_matrix_degrees() {
        let rotation = (90.0f32).to_radians();
        let xr = rotation_x(rotation);
        let yr = rotation_y(rotation);
        let zr = rotation_z(-rotation);
        let p1 = Vector3::new(0.0, 1.0, 0.0);
        let p2 = p1 * xr * yr * zr;
        assert_eq!(Vector3::new(0.0, -1.0, 0.0), p2);
    }

    #[test]
    fn matrix_inverse() {
        let a = Matrix4x4::of([
            [-5.0, 2.0, 6.0, -8.0],
            [1.0, -5.0, 1.0, 8.0],
            [7.0, 7.0, -6.0, -7.0],
            [1.0, -3.0, 7.0, 4.0],
        ]);

        let b = a.inverse();

        let expected = Matrix4x4::of([
            [0.21805, 0.45113, 0.24060, -0.04511],
            [-0.80827, -1.45677, -0.44361, 0.520680],
            [-0.07895, -0.22368, -0.05263, 0.19737],
            [-0.52256, -0.81391, -0.30075, 0.30639],
        ]);

        assert_eq!(b, expected);
    }

    #[test]
    fn matrix_inverse_2() {
        let a = Matrix4x4::of([
            [8.0, -5.0, 9.0, 2.0],
            [7.0, 5.0, 6.0, 1.0],
            [-6.0, 0.0, 9.0, 6.0],
            [-3.0, 0.0, -9.0, -4.0],
        ]);

        let b = a.inverse();

        let expected = Matrix4x4::of([
            [-0.15385, -0.15385, -0.28205, -0.53846],
            [-0.07692, 0.12308, 0.02564, 0.03077],
            [0.35897, 0.35897, 0.43590, 0.92308],
            [-0.69231, -0.69231, -0.76923, -1.92308],
        ]);

        assert_eq!(b, expected);
    }

    #[test]
    fn matrix_inverse_3() {
        let a = Matrix4x4::of([
            [9.0, 3.0, 0.0, 9.0],
            [-5.0, -2.0, -6.0, -3.0],
            [-4.0, 9.0, 6.0, 4.0],
            [-7.0, 6.0, 6.0, 2.0],
        ]);

        let b = a.inverse();

        let expected = Matrix4x4::of([
            [-0.04074, -0.07778, 0.14444, -0.22222],
            [-0.07778, 0.03333, 0.36667, -0.33333],
            [-0.02901, -0.14630, -0.10926, 0.12963],
            [0.17778, 0.06667, -0.26667, 0.33333],
        ]);

        assert_eq!(b, expected);
    }

    #[test]
    fn matrix_multiply_product_by_inverse() {
        let a = Matrix4x4::of([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);

        let b = Matrix4x4::of([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);

        let c = a * b;

        assert_eq!(a, c * b.inverse());
    }

    #[test]
    #[should_panic]
    fn matrix_non_invertible() {
        let a = Matrix4x4::of([
            [-4.0, 2.0, -2.0, -3.0],
            [9.0, 6.0, 2.0, 6.0],
            [-4.0, 4.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        a.inverse();
    }
}
