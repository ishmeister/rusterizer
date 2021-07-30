use crate::EPSILON;
use std::f32;
use std::ops;

#[derive(Copy, Clone, Debug)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl Vector3<f32> {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Vector3::<f32> { x, y, z }
    }

    pub fn from_array(coords: [f32; 3]) -> Self {
        Vector3::<f32> {
            x: coords[0],
            y: coords[1],
            z: coords[2],
        }
    }

    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Vector3<f32> {
        let length = self.dot(*self);
        if length > 0.0 {
            let inv_length = 1.0 / length.sqrt();
            Vector3::new(
                self.x * inv_length,
                self.y * inv_length,
                self.z * inv_length,
            )
        } else {
            *self
        }
    }

    pub fn dot(&self, rhs: Vector3<f32>) -> f32 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn cross(&self, rhs: Vector3<f32>) -> Vector3<f32> {
        Vector3::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }
}

impl PartialEq for Vector3<f32> {
    fn eq(&self, rhs: &Self) -> bool {
        abs_diff_eq!(self.x, rhs.x, epsilon = EPSILON)
            && abs_diff_eq!(self.y, rhs.y, epsilon = EPSILON)
            && abs_diff_eq!(self.z, rhs.z, epsilon = EPSILON)
    }
}

impl ops::Sub<Vector3<f32>> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn sub(self, rhs: Vector3<f32>) -> Vector3<f32> {
        Vector3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl ops::Add<Vector3<f32>> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn add(self, rhs: Vector3<f32>) -> Vector3<f32> {
        Vector3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl ops::Mul<Vector3<f32>> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn mul(self, rhs: Vector3<f32>) -> Vector3<f32> {
        Vector3::new(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
    }
}

impl ops::Mul<f32> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn mul(self, s: f32) -> Vector3<f32> {
        Vector3::new(self.x * s, self.y * s, self.z * s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector3_length_x() {
        let v = Vector3::new(1.0, 0.0, 0.0);
        let l = v.length();
        assert_abs_diff_eq!(l, 1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_length_y() {
        let v = Vector3::new(0.0, 1.0, 0.0);
        let l = v.length();
        assert_abs_diff_eq!(l, 1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_length_z() {
        let v = Vector3::new(0.0, 0.0, 1.0);
        let l = v.length();
        assert_abs_diff_eq!(l, 1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_length_xyz() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let l = v.length();
        assert_abs_diff_eq!(l, 3.74165, epsilon = EPSILON);
    }

    #[test]
    fn vector3_length_minus_xyz() {
        let v = Vector3::new(-1.0, -2.0, -3.0);
        let l = v.length();
        assert_abs_diff_eq!(l, 3.74165, epsilon = EPSILON);
    }

    #[test]
    fn vector3_normalize() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let n = v.normalize();
        assert_abs_diff_eq!(n.x, 0.26726, epsilon = EPSILON);
        assert_abs_diff_eq!(n.y, 0.5345, epsilon = EPSILON);
        assert_abs_diff_eq!(n.z, 0.80178, epsilon = EPSILON);
    }

    #[test]
    fn vector3_normalize_one_axis_is_1() {
        let v = Vector3::new(4.0, 0.0, 0.0);
        let n = v.normalize();
        assert_abs_diff_eq!(n.x, 1.0, epsilon = EPSILON);
        assert_abs_diff_eq!(n.y, 0.0, epsilon = EPSILON);
        assert_abs_diff_eq!(n.z, 0.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_normalize_length_is_1() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let l = v.normalize().length();
        assert_abs_diff_eq!(l, 1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_dot() {
        let v1 = Vector3::new(1.0, 2.0, 3.0);
        let v2 = Vector3::new(1.0, 2.0, 3.0);
        let d = v1.dot(v2);
        assert_abs_diff_eq!(d, 14.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_dot_orthogonal() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);
        let d = v1.dot(v2);
        assert_abs_diff_eq!(d, 0.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_dot_opposite() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(-1.0, 0.0, 0.0);
        let d = v1.dot(v2);
        assert_abs_diff_eq!(d, -1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_cross_product_is_orthogonal() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);
        let v3 = v1.cross(v2);
        assert_abs_diff_eq!(v3.x, 0.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v3.y, 0.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v3.z, 1.0, epsilon = EPSILON);
    }

    #[test]
    fn vector3_cross_product_reverse_is_negated() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);
        let v3 = v2.cross(v1);
        assert_abs_diff_eq!(v3.x, 0.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v3.y, 0.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v3.z, -1.0, epsilon = EPSILON);
    }
}
