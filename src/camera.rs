use crate::matrix::{identity, Matrix4x4};
use crate::vector::Vector3;

pub struct Camera {
    pub fov: f32,
    pub near: f32,
    pub far: f32,
    pub view: Matrix4x4<f32>, // aka world-to-camera matrix
    pub projection: Matrix4x4<f32>,
    pub view_projection: Matrix4x4<f32>,
    pub forward: Vector3<f32>,
}

impl Camera {
    pub fn new(fov: f32, near: f32, far: f32, aspect: f32) -> Self {
        let camera_to_world = identity();
        let view = camera_to_world.inverse();
        let projection = projection(fov, near, far, aspect);
        let view_projection = view * projection;
        let forward = Vector3::new(0.0, 0.0, -1.0);

        Camera {
            fov,
            near,
            far,
            view,
            projection,
            view_projection,
            forward,
        }
    }

    pub fn look_at(&mut self, from: Vector3<f32>, to: Vector3<f32>, up: Vector3<f32>) {
        let forward = (from - to).normalize();
        let right = up.normalize().cross(forward);
        let true_up = forward.cross(right);

        let camera_to_world = Matrix4x4::of([
            [right.x, right.y, right.z, 0.0],
            [true_up.x, true_up.y, true_up.z, 0.0],
            [forward.x, forward.y, forward.z, 0.0],
            [from.x, from.y, from.z, 1.0],
        ]);

        self.forward = forward;
        self.view = camera_to_world.inverse();
        self.view_projection = self.view * self.projection;
    }
}

fn projection(fov: f32, near: f32, far: f32, aspect: f32) -> Matrix4x4<f32> {
    let scale = 1.0 / (fov * 0.5).to_radians().tan();
    Matrix4x4::of([
        [scale / aspect, 0.0, 0.0, 0.0],
        [0.0, scale, 0.0, 0.0],
        [0.0, 0.0, -far / (far - near), -1.0],
        [0.0, 0.0, -far * near / (far - near), 0.0],
    ])
}

pub fn default_camera(width: f32, height: f32) -> Camera {
    let mut camera = Camera::new(45.0, 0.1, 100.0, width / height);
    camera.look_at(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(0.0, 1.0, 0.0),
    );
    camera
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EPSILON;

    #[test]
    fn projection_matrix_90_fov() {
        let m = projection(90.0, 0.1, 100.0, 1.0);
        let expected = Matrix4x4::of([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, -1.001, -1.0],
            [0.0, 0.0, -0.1001, 0.0],
        ]);
        assert_eq!(m, expected);
    }

    #[test]
    fn projection_matrix_45_fov() {
        let m = projection(45.0, 0.1, 100.0, 1.0);
        let expected = Matrix4x4::of([
            [2.41421, 0.0, 0.0, 0.0],
            [0.0, 2.41421, 0.0, 0.0],
            [0.0, 0.0, -1.001, -1.0],
            [0.0, 0.0, -0.1001, 0.0],
        ]);
        assert_eq!(m, expected);
    }

    #[test]
    fn projection_does_perspective_divide() {
        let c = Camera::new(90.0, 0.1, 100.0, 1.0);
        let v = Vector3::new(5.0, 5.0, -5.0) * c.projection;
        assert_abs_diff_eq!(v.x, 1.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v.y, 1.0, epsilon = EPSILON);
        assert_abs_diff_eq!(v.z, 0.981, epsilon = EPSILON);
    }

    #[test]
    fn projection_at_near_clipping_plain_z_is_zero() {
        let c = Camera::new(90.0, 0.1, 100.0, 1.0);
        let v = Vector3::new(5.0, 5.0, -0.1) * c.view_projection;
        assert_abs_diff_eq!(v.z, 0.0, epsilon = EPSILON);
    }

    #[test]
    fn projection_at_far_clipping_plain_z_is_one() {
        let c = Camera::new(90.0, 0.1, 100.0, 1.0);
        let v = Vector3::new(5.0, 5.0, -100.0) * c.view_projection;
        assert_abs_diff_eq!(v.z, 1.0, epsilon = EPSILON);
    }
}
