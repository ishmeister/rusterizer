use crate::camera::Camera;
use crate::matrix::Matrix4x4;
use crate::shape::Shape;
use crate::triangle::*;
use crate::vector::Vector3;

pub const WHITE: Colour = Colour {
    r: 1.0,
    g: 1.0,
    b: 1.0,
    a: 1.0,
};
pub const GREY: Colour = Colour {
    r: 0.5,
    g: 0.5,
    b: 0.5,
    a: 1.0,
};
pub const RED: Colour = Colour {
    r: 1.0,
    g: 0.0,
    b: 0.0,
    a: 1.0,
};
pub const GREEN: Colour = Colour {
    r: 0.0,
    g: 1.0,
    b: 0.0,
    a: 1.0,
};
pub const BLUE: Colour = Colour {
    r: 0.0,
    g: 0.0,
    b: 1.0,
    a: 1.0,
};

pub trait Draw {
    fn draw(&self, camera: &Camera, canvas: &Canvas, frame: &mut [u8], z_buffer: &mut [f32]);
}

#[derive(Clone, Debug)]
pub struct Canvas {
    pub width: u32,
    pub height: u32,
}

#[derive(Clone, Copy, Debug)]
pub struct Colour {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    pub a: f32,
}

impl Canvas {
    pub fn new(width: u32, height: u32) -> Self {
        assert!(width > 0 && height > 0);
        Self { width, height }
    }

    pub fn draw(
        &self,
        camera: &Camera,
        shapes: &[Shape],
        frame: &mut [u8],
        z_buffer: &mut [f32],
    ) {
        for shape in shapes.iter() {
            shape.draw(camera, self, frame, z_buffer);
        }
    }

    pub fn set_pixel(&self, x: u32, y: u32, colour: Colour, frame: &mut [u8]) {
        if x < self.width && y < self.height {
            let index = 4 * (x + y * self.width) as usize;
            let pixel = [
                (255.0 * colour.r) as u8,
                (255.0 * colour.g) as u8,
                (255.0 * colour.b) as u8,
                (255.0 * colour.a) as u8,
            ];
            frame[index..index + 4].copy_from_slice(&pixel);
        }
    }
}

// one step conversion from screen -> NDC -> raster space
pub fn to_raster_space(x: f32, y: f32, width: u32, height: u32) -> (f32, f32) {
    let x_pixel = (x + 1.0) * 0.5 * width as f32;
    let y_pixel = (1.0 - (y + 1.0) * 0.5) * height as f32;
    ((x_pixel + 0.5).floor(), (y_pixel + 0.5).floor())
}

impl Draw for Vector3<f32> {
    fn draw(&self, _: &Camera, canvas: &Canvas, frame: &mut [u8], _: &mut [f32]) {
        let (x, y) = to_raster_space(self.x, self.y, canvas.width, canvas.height);
        canvas.set_pixel(x as u32, y as u32, WHITE, frame);
    }
}

impl Draw for Shape {
    fn draw(&self, camera: &Camera, canvas: &Canvas, frame: &mut [u8], z_buffer: &mut [f32]) {
        for t in self.triangles.iter() {
            draw_triangle(t, self.transform, camera, canvas, frame, z_buffer);
        }
    }
}

fn draw_triangle(
    triangle: &Triangle<f32>,
    to_world: Matrix4x4<f32>,
    camera: &Camera,
    canvas: &Canvas,
    frame: &mut [u8],
    z_buffer: &mut [f32],
) {
    let to_screen = to_world * camera.view_projection;
    let mut colours = triangle.colours;
    let mut raster_points = [Vector3::new(0.0, 0.0, 0.0); 3];
    let mut z_clipped = false;

    for i in 0..3 {
        let screen_point = triangle.points[i] * to_screen;
        let (x, y) = to_raster_space(screen_point.x, screen_point.y, canvas.width, canvas.height);

        raster_points[i].x = x;
        raster_points[i].y = y;
        raster_points[i].z = screen_point.z;

        colours[i].r /= screen_point.z;
        colours[i].g /= screen_point.z;
        colours[i].b /= screen_point.z;

        z_clipped &= screen_point.z < 0.0 || screen_point.z > 1.0;
    }

    if z_clipped {
        // skip triangles that are fully outside near and far clipping planes
        return;
    };

    let raster_triangle = Triangle::new(raster_points);
    let max_w = canvas.width as i32;
    let max_h = canvas.height as i32;

    let (xmin, xmax, ymin, ymax) = get_triangle_2d_bounding_box(&raster_triangle);

    if xmin == xmax || ymin == ymax {
        // skip bounding boxes with zero dimensions
        return;
    }

    if xmin > max_w || xmax < 0 || ymin > max_h || ymax < 0 {
        // skip bounding boxes that are completely off screen
        return;
    }

    if !triangle.double_side && raster_triangle.normal.dot(camera.forward) > 0.0 {
        // ignore triangle back faces (facing away from camera)
        return;
    }

    let x_start = xmin.max(0) as u32;
    let x_end = xmax.min(max_w) as u32;
    let y_start = ymin.max(0) as u32;
    let y_end = ymax.min(max_h) as u32;

    // Initialise the point to test at the top of the bounding box and then
    // step over each pixel testing for overlap. Assumes linear iteration.
    let mut raster_comps =
        TriangleComputations::new(&raster_triangle, x_start as f32 + 0.5, y_start as f32 + 0.5);

    for y in y_start..y_end {
        for x in x_start..x_end {
            let (overlaps, w0, w1, w2, z) = raster_comps.get_overlap();
            let visible = (0.0..=1.0).contains(&z);

            if overlaps && visible {
                let zbuffer_idx = (x + y * canvas.width) as usize;

                if z < z_buffer[zbuffer_idx] {
                    let pixel_colour = Colour {
                        r: (colours[0].r * w0 + colours[1].r * w1 + colours[2].r * w2) * z,
                        g: (colours[0].g * w0 + colours[1].g * w1 + colours[2].g * w2) * z,
                        b: (colours[0].b * w0 + colours[1].b * w1 + colours[2].b * w2) * z,
                        a: 1.0,
                    };

                    canvas.set_pixel(x, y, pixel_colour, frame);
                    z_buffer[zbuffer_idx] = z;
                }
            }

            raster_comps.step_x();
        }

        raster_comps.step_y();
    }
}

fn get_triangle_2d_bounding_box(t: &Triangle<f32>) -> (i32, i32, i32, i32) {
    let xmin = t.points[0].x.min(t.points[1].x.min(t.points[2].x)) as i32;
    let xmax = t.points[0].x.max(t.points[1].x.max(t.points[2].x)) as i32;
    let ymin = t.points[0].y.min(t.points[1].y.min(t.points[2].y)) as i32;
    let ymax = t.points[0].y.max(t.points[1].y.max(t.points[2].y)) as i32;
    (xmin, xmax, ymin, ymax)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_raster_space_screen_origin_is_canvas_center() {
        let (x, y) = to_raster_space(0.0, 0.0, 100, 100);
        assert_eq!(x, 50.0);
        assert_eq!(y, 50.0);
    }

    #[test]
    fn to_raster_space_screen_y_is_inverted_on_canvas() {
        let (x, y) = to_raster_space(-0.5, -0.5, 100, 100);
        assert_eq!(x, 25.0);
        assert_eq!(y, 75.0);
    }

    #[test]
    fn to_raster_space_screen_top_left_is_valid_pixel_index() {
        let (x, y) = to_raster_space(-1.0, 1.0, 100, 100);
        assert_eq!(x, 0.0);
        assert_eq!(y, 0.0);
    }
}
