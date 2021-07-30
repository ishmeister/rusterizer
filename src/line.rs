use super::vector::Vector3;
use super::rasterize::Canvas;
use super::rasterize::Colour;

pub struct Line {
    pub p1: Vector3<f32>,
    pub p2: Vector3<f32>,
    pub colour: Colour
}

impl Line {
    pub fn new(p1: Vector3<f32>, p2: Vector3<f32>, colour: Colour) -> Self {
        Line {
            p1, p2, colour
        }
    }

    pub fn draw_bres(&self, x0: i32, y0: i32, x1: i32, y1: i32, canvas: &Canvas, frame: &mut [u8]) {
        let dx = (x1 - x0).abs();
        let dy = (y1 - y0).abs(); 

        if dy < dx {
            if x0 > x1 {
                self.draw_bres_step_x(x1, y1, x0, y0, canvas, frame);
            }
            else {
                self.draw_bres_step_x(x0, y0, x1, y1, canvas, frame);
            }
        }
        else {
            if y0 > y1 {
                self.draw_bres_step_y(x1, y1, x0, y0, canvas, frame);
            }
            else {
                self.draw_bres_step_y(x0, y0, x1, y1, canvas, frame);
            }
        }
    }

    fn draw_bres_step_x(&self, x0: i32, y0: i32, x1: i32, y1: i32, canvas: &Canvas, frame: &mut [u8]) {
        let dx = x1 - x0;
        let mut dy = y1 - y0;
        let mut y_inc = 1;
        
        if dy < 0 {
            y_inc = -1;
            dy = -dy;
        }

        let two_dy = 2 * dy;
        let two_dy_dx = 2 * (dy - dx);
        let mut p = two_dy - dx;
        let mut y = y0;

        for x in x0..x1 {
            canvas.set_pixel(x as u32, y as u32, self.colour, frame);

            if p > 0 {
                y += y_inc;
                p += two_dy_dx;
            }
            else {
                p += two_dy;
            }
        }
    }

    fn draw_bres_step_y(&self, x0: i32, y0: i32, x1: i32, y1: i32, canvas: &Canvas, frame: &mut [u8]) {
        let mut dx = x1 - x0;
        let dy = y1 - y0;
        let mut x_inc = 1;
        
        if dx < 0 {
            x_inc = -1;
            dx = -dx;
        }

        let two_dx = 2 * dx;
        let two_dx_dy = 2 * (dx - dy);
        let mut p = two_dx - dy;
        let mut x = x0;

        for y in y0..y1 {
            canvas.set_pixel(x as u32, y as u32, self.colour, frame);

            if p > 0 {
                x += x_inc;
                p += two_dx_dy;
            }
            else {
                p += two_dx;
            }
        }
    }
}
