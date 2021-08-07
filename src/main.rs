#![forbid(unsafe_code)]
extern crate rust_gfx;

use log::error;
use pixels::{Error, Pixels, SurfaceTexture};
use std::f32;
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::WindowBuilder;
use winit_input_helper::WinitInputHelper;

use rust_gfx::camera::default_camera;
use rust_gfx::matrix::*;
use rust_gfx::rasterize::Canvas;
use rust_gfx::shape::colour_cube;

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;
const TWO_PI: f32 = f32::consts::PI * 2.0;

fn main() -> Result<(), Error> {
    env_logger::init();
    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();

    let window = {
        let size = LogicalSize::new(WIDTH as f32, HEIGHT as f32);
        WindowBuilder::new()
            .with_title("Rust GFX")
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .unwrap()
    };

    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(WIDTH as u32, HEIGHT as u32, surface_texture)?
    };

    let canvas = Canvas::new(WIDTH, HEIGHT);
    let camera = default_camera(canvas.width as f32, canvas.height as f32);

    let mut z_buffer = vec![f32::MAX; canvas.width as usize * canvas.height as usize];

    let mut shapes = vec![colour_cube()];

    let mut rotation_rad: f32 = 0.0;
    let translation = translation(0.0, 0.0, -5.0);
    let scaling = scaling(1.0, 1.0, 1.0);

    event_loop.run(move |event, _, control_flow| {
        // Draw the current frame
        if let Event::RedrawRequested(_) = event {
            let frame = pixels.get_frame();
            frame.fill(0); // reset the frame from the previous draw

            let rotation =
                rotation_x(rotation_rad) * rotation_y(rotation_rad) * rotation_z(rotation_rad);

            let transform = scaling * rotation * translation;
            shapes.iter_mut().for_each(|s| s.transform = transform);

            canvas.draw(&camera, &shapes, frame, &mut z_buffer);

            z_buffer.fill(f32::MAX); // reset the z-buffer for the next draw

            rotation_rad = (rotation_rad + 0.02) % TWO_PI;

            if pixels
                .render()
                .map_err(|e| error!("pixels.render() failed: {}", e))
                .is_err()
            {
                *control_flow = ControlFlow::Exit;
                return;
            }
        }

        // Handle input events
        if input.update(&event) {
            // Close events
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
                return;
            }

            // Resize the window
            if let Some(size) = input.window_resized() {
                pixels.resize_surface(size.width, size.height);
            }

            window.request_redraw();
        }
    });
}
