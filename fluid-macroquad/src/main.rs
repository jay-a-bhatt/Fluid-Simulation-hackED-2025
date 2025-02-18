use macroquad::prelude::*;
use lib_fluid;
use lib_fluid::{create_particles, FlipFluid, Scene};
use libm::floorf;
use libm::sqrtf;

#[macroquad::main("Fluid Macroquad")]
async fn main() {
    println!("Hello, world!");

    let mut scene = init_scene(1000,1000);
    let zoom = 0.5;
    request_new_screen_size(1000.0,1000.0);
    loop {
        let aspect = screen_width() / screen_height();
        clear_background(Color::new(0.2,0.2, 0.3, 1.0));
        set_camera(&Camera2D {
            zoom: vec2(zoom, zoom),
            target: vec2(0.8, 1.2),
            ..Default::default()
        });
        draw_sim(&scene.fluid);
        scene.fluid.simulate(1.0/60.0,
                             scene.gravity,
                             scene.flip_ratio,
                             scene.num_pressure_iters,
                             scene.num_particle_iters,
                             scene.over_relaxation,
                             scene.compensate_drift,
                             scene.separate_particles,
                             scene.obstacle_x,
                             scene.obstacle_y,
                             scene.obstacle_radius
        );
        next_frame().await
    }
}

pub fn draw_sim(fluid: &FlipFluid)
{
    for i in 0..fluid.num_particles
    {
        let index = i as usize;
        let pos_x = fluid.particle_pos[(index * 2) + 0];
        let pos_y = fluid.particle_pos[(index * 2) + 1];

        let r = fluid.particle_colour[(index * 3) + 0];
        let g = fluid.particle_colour[(index * 3) + 1];
        let b = fluid.particle_colour[(index * 3) + 2];

        draw_rectangle(pos_x, pos_y, 0.009*2.0, 0.009*2.0,Color::new(0.0,0.5,1.0,1.0));
    }

}


pub fn init_scene(canvas_height: i32, canvas_width: i32) -> Scene
{
    let simHeight: f32 = 3.0;
    let canvas_scale:f32 = canvas_height as f32 / simHeight;
    let simWidth:f32 = canvas_width as f32 / canvas_scale;

    // NOTE(Doesnt change)
    let tank_width = 1.0 * simWidth;
    let tank_height = 1.0 * simHeight;

    let relWaterHeight = 0.8;
    let relWaterWidth = 0.6;
    let density = 1000.0;

    let res = 100.0;
    let cell_size = tank_height / res;

    // Particle Radius is with respect to cell size.
    let particle_radius = 0.3 * cell_size;

    let dx = 2.0 * particle_radius;
    let dy = sqrtf(3.0) / 2.0 * dx;

    // Number of potential particles on the X axis.
    let num_x = floorf((relWaterWidth * tank_width - 2.0 * cell_size - 2.0 * particle_radius) / dx) as i32;
    // Number of potential particles on the y axis.
    let num_y = floorf((relWaterHeight * tank_height - 2.0 * cell_size - 2.0 * particle_radius) / dy) as i32;

    // Maximum possible particles in our simulation.
    let max_particles = num_x * num_y;

    let mut scene = Scene::new(density, tank_width, tank_height, cell_size, particle_radius, max_particles);
    create_particles(&mut scene.fluid, num_x, num_y);
    return scene;

}