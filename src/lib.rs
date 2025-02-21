extern crate wasm_bindgen;
extern crate console_error_panic_hook;

use std::panic;

use rand::{Rng, RngCore};
use wasm_bindgen::prelude::*;
// BOUNCY BALL SIMULATOR!!!!
use lib_fluid::*;
// GRAPHICS ---------------------------------------------------------------

#[wasm_bindgen]
extern "C"
{
    #[wasm_bindgen(js_namespace = exports)]
    fn hello();
}

struct InstanceData
{
    max_instances: usize,
    instance_size_F32: usize,
    current_instance: usize,
    instance_buffer: Vec<f32>
}
static mut CIRCLE_INSTANCE_DATA: InstanceData = InstanceData {max_instances: 0, instance_size_F32: 0, current_instance: 0, instance_buffer: Vec::new() };
static mut SQUARE_INSTANCE_DATA: InstanceData = InstanceData {max_instances: 0, instance_size_F32: 0, current_instance: 0, instance_buffer: Vec::new() };

#[wasm_bindgen]
pub fn init_instance_gfx_buffers(maxCircleInstances: usize, maxSquareInstances: usize)
{
    panic::set_hook(Box::new(console_error_panic_hook::hook));

    unsafe
        {
            CIRCLE_INSTANCE_DATA.max_instances = maxCircleInstances;
            CIRCLE_INSTANCE_DATA.instance_size_F32 = 4 + 2 + 2; // FOUR FOR COLOR (R, G, B, A) TWO FOR POSITION (X, Y) TWO FOR SIZE (WIDTH, HEIGHT) 8 TOTAL FLOATS
            CIRCLE_INSTANCE_DATA.instance_buffer = vec![0.0; CIRCLE_INSTANCE_DATA.max_instances * CIRCLE_INSTANCE_DATA.instance_size_F32];
            CIRCLE_INSTANCE_DATA.current_instance = 0;

            SQUARE_INSTANCE_DATA.max_instances = maxSquareInstances;
            SQUARE_INSTANCE_DATA.instance_size_F32 = 4 + 2 + 1; // FOUR FOR COLOR (R, G, B, A) TWO FOR POS (X, Y), 1 FOR SCALE (S)
            SQUARE_INSTANCE_DATA.instance_buffer = vec![0.0; SQUARE_INSTANCE_DATA.max_instances * SQUARE_INSTANCE_DATA.instance_size_F32];
            SQUARE_INSTANCE_DATA.current_instance = 0;
        }
}

pub fn draw_circle(r: f32, g: f32, b: f32, x: f32, y: f32, s_x: f32, s_y: f32, instanceData: &mut InstanceData)
{
    let layout: [f32; 8] = [r, g, b, 1.0, x, y, s_x, s_y];
    let instance_offset: usize = instanceData.current_instance * instanceData.instance_size_F32; // instance offset corresponds to total # of floats in the buffer
    for i in 0.. 8
    {
        let index = instance_offset + i;
        instanceData.instance_buffer[index] = layout[i];
    }
    instanceData.current_instance += 1;
    if instanceData.current_instance == instanceData.max_instances
    {
        // TODO(rordon): Force draw call when buffer reaches max capacity.
        // Wrap around back to instance zero.
        instanceData.current_instance = 0;
    }
}
pub fn draw_square(r: f32, g: f32, b: f32, x: f32, y: f32, s:f32, data: &mut InstanceData)
{
    let layout: [f32; 7] = [r, g, b, 1.0, x, y, s];
    let instance_offset: usize = data.current_instance * data.instance_size_F32; // instance offset corresponds to total # of floats in the buffer
    for i in 0.. 7
    {
        let index = instance_offset + i;
        data.instance_buffer[index] = layout[i];
    }
    data.current_instance += 1;
    if data.current_instance == data.max_instances
    {
        // TODO(rordon): Force draw call when buffer reaches max capacity.
        // Wrap around back to instance zero.
        data.current_instance = 0;
    }
}

#[wasm_bindgen]
pub unsafe fn get_circle_instance_buf_ptr() -> *const f32
{
    return CIRCLE_INSTANCE_DATA.instance_buffer.as_ptr();
}
#[wasm_bindgen]
pub unsafe fn get_square_instance_buf_ptr() -> *const f32
{
    return SQUARE_INSTANCE_DATA.instance_buffer.as_ptr();
}

const INSTANCE_DATA_SIZE: usize = 4 + 2 + 2; //4+2+2

// END OF GFX. THAT'S LITERALLY IT. -----------------------------------------------------------------------------

struct Ball
{
    position: [f32; 2], // X Y
    color:    [f32; 3], // R G B
    size:     [f32; 2], // Width Height
    velocity: [f32; 2]  // VelocityX VelocityY
}
struct BounceSim
{
    areaWidth: f32,
    areaHeight: f32,
    ballz: Vec<Ball>
}

impl BounceSim {
    fn integrate_balls(&mut self, delta_time: f32)
    {
        for i in 0..self.ballz.len()
        {
            self.ballz[i].position[0] += self.ballz[i].velocity[0] * delta_time;
            self.ballz[i].position[1] += self.ballz[i].velocity[1] * delta_time;
        }
    }

    fn handle_ball_wall_collision(&mut self)
    {
        for i in 0..self.ballz.len()
        {
            // check if out of bounds on x axis
            if self.ballz[i].position[0] > (self.areaWidth/2.0) || self.ballz[i].position[0] < -(self.areaWidth/2.0)
            {
                self.ballz[i].velocity[0] *= -1.0;
            }
            if self.ballz[i].position[1] > (self.areaHeight/2.0) || self.ballz[i].position[1] < -(self.areaHeight/2.0)
            {
                self.ballz[i].velocity[1] *= -1.0;
            }
        }
    }

    unsafe fn draw_balls(&mut self)
    {
        for ball in self.ballz.iter_mut()
        {
            draw_circle(ball.color[0],
                        ball.color[1],
                        ball.color[2],
                        ball.position[0],
                        ball.position[1],
                        ball.size[0],
                        ball.size[1],
                        &mut CIRCLE_INSTANCE_DATA,
            );
        }
    }
}

// Generate a random number from 0.0 to 1.0.
pub fn randNormF32() -> f32
{
    let r: f32 = rand::thread_rng().gen_range(0.0..1.0);
    return r;
}
// BOUNCY BALL TEST SIMULATION.
static mut SIM: BounceSim = BounceSim { areaWidth:2.5f32, areaHeight:1.0f32, ballz: Vec::new() };

#[wasm_bindgen]
pub fn init_test_simulation()
{
    unsafe
    {
        // make some ballz
        for i in 1..100
        {
            SIM.ballz.push(Ball{ position: [(SIM.areaWidth / 2.0)*randNormF32(), (SIM.areaHeight / 2.0)*randNormF32()],
                                       color:[randNormF32(), randNormF32(), randNormF32()],
                                       size: [0.1,0.1],
                                       velocity: [(randNormF32()*2.0)-1.0, (randNormF32()*2.0)-1.0]});
        }
    }
}

// ---------------------------------------------------------------------------

unsafe fn draw_simulation(fluid: &FlipFluid, particle_radius: f32)
{
    // Draw particles
    // Reset Instance Buffer
    CIRCLE_INSTANCE_DATA.current_instance = 0;
    for i in 0..fluid.num_particles
    {
        let index = i as usize;
        let pos_x = fluid.particle_pos[(index * 2) + 0];
        let pos_y = fluid.particle_pos[(index * 2) + 1];

        let r = fluid.particle_colour[(index * 3) + 0];
        let g = fluid.particle_colour[(index * 3) + 1];
        let b = fluid.particle_colour[(index * 3) + 2];

        draw_circle(r, 0.5, 1.0, pos_x, pos_y, particle_radius * 2.0, particle_radius * 2.0, &mut CIRCLE_INSTANCE_DATA)
    }
    SQUARE_INSTANCE_DATA.current_instance = 0;
    for x in 0..fluid.f_num_x{
        for y in 0..fluid.f_num_y
        {
            let index = ((x * fluid.f_num_y + y) * 3) as usize;
            let r = fluid.cell_colour[index + 0];
            let g = fluid.cell_colour[index + 1];
            let b = fluid.cell_colour[index + 2];

            draw_square(r, g, b,
                        (x as f32 + 0.5) * fluid.cell_size,
                        (y as f32 + 0.5) * fluid.cell_size,
                        fluid.cell_size * 0.8,
                        &mut SQUARE_INSTANCE_DATA
            );
        }
    }
}

#[wasm_bindgen]
pub struct SimulationHandler
{
    scene: lib_fluid::Scene,
}

#[wasm_bindgen]
impl SimulationHandler
{
    // Create scene struct.
    #[wasm_bindgen(constructor)]
    pub fn new(num_x: i32, num_y: i32, density: f32, width: f32, height: f32, cell_size: f32, particle_radius: f32, max_particles: i32) -> Self
    {
        let mut scene = lib_fluid::Scene::new(
            density,
            width,
            height,
            cell_size,
            particle_radius,
            max_particles
        );

        create_particles(&mut scene.fluid, num_x, num_y);
        setup_grid(&mut scene.fluid);
        // TODO: finish this if we have time.
        // params are for mouse_x, mouse_y, and reset
        //set_obstacle(&mut scene, 3.0, 2.0, true, num_x, num_y);

        return SimulationHandler { scene };
    }

    #[wasm_bindgen]
    pub fn update(&mut self, delta_time: f32, mouse_x:f32, mouse_y:f32, gravity:f32)
    {
        // TODO: add pausing to scene
        if (true)
        {
            self.scene.fluid.simulate(
                self.scene.dt,
                gravity,
                self.scene.flip_ratio,
                self.scene.num_pressure_iters,
                self.scene.num_particle_iters,
                self.scene.over_relaxation,
                self.scene.compensate_drift,
                self.scene.separate_particles,
                self.scene.obstacle_x,
                self.scene.obstacle_y,
                self.scene.obstacle_radius
            );
        }

        unsafe { draw_simulation(&self.scene.fluid, 0.018/2.0); }

        // hello();
        unsafe
        {
            SIM.integrate_balls(delta_time);
            SIM.handle_ball_wall_collision();
            //SIM.draw_balls();
        }
    }
}