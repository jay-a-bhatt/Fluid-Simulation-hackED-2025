extern crate wasm_bindgen;
use rand::{Rng, RngCore};
use wasm_bindgen::prelude::*;
// BOUNCY BALL SIMULATOR!!!!
use lib_fluid::*;
// GRAPHICS ---------------------------------------------------------------
const MAX_INSTANCES: usize = 16384;
// NUMBER OF FLOATS PER INSTANCE

#[wasm_bindgen]
extern "C"
{
    #[wasm_bindgen(js_namespace = exports)]
    fn hello();
}

const INSTANCE_SIZE_F32: usize = 4+2+2; // FOUR FOR COLOR (R, G, B, A) TWO FOR POSITION (X, Y) TWO FOR SIZE (WIDTH, HEIGHT) 8 TOTAL FLOATS
static mut INSTANCE_GFX_BUFFER: [f32; MAX_INSTANCES* INSTANCE_SIZE_F32] = [0.0; INSTANCE_SIZE_F32 * MAX_INSTANCES];
static mut CURRENT_INSTANCE: usize = 0;
pub fn draw_circle(r: f32, g: f32, b: f32, x: f32, y: f32, s_x: f32, s_y: f32)
{
    unsafe
    {
        let array: [f32; 8] = [r, g, b, 1.0, x, y, s_x, s_y];
        let instance_offset: usize = CURRENT_INSTANCE * INSTANCE_SIZE_F32; // instance offset corresponds to total # of floats in the buffer
        for i in 0.. 8
        {
            let index = instance_offset + i;
            INSTANCE_GFX_BUFFER[index] = array[i];
        }
        CURRENT_INSTANCE += 1;
        // RESET BUFFER FOR NOW
        // TODO(rordon): Force draw call when buffer reaches max capacity.
        if CURRENT_INSTANCE == MAX_INSTANCES
        {
            // force draw call
            CURRENT_INSTANCE = 0;
        }
    }
}
#[wasm_bindgen]
pub fn get_instance_buffer_ptr() -> *const f32 {
    let p: *const f32;
    unsafe
    {
        p = INSTANCE_GFX_BUFFER.as_ptr();
    }
    return p;
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

    fn draw_balls(&mut self)
    {
        for ball in self.ballz.iter_mut()
        {
            draw_circle(ball.color[0],
                        ball.color[1],
                        ball.color[2],
                        ball.position[0],
                        ball.position[1],
                        ball.size[0],
                        ball.size[1]
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

fn draw_simulation(fluid: &FlipFluid, particle_radius: f32)
{
    unsafe { CURRENT_INSTANCE = 0 };
    for i in 0..fluid.num_particles
    {
        let index = i as usize;
        let pos_x = fluid.particle_pos[(index * 2) + 0];
        let pos_y = fluid.particle_pos[(index * 2) + 1];

        let r = fluid.particle_colour[(index * 3) + 0];
        let g = fluid.particle_colour[(index * 3) + 1];
        let b = fluid.particle_colour[(index * 3) + 2];

        draw_circle(1.0, 0.0, 0.0, pos_x, pos_y, particle_radius * 2.0, particle_radius * 2.0)
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
            max_particles,
        );

        create_particles(&mut scene.fluid, num_x, num_y);
        setup_grid(&mut scene.fluid);
        // TODO: finish this if we have time.
        // set_obstacle(&mut scene.fluid);

        return SimulationHandler { scene };
    }

    #[wasm_bindgen]
    pub fn update(&mut self, delta_time: f32)
    {
        // TODO: add pausing to scene
        if (true)
        {
            self.scene.fluid.simulate(
                delta_time,
                self.scene.gravity,
                self.scene.flip_ratio,
                self.scene.num_pressure_iters,
                self.scene.num_particle_iters,
                self.scene.over_relaxation,
                Some(self.scene.compensate_drift),
                self.scene.separate_particles,
                self.scene.obstacle_x,
                self.scene.obstacle_y,
                self.scene.obstacle_radius
            );
        }

        draw_simulation(&self.scene.fluid, 0.018/2.0);

        // hello();
        unsafe
        {
            CURRENT_INSTANCE = 0;
            SIM.integrate_balls(delta_time);
            SIM.handle_ball_wall_collision();
            //SIM.draw_balls();
        }
    }
}