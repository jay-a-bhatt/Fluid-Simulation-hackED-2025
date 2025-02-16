extern crate wasm_bindgen;
use rand::{Rng, RngCore};
use std::path::Component::ParentDir;
use wasm_bindgen::prelude::*;
const WASM_MEMORY_BUFFER_SIZE: usize = 8;
static mut WASM_MEMORY_BUFFER: [f32; WASM_MEMORY_BUFFER_SIZE] = [0.0; INSTANCE_DATA_SIZE];

// BOUNCY BALL SIMULATOR!!!!

// GRAPHICS ---------------------------------------------------------------
const MAX_INSTANCES: usize = 100;
// NUMBER OF FLOATS PER INSTANCE
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
}
const INSTANCE_DATA_SIZE: usize = 4 + 2 + 2; //4+2+2
const max_instances: u32 = 100;
static mut current_instance: usize = 0; // increments 0 - 799 = 800x before being reset

pub unsafe fn draw_circle(r: f32, g: f32, b: f32, a: f32, x: f32, y: f32, s_x: f32, s_y: f32) {
    let mut array: [f32; 8] = [r, g, b, a, x, y, s_x, s_y];
    let instance_offset: usize = current_instance * 8; // instance offset corresponds to total # of floats in the buffer
    for i in 0..8 {
        let index = instance_offset + i;
        WASM_MEMORY_BUFFER[index] = array[i];
        println!("{}", i);
    }
    if current_instance == 100 {
        current_instance = 0;
    }
}

// provide js a pointer
#[wasm_bindgen]
pub fn return_pointer() -> *const f32 {
    let p: *const f32;
    unsafe {
        p = WASM_MEMORY_BUFFER.as_ptr();
    }

    return p;
}
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

static mut SIM: BounceSim = BounceSim { areaWidth:2.5f32, areaHeight:1.0f32, ballz: Vec::new() };

#[wasm_bindgen]
pub fn init_simulation()
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

// UPDATE FUNCTION
#[wasm_bindgen]
pub fn update(delta_time: f32)
{
    unsafe
    {
        CURRENT_INSTANCE = 0;
        SIM.integrate_balls(delta_time);
        SIM.handle_ball_wall_collision();
        SIM.draw_balls();
    }

    // Update simulation
    // Draw particles
}

struct context
{
    canvas_x: i32,
    canvas_y: i32,
    mouse_x: i32,
    mouse_y: i32
}

fn init_context(canvas_x: i32, canvas_y: i32, mouse_x: i32, mouse_y: i32) -> context
{
    context
    {
        canvas_x,
        canvas_y,
        mouse_x,
        mouse_y
    }

}
