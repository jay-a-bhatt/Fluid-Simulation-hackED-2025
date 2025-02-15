extern crate wasm_bindgen;
use wasm_bindgen::prelude::*;
const WASM_MEMORY_BUFFER_SIZE: usize = 8;
static mut WASM_MEMORY_BUFFER: [f32; WASM_MEMORY_BUFFER_SIZE] = [0.0; INSTANCE_DATA_SIZE]; 

/* Initialize the buffer with values manually
#[wasm_bindgen]
pub fn init_buffer()
{
    unsafe
    {
        // WASM_MEMORY_BUFFER[0] = 1.1;
    }
}
*/

const INSTANCE_DATA_SIZE: usize = 4+2+2;
const max_instances: u32 = 100;
static mut current_instance: usize = 0; // increments 0 - 799 = 800x before being reset

pub unsafe fn draw_circle(r: f32, g: f32, b: f32, a: f32, x: f32, y: f32, s_x: f32, s_y: f32)
{
    let mut array: [f32; 8] = [r, g, b, a, x, y, s_x, s_y];
    let instance_offset: usize = current_instance * 8; // instance offset corresponds to total # of floats in the buffer
    for i in 0.. 8
    {
        let index = instance_offset + i;
        WASM_MEMORY_BUFFER[index] = array[i];
    }
    if current_instance == 100
    {
        current_instance = 0;
    }
}

#[wasm_bindgen]
pub fn update()
{
    //
}

fn main()
{
    unsafe { draw_circle(50.5, 100.1, 120.5, 200.0, 123.2, 54.3, 40.12, 43.54);}
}