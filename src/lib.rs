extern crate wasm_bindgen;
use wasm_bindgen::prelude::*;
//
const WASM_MEMORY_BUFFER_SIZE: usize = 4;
static mut WASM_MEMORY_BUFFER: [f32; WASM_MEMORY_BUFFER_SIZE] = [0.0; WASM_MEMORY_BUFFER_SIZE];

#[wasm_bindgen]
pub fn init_buffer() {
    unsafe {
        // Initialize the buffer with values
        WASM_MEMORY_BUFFER[0] = 1.1;
        WASM_MEMORY_BUFFER[1] = 2.2;
        WASM_MEMORY_BUFFER[2] = 3.3;
        WASM_MEMORY_BUFFER[3] = 4.4;
    }
}

#[wasm_bindgen]
pub fn buffer_pointer() -> *const f32 {
    let p: *const f32;
    unsafe {
        p = WASM_MEMORY_BUFFER.as_ptr();
    }
    return p;
}

fn context()
{
    println!("Hello, world!");
}




// destructure the tuple into two variables names

