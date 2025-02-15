struct FlipFluid {
    density: i32,
    f_num_x: f32,
    f_num_y: f32,
    h: i32,
    f_inv_spacing: f32,
    f_num_cells: i32,

    u: Vec<f32>,
    v: Vec<f32>,
    du: Vec<f32>,
    dv: Vec<f32>,
    prev_u: Vec<f32>,
    prev_v: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>,
    cell_type: Vec<f32>,
    cell_colour: Vec<f32>,

    max_particles: i32,

    particle_pos: Vec<f32>,
    particle_colour: Vec<f32>,

    particle_vel: Vec<f32>,
    particle_density: Vec<f32>,
    particle_rest_density: Vec<f32>,

    particle_radius: i32,
    p_inv_spacing: f32,
    p_num_x: f32,
    p_num_y: f32,
    p_num_cells: f32,

    num_cell_particles: Vec<i32>,
    first_cell_particles: Vec<i32>,
    cell_particle_ids: Vec<i32>,

    num_particles: i32,
}

impl FlipFluid {
    fn new(
        density: i32,
        width: i32,
        height: i32,
        spacing: i32,
        particle_radius: i32,
        max_particles: i32,
    ) {
    }
}

fn main() {
    println!("Hello, world!");
}
fn integrate_particles(dt: f32 , gravity: f32){
    let num_particles = 0;
    for i in 0..num_particles{

    }
    //Replace arrays with particles
    todo!()
}