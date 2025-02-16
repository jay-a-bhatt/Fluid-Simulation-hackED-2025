use std::{cell, cmp};

fn clamp(x: f32, min: f32, max: f32) -> i32 {
    if (x < min) {
        return (min) as i32;
    } else if (x > max) {
        return (max) as i32;
    } else {
        return (x) as i32;
    }
}

struct Scene {
    gravity: f32,
    dt: f32,
    flip_ratio: f32,
    num_pressure_iters: i32,
    num_particle_iters: i32,
    frame_nr: i32,
    over_relaxation: f32,
    compensate_drift: bool,
    separate_particles: bool,
    obstacle_x: f32,
    obstacle_y: f32,
    obstacle_radius: f32,
    paused: bool,
    show_obstacle: bool,
    obstacle_vel_x: f32,
    obstacle_vel_y: f32,
    show_particles: bool,
    show_grid: bool,
    fluid: FlipFluid,
}

impl Scene {
    fn new(
        self,
        density: f32,
        width: f32,
        height: f32,
        spacing: f32,
        particle_radius: f32,
        max_particles: i32,
    ) -> Self {
        Scene {
            gravity: -9.81,
            dt: 1.0 / 120.0,
            flip_ratio: 0.9,
            num_pressure_iters: 100,
            num_particle_iters: 2,
            frame_nr: 0,
            over_relaxation: 1.9,
            compensate_drift: true,
            separate_particles: true,
            obstacle_x: 0.0,
            obstacle_y: 0.0,
            obstacle_radius: 0.15,
            paused: true,
            show_obstacle: true,
            obstacle_vel_x: 0.0,
            obstacle_vel_y: 0.0,
            show_particles: true,
            show_grid: false,
            fluid: FlipFluid::new(
                density,
                width,
                height,
                spacing,
                particle_radius,
                max_particles,
            ),
        }
    }
}
struct FlipFluid {
    density: f32,
    f_num_x: i32,
    f_num_y: i32,
    h: f32,
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
    cell_type: Vec<i32>,
    cell_colour: Vec<f32>,

    max_particles: i32,

    particle_pos: Vec<f32>,
    particle_colour: Vec<f32>,
    particle_vel: Vec<f32>,
    particle_density: Vec<f32>,
    particle_rest_density: f32,

    particle_radius: f32,
    p_inv_spacing: f32,
    p_num_x: i32,
    p_num_y: i32,
    p_num_cells: i32,

    num_cell_particles: Vec<i32>,
    first_cell_particles: Vec<i32>,
    cell_particle_ids: Vec<i32>,

    num_particles: i32,
}

impl FlipFluid {
    fn new(
        density: f32,
        width: f32,
        height: f32,
        spacing: f32,
        particle_radius: f32,
        max_particles: i32,
    ) -> Self {
        let f_num_x = (width / spacing).floor() as i32 + 1;
        let f_num_y = (height / spacing).floor() as i32 + 1;
        let h = (width / f_num_x as f32).max(height / f_num_y as f32);
        let f_inv_spacing = 1.0 / h;
        let f_num_cells = f_num_x * f_num_y;

        let p_inv_spacing = 1.0 / (2.2 * particle_radius);
        let p_num_x = (width * p_inv_spacing).floor() as i32 + 1;
        let p_num_y = (height * p_inv_spacing).floor() as i32 + 1;
        let p_num_cells = p_num_x * p_num_y;

        FlipFluid {
            density,
            f_num_x,
            f_num_y,
            h,
            f_inv_spacing,
            f_num_cells,

            u: vec![0.0; f_num_cells as usize],
            v: vec![0.0; f_num_cells as usize],
            du: vec![0.0; f_num_cells as usize],
            dv: vec![0.0; f_num_cells as usize],
            prev_u: vec![0.0; f_num_cells as usize],
            prev_v: vec![0.0; f_num_cells as usize],
            p: vec![0.0; f_num_cells as usize],
            s: vec![0.0; f_num_cells as usize],
            cell_type: vec![0; f_num_cells as usize],
            cell_colour: vec![0.0; 3 * f_num_cells as usize],

            max_particles,

            particle_pos: vec![0.0; 2 * max_particles as usize],
            particle_colour: vec![0.0; 3 * max_particles as usize],
            particle_vel: vec![0.0; 2 * max_particles as usize],
            particle_density: vec![0.0; f_num_cells as usize],
            particle_rest_density: 0.0,

            particle_radius,
            p_inv_spacing,
            p_num_x,
            p_num_y,
            p_num_cells,

            num_cell_particles: vec![0; p_num_cells as usize],
            first_cell_particles: vec![0; p_num_cells as usize + 1],
            cell_particle_ids: vec![0; max_particles as usize],

            num_particles: 0,
        }
    }

    fn integrate_particles(mut self, dt: f32, gravity: f32) {
        for i in 0..self.num_particles {
            self.particle_vel[(2 * i + 1) as usize] += dt * gravity;
            self.particle_pos[(2 * i) as usize] += self.particle_vel[(2 * 1) as usize] * dt;
            self.particle_pos[(2 * i + 1) as usize] += self.particle_vel[(2 * i + 1) as usize] * dt;
        }
    }

    fn handle_particle_collisions(
        mut self,
        obstacle_x: f32,
        obstacle_y: f32,
        obstacle_radius: f32,
    ) {
        let h: f32 = 1.0 / self.f_inv_spacing;
        let r: f32 = self.particle_radius as f32;
        let or2: f32 = obstacle_radius * obstacle_radius;
        let min_dist: f32 = obstacle_radius + r;
        let min_dist2: f32 = min_dist * min_dist;
        let min_x: f32 = h + r;
        let min_y: f32 = h + r;
        let max_x: f32 = (self.f_num_x - 1) as f32 * h - r;
        let max_y: f32 = (self.f_num_y - 1) as f32 * h - r;

        for i in 0..self.num_particles {
            let mut x = self.particle_pos[(2 * i) as usize];
            let mut y: f32 = self.particle_pos[(2 * i + 1) as usize];
            let dx: f32 = x - obstacle_x;
            let dy: f32 = y - obstacle_y;
            let d2: f32 = dx * dx + dy * dy;

            //obstacle collision

            if d2 < min_dist2 {
                self.particle_vel[(2 * i) as usize] = 0.0; //scene.obstacle_vel_x
                self.particle_vel[(2 * i + 1) as usize] = 0.0; //scene.obstacle_vel_y
            }

            //wall collisions

            if x < min_x {
                x = min_x;
                self.particle_vel[(2 * i) as usize] = 0.0;
            }

            if x > max_x {
                x = max_x;
                self.particle_vel[(2 * i) as usize] = 0.0;
            }

            if y < min_y {
                y = min_y;
                self.particle_vel[(2 * i + 1) as usize] = 0.0;
            }

            if y > max_y {
                y = max_y;
                self.particle_vel[(2 * i + 1) as usize] = 0.0;
            }
            self.particle_pos[(2 * i) as usize] = x;
            self.particle_pos[(2 * i + 1) as usize] = y;
        }
    }

    fn push_particles_apart(mut self, num_iters: i32) {
        let colour_diffusion_coeff: f32 = 0.001;

        // count particles per cell

        self.num_cell_particles.fill(0);

        for i in 0..self.num_particles {
            let x: f32 = self.particle_pos[(2 * i) as usize];
            let y: f32 = self.particle_pos[(2 * i + 1) as usize];
            let xi: i32 = clamp(
                (x * self.p_inv_spacing).floor(),
                0.0,
                (self.p_num_x - 1) as f32,
            );
            let yi: i32 = clamp(
                (y * self.p_inv_spacing).floor(),
                0.0,
                (self.p_num_y - 1) as f32,
            );
            let cell_nr: i32 = xi * self.p_num_y + yi;
            self.num_cell_particles[(cell_nr) as usize] += 1;
        }

        //partial sums

        let mut first: i32 = 0;

        for i in 0..self.p_num_cells {
            first += self.num_cell_particles[(i) as usize];
            self.first_cell_particles[(i) as usize] = first;
        }
        self.first_cell_particles[(self.p_num_cells) as usize] = first; //guard

        // fill particles into cells

        for i in 0..self.num_particles {
            let x: f32 = self.particle_pos[(2 * i) as usize];
            let y: f32 = self.particle_pos[(2 * i + 1) as usize];

            let xi: i32 = clamp(
                (x * self.p_inv_spacing).floor(),
                0.0,
                (self.p_num_x - 1) as f32,
            );
            let yi: i32 = clamp(
                (y * self.p_inv_spacing).floor(),
                0.0,
                (self.p_num_y - 1) as f32,
            );
            let cell_nr: i32 = xi * self.p_num_y + yi;
            self.num_cell_particles[(cell_nr) as usize] -= 1;
            self.cell_particle_ids[(self.first_cell_particles[(cell_nr) as usize]) as usize] = i;
        }

        //push particles apart

        let min_dist: f32 = 2.0 * self.particle_radius;
        let min_dist_2: f32 = min_dist * min_dist;

        for iter in 0..num_iters {
            for i in 0..self.num_particles {
                let px: f32 = self.particle_pos[(2 * i) as usize];
                let py: f32 = self.particle_pos[(2 * i + 1) as usize];

                let pxi: i32 = ((px * self.p_inv_spacing).floor()) as i32;
                let pyi: i32 = ((py * self.p_inv_spacing).floor()) as i32;
                let x0: i32 = cmp::max(pxi - 1, 0);
                let y0: i32 = cmp::max(pyi - 1, 0);
                let x1: i32 = cmp::min(pxi + 1, self.p_num_x - 1);
                let y1: i32 = cmp::min(pyi + 1, self.p_num_y - 1);

                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cell_nr: i32 = xi * self.p_num_y + yi;
                        let first: i32 = self.first_cell_particles[(cell_nr) as usize];
                        let last: i32 = self.first_cell_particles[(cell_nr + 1) as usize];
                        for j in first..last {
                            let id: i32 = self.cell_particle_ids[j as usize];
                            if id == i {
                                continue;
                            }
                            let qx: f32 = self.particle_pos[(2 * id) as usize];
                            let qy: f32 = self.particle_pos[(2 * id + 1) as usize];

                            let mut dx: f32 = qx - px;
                            let mut dy: f32 = qy - py;
                            let d2: f32 = dx * dx + dy * dy;
                            if d2 > min_dist_2 || d2 == 0.0 {
                                continue;
                            }
                            let d: f32 = libm::sqrtf(d2);
                            let s: f32 = 0.5 * (min_dist - d) / d;
                            dx *= s;
                            dy *= s;
                            self.particle_pos[(2 * i) as usize] -= dx;
                            self.particle_pos[(2 * i + 1) as usize] -= dy;
                            self.particle_pos[(2 * id) as usize] += dx;
                            self.particle_pos[(2 * id + 1) as usize] -= dy;

                            // diffuse colours

                            for k in 0..3 {
                                let colour0: f32 = self.particle_colour[(3 * i + k) as usize];
                                let colour1: f32 = self.particle_colour[(3 * id + k) as usize];
                                let colour: f32 = (colour0 + colour1) * 0.5;
                                self.particle_colour[(3 * i + k) as usize] =
                                    colour0 + (colour - colour0) * colour_diffusion_coeff;
                                self.particle_colour[(3 * id + k) as usize] =
                                    colour1 + (colour - colour1) * colour_diffusion_coeff;
                            }
                        }
                    }
                }
            }
        }
    }

    fn transfer_velocities(mut self, to_grid: bool, flip_ratio: f32) {
        let n: f32 = self.f_num_y as f32;
        let h: f32 = self.h;
        let h1: f32 = self.f_inv_spacing;
        let h2: f32 = h * 0.5;

        if (to_grid) {
            self.prev_u = self.u.clone(); //Using .clone() so might use more memory
            self.prev_v = self.v.clone();

            self.du.fill(0.0);
            self.dv.fill(0.0);
            self.u.fill(0.0);
            self.v.fill(0.0);

            for i in 0..self.f_num_cells {
                if self.s[i as usize] == 0.0 {
                    self.cell_type[i as usize] = var::SOLID_CELL;
                } else {
                    self.cell_type[i as usize] = var::AIR_CELL;
                }
            }

            for i in 0..self.num_particles {
                let x: f32 = self.particle_pos[(2 * i) as usize];
                let y: f32 = self.particle_pos[(2 * i + 1) as usize];
                let xi: i32 = clamp((x * h1).floor(), 0.0, (self.f_num_x - 1) as f32);
                let yi: i32 = clamp((y * h1).floor(), 0.0, (self.f_num_y - 1) as f32);
                let cell_nr: i32 = xi * n as i32+ yi;

                if self.cell_type[cell_nr as usize] == var::AIR_CELL {
                    self.cell_type[cell_nr as usize] = var::FLUID_CELL;
                }
            }
        }

        for component in 0..2 {
            let mut dx: f32;
            let mut dy: f32;

            let mut f: Vec<f32>;
            let mut prev_f: Vec<f32>;
            let mut d: Vec<f32>; 

            if component == 0 {
                dx = 0.0;
                dy = h2;

                f = self.u;
                prev_f = self.prev_u;
                d = self.du;
            } else {
                dx = h2;
                dy = 0.0;

                f = self.v;
                prev_f = self.prev_v;
                d = self.dv;
            }

            for i in 0..self.num_particles{
                let mut x: f32 = self.particle_pos[(2*i) as usize];
                let mut y: f32 = self.particle_pos[(2*i+1) as usize];

                x = clamp(x,h,(self.f_num_x as f32 - 1.0) * h) as f32;
                y = clamp(y,h,(self.f_num_y as f32 - 1.0) * h) as f32;

                let x0: f32 = f32::min(((x-dx)*h1).floor(),self.f_num_x as f32-2.0);
                let tx: f32 = ((x-dx)- x0 * h) * h1;
                let x1: f32 = f32::min((x0+1.0), self.f_num_x as f32 -2.0);

                let y0: f32 = f32::min(((y-dy)*h1).floor(), self.f_num_y as f32-2.0);
                let ty: f32 = ((y - dy)- y0*h)*h1;
                let y1: f32 = f32::min(y0+1.0, self.f_num_y as f32 - 2.0);

                let sx: f32 = 1.0 - tx;
                let sy: f32 = 1.0 - ty;

                let d0: f32 = sx*sy;
                let d1: f32 = tx*sy;
                let d2: f32 = tx*ty;
                let d3: f32 = sx*ty;

                let nr0: f32 = x0*n+y0;
                let nr1: f32 = x1*n+y0;
                let nr2: f32 = x1*n+y1;
                let nr3: f32 = x0*n+y1;

                if to_grid{
                    let pv: f32 = self.particle_vel[(2*i+component) as usize];
                    f[nr0 as usize] += pv * d0; d[nr0 as usize] += d0;
                    f[nr1 as usize] += pv * d1; d[nr1 as usize] += d1;
                    f[nr2 as usize] += pv * d2; d[nr2 as usize] += d2;
                    f[nr3 as usize] += pv * d3; d[nr3 as usize] += d3;
                }else{
                    let mut offset: f32 = 1.0;
                    let mut valid0: f32 = 0.0;
                    let mut valid1: f32 = 0.0;
                    let mut valid2: f32 = 0.0;
                    let mut valid3: f32 = 0.0;
                    if component == 0{
                        offset = n;
                    }
                    if (self.cell_type[nr0 as usize] != var::AIR_CELL || self.cell_type[nr0 as usize - offset as usize] != var::AIR_CELL){
                        valid0 = 1.0;
                    }
                    if (self.cell_type[nr1 as usize] != var::AIR_CELL || self.cell_type[nr1 as usize - offset as usize] != var::AIR_CELL){
                        valid1 = 1.0;
                    }
                    if (self.cell_type[nr2 as usize] != var::AIR_CELL || self.cell_type[nr2 as usize - offset as usize] != var::AIR_CELL){
                        valid2 = 1.0;
                    }
                    if (self.cell_type[nr3 as usize] != var::AIR_CELL || self.cell_type[nr3 as usize - offset as usize] != var::AIR_CELL){
                        valid3 = 1.0;
                    }

                    let v: f32 = self.particle_vel[(2*i+component) as usize];
                    let d: f32 = valid0 * d0 + valid1 *d1 + valid2 * d2 + valid3 * d3;

                    if d > 0.0 {
                        let pic_v: f32 = (valid0*d0*f[nr0 as usize]+valid1*d1*f[nr1 as usize]+valid2*d2*f[nr2 as usize]+valid3*d3*f[nr3 as usize])/d;
                        let corr: f32 = (valid0*d0*(f[nr0 as usize])-prev_f[nr0 as usize]+valid1*d1*(f[nr1 as usize])-prev_f[nr1 as usize]+valid2*d2*(f[nr2 as usize])-prev_f[nr2 as usize]+valid3*d3*(f[nr3 as usize])-prev_f[nr3 as usize])/d;
                        let flip_v: f32 = v + corr;

                        self.particle_vel[(2*i+component) as usize] = (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v;
                    }
                }
            }
            if to_grid {
                for i in 0..f.len(){
                    if d[i as usize] > 0.0{
                        f[i as usize] /= d[i as usize];
                    }
                }
            
                // restore solid cells
            
                for i in 0..self.f_num_x{
                    for j in 0..self.f_num_y{
                        let solid = self.cell_type[(i as f32 * n + j as f32) as usize] == var::SOLID_CELL;
                        if solid || (i > 0 && self.cell_type[((i - 1) as f32 * n + j as f32) as usize] == var::SOLID_CELL){
                            self.u[(i as f32 * n + j as f32) as usize] = self.prev_u[(i as f32 * n + j as f32) as usize];
                        }
                        if solid || (j > 0 && self.cell_type[(i as f32 * n + j as f32 - 1.0) as usize] == var::SOLID_CELL){
                            self.v[(i as f32 * n + j as f32) as usize] = self.prev_v[(i as f32 * n + j as f32) as usize];
                        }
                    }
                }
            }
        }
    }
}
fn main() {
    println!("Hello, world!");
}
