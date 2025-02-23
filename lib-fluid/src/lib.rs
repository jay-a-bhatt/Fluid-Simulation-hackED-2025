
use std::{cell, cmp};
use libm::sqrtf;
use var::{FLUID_CELL, SOLID_CELL};
mod var;

fn clampF32(x: f32, min: f32, max: f32) -> f32 {
    if (x < min) { return min; }
    else if (x > max) { return max; }
    else { return x; }
}

pub struct Scene {
    pub gravity: f32,
    pub dt: f32,
    pub flip_ratio: f32,
    pub num_pressure_iters: i32,
    pub num_particle_iters: i32,
    frame_nr: i32,
    pub over_relaxation: f32,
    pub compensate_drift: bool,
    pub separate_particles: bool,
    pub obstacle_x: f32,
    pub obstacle_y: f32,
    pub obstacle_radius: f32,
    paused: bool,
    show_obstacle: bool,
    pub obstacle_vel_x: f32,
    pub obstacle_vel_y: f32,
    show_particles: bool,
    show_grid: bool,
    pub fluid: FlipFluid,
}

impl Scene {
    pub fn new(
        density: f32,
        width: f32,
        height: f32,
        cell_size: f32,
        particle_radius: f32,
        max_particles: i32,
    ) -> Self {
        Scene {
            gravity: -9.81,
            dt: 1.0 / 120.0,
            flip_ratio: 0.9,
            num_pressure_iters: 1,
            num_particle_iters: 2,
            frame_nr: 0,
            over_relaxation: 1.9,
            compensate_drift: true,
            separate_particles: true,
            obstacle_x: 1.0,
            obstacle_y: 1.5,
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
                cell_size,
                particle_radius,
                max_particles,
            ),
        }
    }
}
pub struct FlipFluid {
    density: f32,
    pub f_num_x: i32,
    pub f_num_y: i32,
    pub cell_size: f32,
    f_inv_spacing: f32,
    pub f_num_cells: i32,

    u: Vec<f32>,
    v: Vec<f32>,
    du: Vec<f32>,
    dv: Vec<f32>,
    prev_u: Vec<f32>,
    prev_v: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>,
    cell_type: Vec<i32>,
    pub cell_colour: Vec<f32>,

    max_particles: i32,

    pub particle_pos: Vec<f32>,
    pub particle_colour: Vec<f32>,
    particle_vel: Vec<f32>,
    particle_density: Vec<f32>,
    particle_rest_density: f32,

    pub particle_radius: f32,
    p_inv_spacing: f32,
    p_num_x: i32,
    p_num_y: i32,
    p_num_cells: i32,

    num_cell_particles: Vec<i32>,
    first_cell_particles: Vec<i32>,
    cell_particle_ids: Vec<i32>,

    pub num_particles: i32,
}

impl FlipFluid {
    fn new(
        density: f32,
        width: f32,
        height: f32,
        cell_size: f32,
        particle_radius: f32,
        max_particles: i32,
    ) -> Self {
        let f_num_x = (width / cell_size).floor() as i32 + 1;
        let f_num_y = (height / cell_size).floor() as i32 + 1;
        let h = (width / f_num_x as f32).max(height / f_num_y as f32);
        let f_inv_spacing = 1.0 / h;
        let f_num_cells = f_num_x * f_num_y;

        let p_inv_spacing = 1.0 / (2.2 * particle_radius);
        let p_num_x =  (width * p_inv_spacing).floor() as i32 + 1;
        let p_num_y = (height * p_inv_spacing).floor() as i32 + 1;
        let p_num_cells = p_num_x * p_num_y;

        FlipFluid {
            density,
            f_num_x,
            f_num_y,
            cell_size: h,
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

            num_cell_particles:   vec![0; p_num_cells as usize],
            first_cell_particles: vec![0; p_num_cells as usize + 1],
            cell_particle_ids:    vec![0; max_particles as usize],

            num_particles: 0,
        }
    }

    fn integrate_particles(&mut self, dt: f32, gravity: f32) {
        for i in 0..self.num_particles {
            self.particle_vel[(2 * i + 1) as usize] += dt * gravity;
            self.particle_pos[(2 * i + 0) as usize] += self.particle_vel[(2 * i + 0) as usize] * dt;
            self.particle_pos[(2 * i + 1) as usize] += self.particle_vel[(2 * i + 1) as usize] * dt;
        }
    }

    fn handle_particle_collisions(
        &mut self,
        obstacle_x: f32,
        obstacle_y: f32,
        obstacle_radius: f32,
        obstacle_vel_x: f32,
        obstacle_vel_y: f32
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

            if d2 < min_dist2
            {
                self.particle_vel[(2 * i) as usize] = obstacle_vel_x;
                self.particle_vel[(2 * i + 1) as usize] = obstacle_vel_y;
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

    fn push_particles_apart(&mut self, num_iters: i32) {
        let colour_diffusion_coeff: f32 = 0.001;

        // count particles per cell
        self.num_cell_particles.fill(0);

        for i in 0..self.num_particles
        {
            let x: f32 = self.particle_pos[(2 * i + 0) as usize];
            let y: f32 = self.particle_pos[(2 * i + 1) as usize];

            let xi: f32 = clampF32((x * self.p_inv_spacing).floor(), 0.0, (self.p_num_x - 1) as f32);
            let yi: f32 = clampF32((y * self.p_inv_spacing).floor(), 0.0, (self.p_num_y - 1) as f32);

            let cell_index = (xi * self.p_num_y as f32 + yi) as usize;
            self.num_cell_particles[cell_index] += 1;
        }

        //partial sums

        let mut first: i32 = 0;

        for i in 0..self.p_num_cells
        {
            first += self.num_cell_particles[i as usize];
            self.first_cell_particles[i as usize] = first;
        }
        self.first_cell_particles[self.p_num_cells as usize] = first; //guard

        // fill particles into cells
        for i in 0..self.num_particles
        {
            let x: f32 = self.particle_pos[(2 * i + 0) as usize];
            let y: f32 = self.particle_pos[(2 * i + 1) as usize];

            let xi: f32 = clampF32((x * self.p_inv_spacing).floor(), 0.0, (self.p_num_x - 1) as f32);
            let yi: f32 = clampF32((y * self.p_inv_spacing).floor(), 0.0, (self.p_num_y - 1) as f32);

            let cell_nr  = (xi * self.p_num_y as f32 + yi) as usize;
            self.first_cell_particles[cell_nr] -= 1;

            let fcp:usize = self.first_cell_particles[cell_nr] as usize;
            self.cell_particle_ids[fcp] = i;
        }

        // push particles apart

        let min_dist: f32 = 2.0 * self.particle_radius;
        let min_dist_2: f32 = min_dist * min_dist;

        for _ in 0..num_iters
        {
            for i in 0..self.num_particles
            {
                let px: f32 = self.particle_pos[(2 * i + 0) as usize];
                let py: f32 = self.particle_pos[(2 * i + 1) as usize];

                let pxi: i32 = (px * self.p_inv_spacing).floor() as i32;
                let pyi: i32 = (py * self.p_inv_spacing).floor() as i32;

                let x0: i32 = cmp::max(pxi - 1, 0);
                let y0: i32 = cmp::max(pyi - 1, 0);
                let x1: i32 = cmp::min(pxi + 1, self.p_num_x - 1);
                let y1: i32 = cmp::min(pyi + 1, self.p_num_y - 1);

                for xi in x0..=x1
                {
                    for yi in y0..=y1
                    {
                        let cell_nr: i32 = xi * self.p_num_y + yi;
                        let first: i32 = self.first_cell_particles[cell_nr as usize];
                        let last: i32  = self.first_cell_particles[(cell_nr + 1) as usize];

                        for j in first..last
                        {
                            let id: i32 = self.cell_particle_ids[j as usize];

                            if id == i { continue; }

                            let qx: f32 = self.particle_pos[(2 * id + 0) as usize];
                            let qy: f32 = self.particle_pos[(2 * id + 1) as usize];

                            let mut dx: f32 = qx - px;
                            let mut dy: f32 = qy - py;
                            let d2: f32 = (dx * dx) + (dy * dy);

                            if d2 > min_dist_2 || d2 == 0.0 {  continue;  }

                            let d: f32 = sqrtf(d2);
                            let s: f32 = 0.5 * (min_dist - d) / d;
                            dx *= s;
                            dy *= s;
                            self.particle_pos[(2 * i  + 0) as usize] -= dx;
                            self.particle_pos[(2 * i  + 1) as usize] -= dy;
                            self.particle_pos[(2 * id + 0) as usize] += dx;
                            self.particle_pos[(2 * id + 1) as usize] += dy;

                            // diffuse colours
                            /*
                            for k in 0..3
                            {
                                let colour0: f32 = self.particle_colour[(3 * i + k) as usize];
                                let colour1: f32 = self.particle_colour[(3 * id + k) as usize];
                                let colour: f32 = (colour0 + colour1) * 0.5;
                                self.particle_colour[(3 * i + k) as usize] =
                                    colour0 + (colour - colour0) * colour_diffusion_coeff;
                                self.particle_colour[(3 * id + k) as usize] =
                                    colour1 + (colour - colour1) * colour_diffusion_coeff;
                            }
                            */
                        }
                    }
                }
            }
        }
    }

    fn transfer_velocities(&mut self, to_grid: bool, flip_ratio: f32)
    {
        let n: f32 = self.f_num_y as f32;
        let h: f32 = self.cell_size;
        let h1: f32 = self.f_inv_spacing;
        let h2: f32 = h * 0.5;

        if to_grid
        {
            // Make a copy of the grid's previous velocities
            self.prev_u.copy_from_slice(&self.u);
            self.prev_v.copy_from_slice(&self.v);

            self.du.fill(0.0);
            self.u.fill(0.0);
            self.v.fill(0.0);
            self.dv.fill(0.0);

            for i in 0..self.f_num_cells
            {
                let ind = i as usize;
                if self.s[ind] == 0.0
                {
                    self.cell_type[ind] = var::SOLID_CELL;
                }
                else
                {
                    self.cell_type[ind] = var::AIR_CELL;
                }
            }

            for i in 0..self.num_particles
            {
                // X and Y position of particle.
                let x: f32 = self.particle_pos[(2 * i + 0) as usize];
                let y: f32 = self.particle_pos[(2 * i + 1) as usize];

                // X and Y index of particle
                let xi: f32 = clampF32((x * h1).floor(), 0.0, (self.f_num_x - 1) as f32);
                let yi: f32 = clampF32((y * h1).floor(), 0.0, (self.f_num_y - 1) as f32);

                let cell_nr = (xi * n + yi) as usize;

                // If the particle is inside an air cell, make it a fluid cell.
                if self.cell_type[cell_nr] == var::AIR_CELL
                {
                    self.cell_type[cell_nr] = var::FLUID_CELL;
                }
            }
        }

        // Transfer by x and y components 0 for x and 1 for y
        for component in 0..2
        {
            let dx: f32;
            let dy: f32;

            let mut f: &mut Vec<f32>;
            let mut prev_f: &mut Vec<f32>;
            let mut d: &mut Vec<f32>;

            if component == 0
            {
                dx = 0.0;
                dy = h2;
                // Someone left borrowing errors here, so imma just .clone() them
                // mfs for now and fix it later - Jay
                // It is now later - Rordon
                f = &mut self.u;
                prev_f = &mut self.prev_u;
                d = &mut self.du;
            }
            else
            {
                dx = h2;
                dy = 0.0;

                f = &mut self.v;
                prev_f = &mut self.prev_v;
                d = &mut self.dv;
            }

            for i in 0..self.num_particles
            {
                let mut x: f32 = self.particle_pos[(2 * i + 0) as usize];
                let mut y: f32 = self.particle_pos[(2 * i + 1) as usize];

                x = clampF32(x, h, (self.f_num_x as f32 - 1.0) * h);
                y = clampF32(y, h, (self.f_num_y as f32 - 1.0) * h);

                let x0: f32 = f32::min(((x - dx) * h1).floor(), self.f_num_x as f32 - 2.0);
                let tx: f32 = ((x - dx) - x0 * h) * h1;
                let x1: f32 = f32::min(x0 + 1.0, self.f_num_x as f32 - 2.0);

                let y0: f32 = f32::min(((y - dy) * h1).floor(), self.f_num_y as f32 - 2.0);
                let ty: f32 = ((y - dy) - y0 * h) * h1;
                let y1: f32 = f32::min(y0 + 1.0, self.f_num_y as f32 - 2.0);

                let sx: f32 = 1.0 - tx;
                let sy: f32 = 1.0 - ty;

                let d0: f32 = sx * sy;
                let d1: f32 = tx * sy;
                let d2: f32 = tx * ty;
                let d3: f32 = sx * ty;

                let nr0 =(x0 * n + y0) as usize;
                let nr1 =(x1 * n + y0) as usize;
                let nr2 =(x1 * n + y1) as usize;
                let nr3 =(x0 * n + y1) as usize;

                if to_grid
                {
                    let pv: f32 = self.particle_vel[(2 * i + component) as usize];
                    f[nr0] += pv * d0;
                    d[nr0] += d0;
                    f[nr1] += pv * d1;
                    d[nr1] += d1;
                    f[nr2] += pv * d2;
                    d[nr2] += d2;
                    f[nr3] += pv * d3;
                    d[nr3] += d3;
                }
                else
                {
                    let mut offset:usize = 1;
                    if component == 0 { offset = n as usize; }

                    let mut valid0: f32 = 0.0;
                    let mut valid1: f32 = 0.0;
                    let mut valid2: f32 = 0.0;
                    let mut valid3: f32 = 0.0;

                    if self.cell_type[nr0] != var::AIR_CELL || self.cell_type[nr0 - offset] != var::AIR_CELL
                    {
                        valid0 = 1.0;
                    }
                    if self.cell_type[nr1] != var::AIR_CELL || self.cell_type[nr1 - offset] != var::AIR_CELL
                    {
                        valid1 = 1.0;
                    }
                    if self.cell_type[nr2] != var::AIR_CELL || self.cell_type[nr2 - offset] != var::AIR_CELL
                    {
                        valid2 = 1.0;
                    }
                    if self.cell_type[nr3] != var::AIR_CELL || self.cell_type[nr3 - offset] != var::AIR_CELL
                    {
                        valid3 = 1.0;
                    }

                    let v: f32 = self.particle_vel[(2 * i + component) as usize];
                    let d: f32 = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d > 0.0
                    {
                        let pic_v: f32 = (
                              valid0 * d0 * f[nr0]
                            + valid1 * d1 * f[nr1]
                            + valid2 * d2 * f[nr2]
                            + valid3 * d3 * f[nr3]
                        ) / d;
                        let corr: f32 = (
                              valid0 * d0 * (f[nr0] - prev_f[nr0])
                            + valid1 * d1 * (f[nr1] - prev_f[nr1])
                            + valid2 * d2 * (f[nr2] - prev_f[nr2])
                            + valid3 * d3 * (f[nr3] - prev_f[nr3])
                        ) / d;

                        let flip_v: f32 = v + corr;
                        self.particle_vel[(2 * i + component) as usize] = (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v;
                    }
                }
            }

            if to_grid
            {
                // No idea what this is..,
                for i in 0..f.len()
                {
                    if d[i] > 0.0
                    {
                        f[i] /= d[i];
                    }
                }

                // Restore solid cells whatever that means
                for i in 0..self.f_num_x
                {
                    for j in 0..self.f_num_y
                    {
                        let ny = self.f_num_y; // Number of cells on y axis.
                        let center= (i * ny + j) as usize; // Index of center cell
                        let left  = ((i - 1) * ny + j) as usize;
                        let bottom= (i * ny + (j - 1)) as usize;

                        let isSolidCell = self.cell_type[center] == SOLID_CELL;

                        if isSolidCell || (i > 0 && self.cell_type[left] == SOLID_CELL)
                        {
                            self.u[center] = self.prev_u[center];
                        }

                        if isSolidCell || (j > 0 && self.cell_type[bottom] == SOLID_CELL)
                        {
                            self.v[center] = self.prev_v[center];
                        }
                    }
                }
            }
        }
    }

    fn update_particle_density(&mut self)
    {
        let n: usize = self.f_num_y as usize;

        // Cell Size
        let h  = self.cell_size;
        let h1 = self.f_inv_spacing;
        // 1/2 of Cell Size
        let h2 = self.cell_size * 0.5;


        let nx = self.f_num_x as f32; // Number of cells on x as float32
        let ny = self.f_num_y as f32; // Number of cells on y as float32

        self.particle_density.fill(0.0);

        for i in 0..self.num_particles
        {
            // Particle X and Y positions.
            let mut px = self.particle_pos[(2 * i) as usize];
            let mut py = self.particle_pos[(2 * i + 1) as usize];

            // Clamp position on each axis between Cell Size and Cell Size * Number of Cells
            px = clampF32(px, h, (nx - 1.0) * h);
            py = clampF32(py, h, (ny - 1.0) * h);

            let x0 = ((px - h2) * h1).floor();
            let tx = ((px - h2) - x0 * h) * h1;
            let x1 = f32::min(x0 + 1.0, nx - 2.0);

            let y0 = ((py - h2) * h1).floor();
            let ty = ((py - h2) - y0 * h) * h1;
            let y1 = f32::min(y0 + 1.0, ny - 2.0);

            let sx = 1.0 - tx;
            let sy = 1.0 - ty;

            if x0 < nx && y0 < ny { self.particle_density[(x0 as usize) * n + (y0 as usize)] += sx * sy; }
            if x1 < nx && y0 < ny { self.particle_density[(x1 as usize) * n + (y0 as usize)] += tx * sy; }
            if x1 < nx && y1 < ny { self.particle_density[(x1 as usize) * n + (y1 as usize)] += tx * ty; }
            if x0 < nx && y1 < ny { self.particle_density[(x0 as usize) * n + (y1 as usize)] += sx * ty; }
        }

        if self.particle_rest_density == 0.0
        {
            let mut sum = 0.0;
            let num_fluid_cells = 0;

            for i in 0..self.f_num_cells
            {
                if self.cell_type[i as usize] == FLUID_CELL
                {
                    sum += self.particle_density[i as usize];
                }
            }

            if num_fluid_cells > 0
            {
                self.particle_rest_density = sum / num_fluid_cells as f32;
            }
        }
    }
    fn solve_incompressibility(&mut self, num_iters: i32, dt: f32, over_relaxation: f32, compensate_drift: bool)
    {
        self.p.fill(0.0);

        self.prev_u.copy_from_slice(&self.u);
        self.prev_v.copy_from_slice(&self.v);

        let n:  i32 = self.f_num_y;
        let cp: f32 = (self.density * self.cell_size) / dt;

        /* NOTE(rordon): these values are never used in reference sim
        for i in 0..self.f_num_cells
        {
            let u: f32 = self.u[i as usize];
            let v: f32 = self.v[i as usize];
        }
        */

        for _ in 0..num_iters
        {
            for i in 1..self.f_num_x - 1
            {
                for j in 1..self.f_num_y - 1
                {
                    if self.cell_type[((i * n) + j) as usize] != FLUID_CELL
                    {
                        continue;
                    }

                    // Indices
                    let center= ((i * n) + j) as usize;
                    let left  = ((i - 1) * n + j) as usize;
                    let right = ((i + 1) * n + j) as usize;
                    let bottom= (i * n + j - 1) as usize;
                    let top   = (i * n + j + 1) as usize;

                    let s:   f32 = self.s[center];
                    let sx0: f32 = self.s[left];
                    let sx1: f32 = self.s[right];
                    let sy0: f32 = self.s[bottom];
                    let sy1: f32 = self.s[top];

                    let s: f32 = sx0 + sx1 + sy0 + sy1;
                    if s == 0.0 { continue; }

                    let mut div: f32 = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                    if self.particle_rest_density > 0.0 && compensate_drift
                    {
                        let k: f32 = 1.0;
                        let compression = self.particle_density[center] - self.particle_rest_density;
                        if compression > 0.0 { div = div - k * compression; }
                    }

                    let mut p: f32 = -div / s;
                    p *= over_relaxation;
                    self.p[center] += cp * p;

                    self.u[center] -= sx0 * p;
                    self.u[right]  += sx1 * p;
                    self.v[center] -= sy0 * p;
                    self.v[top]    += sy1 * p;
                }
            }
        }
    }

    fn update_particle_colours(&mut self) {
        let h1: f32 = self.f_inv_spacing;

        for i in 0..self.num_particles {
            let s: f32 = 0.01;

            self.particle_colour[(3 * i) as usize] =
                (clampF32(self.particle_colour[(3 * i) as usize] - s, 0.0, 1.0)) as f32;
            self.particle_colour[(3 * i + 1) as usize] =
                (clampF32(self.particle_colour[(3 * i + 1) as usize] - s, 0.0, 1.0)) as f32;
            self.particle_colour[(3 * i + 2) as usize] =
                (clampF32(self.particle_colour[(3 * i + 2) as usize] - s, 0.0, 1.0)) as f32;

            let x: f32 = self.particle_pos[(2 * i) as usize];
            let y: f32 = self.particle_pos[(2 * i + 1) as usize];
            let xi: f32 = clampF32((x * h1).floor(), 1.0, (self.f_num_x - 1) as f32);
            let yi: f32 = clampF32((y * h1).floor(), 1.0, (self.f_num_y - 1) as f32);
            let cell_nr: f32 = xi * self.f_num_y as f32 + yi;

            let d0: f32 = self.particle_rest_density;

            if d0 > 0.0 {
                let rel_density: f32 = self.particle_density[(cell_nr) as usize] / d0;
                if rel_density < 0.7 {
                    let s: f32 = 0.8;
                    self.particle_colour[(3 * i) as usize] = s;
                    self.particle_colour[(3 * i + 1) as usize] = s;
                    self.particle_colour[(3 * i + 2) as usize] = 1.0;
                }
            }
        }
    }

    fn update_cell_colours(&mut self) {
        self.cell_colour.fill(1.0);
        for i in 0..self.f_num_cells
        {
            if self.cell_type[i as usize] == SOLID_CELL
            {
                let index: usize = 3 * i as usize;
                self.cell_colour[index + 0] = 0.5;
                self.cell_colour[index + 1] = 0.5;
                self.cell_colour[index + 2] = 0.5;
            }
            else if self.cell_type[i as usize] == FLUID_CELL
            {
                let mut d = self.particle_density[i as usize];
                if self.particle_rest_density > 0.0
                {
                    d /= self.particle_rest_density;
                }
                self.set_sci_colour(i, d, 0.0, 2.0);
            }
        }
    }

    fn set_sci_colour(&mut self, cell_nr: i32, mut val: f32, min_val: f32, max_val: f32)
    {
        val = f32::min(f32::max(val, min_val), max_val - 0.0001);
        let d: f32 = max_val - min_val;

        if d == 0.0
        {
            val = 0.5;
        }
        else
        {
            val = (val - min_val) / d;
        }

        let m: f32 = 0.25;
        let num: f32 = (val / m).floor();
        let s: f32 = (val - num * m) / m;
        let r: f32;
        let g: f32;
        let b: f32;

        match num as i32 {
            0 => { r = 0.0; g = s;   b = 1.0; }

            1 => { r = 0.0; g = 1.0; b = 1.0 - s; }

            2 => { r = s;   g = 1.0; b = 0.0 }
            3 => { r = 1.0; g = 1.0 - s; b = 0.0 }
            _ => { r = 0.0; g = 0.0; b = 0.0 }
        }

        self.cell_colour[(3 * cell_nr + 0) as usize] = r;
        self.cell_colour[(3 * cell_nr + 1) as usize] = g;
        self.cell_colour[(3 * cell_nr + 2) as usize] = b;
    }

    pub fn simulate(
        &mut self,
        dt: f32,
        gravity: f32,
        flip_ratio: f32,
        num_pressure_iters: i32,
        num_particle_iters: i32,
        over_relaxation: f32,
        compensate_drift: bool,
        separate_particles: bool,
        obstacle_x: f32,
        obstacle_y: f32,
        obstacle_vel_x: f32,
        obstacle_vel_y: f32,
        obstacle_radius: f32,
    ) {
        let num_sub_steps = 1;
        let sdt = dt / num_sub_steps as f32;

        for _ in 0..num_sub_steps {
            self.integrate_particles(sdt, gravity); // NOTE(rordon): THIS IS GOOD!
            if separate_particles { self.push_particles_apart(num_particle_iters); }

            self.handle_particle_collisions(obstacle_x, obstacle_y, obstacle_radius, obstacle_vel_x, obstacle_vel_y);

            self.transfer_velocities(true, 0.0);
            self.update_particle_density();
            self.solve_incompressibility(num_pressure_iters, sdt, over_relaxation, compensate_drift);
            self.transfer_velocities(false, 0.9);
        }

        //self.update_particle_colours();
        self.update_cell_colours();
    }
}

pub fn create_particles(fluid: &mut FlipFluid, num_x: i32, num_y: i32)
{
    fluid.num_particles = num_x * num_y;

    let dist_x: f32 = 2.0 * fluid.particle_radius;
    let dist_y: f32 = sqrtf(3.0) / 2.0 * dist_x;

    let cell_size = fluid.cell_size;
    let particle_rad = fluid.particle_radius;
    let mut p:usize = 0;
    for i in 0..num_x
    {
        for j in 0..num_y
        {
            let mut offset_x: f32 = particle_rad;
            if (j % 2 == 0) { offset_x = 0.0 }

            fluid.particle_pos[p] = cell_size + particle_rad + dist_x * (i as f32) + offset_x;
            p += 1;
            fluid.particle_pos[p] = cell_size + particle_rad + dist_y * (j as f32);
            p += 1;
        }
    }
}

pub fn setup_grid(fluid: &mut FlipFluid)
{
    for x in 0..fluid.f_num_x // Columns
    {
        for y in 0..fluid.f_num_y // Rows
        {
            let mut state = 1.0; // Default fluid state
            if (x == 0 || x == fluid.f_num_x - 1 || y == 0) { state = 0.0; }

            let cell_index = (x * fluid.f_num_y + y) as usize;
            fluid.s[cell_index] = state;

            let r = x as f32 / fluid.f_num_x as f32;
            let g = y as f32 / fluid.f_num_y as f32;

            fluid.cell_colour[cell_index * 3 + 0] = r;
            fluid.cell_colour[cell_index * 3 + 1] = g;
            fluid.cell_colour[cell_index * 3 + 2] = 0.0;
        }
    }
}

pub fn set_obstacle(scene: &mut Scene, mouse_x: f32, mouse_y: f32, reset: bool)
{
    let mut vel_x: f32 = 0.0;
    let mut vel_y: f32 = 0.0;

    if !reset
    {
        vel_x = (mouse_x - scene.obstacle_x) / scene.dt;
        vel_y = (mouse_y - scene.obstacle_y) / scene.dt;
    }

    scene.obstacle_x = mouse_x;
    scene.obstacle_y = mouse_y;
    let obs_rad: f32 = scene.obstacle_radius;
    let mut fluid = &mut scene.fluid;

    for i in 1..fluid.f_num_x-2
    {
        for j in 1..fluid.f_num_y-2
        {
            // Not confusing at all...
            let index_1 = (i * fluid.f_num_y + j) as usize;
            let index_2 = ((i+1) * fluid.f_num_y + j) as usize;
            let index_3 = (i * fluid.f_num_y + (j+1)) as usize;

            let dist_x = (i as f32 + 0.5) * fluid.cell_size - mouse_x;
            let dist_y = (j as f32 + 0.5) * fluid.cell_size - mouse_y;

            fluid.s[index_1] = 1.0;

            if (dist_x * dist_x + dist_y * dist_y <= (obs_rad * obs_rad))
            {
                fluid.s[index_1] = 0.0;

                fluid.u[index_1] = vel_x;
                fluid.u[index_2] = vel_x;

                fluid.v[index_1] = vel_y;
                fluid.v[index_3] = vel_y;
            }
        }
    }

    scene.show_obstacle = true;
    scene.obstacle_vel_x = vel_x;
    scene.obstacle_vel_y = vel_y;
}