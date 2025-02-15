use std::cmp;

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

    fn update_particle_density(self) {
        let n = self.f_num_y;
        let h = self.h;
        let h1 = self.f_inv_spacing;
        let h2 = self.h * 0.5;

        let mut d = self.particle_density;
        d.fill(0.0);

        for i in 0..self.num_particles {
            let mut x = self.particle_pos[(2 * i) as usize];
            let mut y = self.particle_pos[(2 * i + 1) as usize];

            x = clamp(x, h, (self.f_num_x as f32 - 1.0) * h) as f32;
            y = clamp(x, h, (self.f_num_y as f32 - 1.0) * h) as f32;

            let x0 = ((x - h2) * h1).floor();
            let tx = ((x - h2) - x0 * h) * h1;
            let x1 = f32::min(x0 + 1.0, (self.f_num_x - 2) as f32);

            let y0 = ((y - h2) * h1).floor();
            let ty = ((y - h2) - y0 * h) * h1;
            let y1 = f32::min(y0 + 1.0, self.f_num_y as f32 - 2.0);

            let sx = 1.0 - tx;
            let sy = 1.0 - ty;

            if x0 < self.f_num_x as f32 && y0 < self.f_num_y as f32 {
                d[(x0 as usize) * (n as usize) + (y0 as usize)] += sx * sy;
            }
            if x1 < self.f_num_x as f32 && y0 < self.f_num_y as f32 {
                d[(x1 as usize) * (n as usize) + (y0 as usize)] += tx * sy;
            }
            if x1 < self.f_num_x as f32 && y1 < self.f_num_y as f32 {
                d[(x1 as usize) * (n as usize) + (y1 as usize)] += tx * ty;
            }
            if x0 < self.f_num_x as f32 && y1 < self.f_num_y as f32 {
                d[(x0 as usize) * (n as usize) + (y1 as usize)] += sx * ty;
            }
        }

        if self.particle_rest_density == 0.0 {
            let sum = 0.0;
            let num_fluid_cells = 0;
        }
    }
}
fn main() {
    println!("Hello, world!");
}
