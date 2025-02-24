use rand::{Rng, RngCore};
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

// BOUNCY BALL TEST SIMULATION.
static mut SIM: BounceSim = BounceSim { areaWidth:2.5f32, areaHeight:1.0f32, ballz: Vec::new() };

const INSTANCE_DATA_SIZE: usize = 4 + 2 + 2; //4+2+2
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

// BOUNCY BALL TEST SIMULATION.Vj
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

#[wasm_bindgen]
pub fn update_test_simulation()
{
    SIM.integrate_balls(delta_time);
    SIM.handle_ball_wall_collision();
    SIM.draw_balls();
}