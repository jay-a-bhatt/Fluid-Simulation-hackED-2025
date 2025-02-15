struct FlipFluid {
    density: i32,
    fNumX: f32,
    fNumY: f32,
    h: i32,
    fInvSpacing: f32,
    fNumCells: i32,

    u: Vec<f32>,
    v: Vec<f32>,
    du: Vec<f32>,
    dv: Vec<f32>,
    prevU: Vec<f32>,
    prevV: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>,
    cellType: Vec<f32>,
    cellColour: Vec<f32>,

    maxParticles: i32,

    particlePos: Vec<f32>,
    particleColour: Vec<f32>,

    particleVel: Vec<f32>,
    particleDensity: Vec<f32>,
    particleRestDensity: Vec<f32>,

    particleRadius: i32,
    pInvSpacing: f32,
    pNumX: f32,
    pNumY: f32,
    pNumCells: f32,

    numCellParticles: Vec<i32>,
    firstCellParticles: Vec<i32>,
    cellParticleIds: Vec<i32>,

    numParticles: i32,
}

impl FlipFluid {}

fn main() {
    println!("Hello, world!");
}
