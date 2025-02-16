import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";

const runWasm = async () =>
{
    const rustWasm = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");

    // providing mouse movement to WASM
    window.addEventListener('mousemove', function (e)
    {
        console.log(e.x, e.y);
        console.log(canvas.width, canvas.height);
    });
};
runWasm();