import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";

const runWasm = async () => {
  // instantiate wasm module as ?
  const rustWasm = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");

// referring to canvas, calling wasm with width & height

function context()
{
  return [canvas.width, canvas.height]
}


};
runWasm();

export default context;