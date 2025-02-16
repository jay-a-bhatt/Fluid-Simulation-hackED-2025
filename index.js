import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";

const runWasm = async () => {
// instantiate wasm module as ?
const rustWasm = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");

// referring to canvas, calling wasm with width & height
function context()
{
  return [canvas.width, canvas.height]
}

const buffer_index = rustWasm.return_pointer() / Float32Array.BYTES_PER_ELEMENT; // ?????
for (let i = 0; i < 12; i++)
{
  console.log([buffer_index + i]);
}

};

runWasm();

//export default context;