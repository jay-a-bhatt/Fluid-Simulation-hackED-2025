import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";

const runWasm = async () => {
  // instantiate wasm module
  const rustWasm = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");
const canvas = document.querySelector('canvas');


};
runWasm();