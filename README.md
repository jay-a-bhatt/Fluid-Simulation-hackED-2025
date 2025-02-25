<div align="center">

  <h1>🌊 Fluid Simulation🌊</h1>

  <p>
    <strong>PIC/FLIP fluid simulation written in Rust, targeting WebAssembly and rendered with WebGPU</strong>
  </p>
</div>

## About

This project was created for the HackED 2025 hackathon and since been expanded upon. It is a Rust + Wasm implementation of a PIC/FLIP fluid simulation created by [Matthias Müller](https://youtu.be/XmzBREkK8kY?si=5a8RqvdLVUWC9ErL). The simulation is interactable and features different sliders and toggles to change the behavior of the fluid.
<div align="center">
<img aligne = "center" src="./res/fluid_demo.gif"/>
</div>

## 📝 Prerequisites
- Rust 1.30.0 or later.
- [wasm-pack](https://rustwasm.github.io/wasm-pack/installer)
- Python3 (for http server)
- Browser with WebGPU support.

## ⚙ Building and Running

- Clone the repository
- Inside the project directory run `wasm-pack build --target web`
- After building, create a local http server by running `python -m http.server 8000`
- In your browser go to http://localhost:8000
