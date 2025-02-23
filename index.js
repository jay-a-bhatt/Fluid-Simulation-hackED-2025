import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";
import {SimulationHandler} from "./pkg/Fluid_Simulation_hackED_2025.js";
import {
    ortho,
    initObjects,
    updateObjects,
    updateInstanceValues,
    view,
    mat4MultV2,
    getMousePosNDC,
    mat4Identity,
    invertMat4, mulMat4
} from './util.js'
import {initWebGPU, compileShaders, createCircleBuffers, createSquareBuffers} from './webgpu.js'


// render global state
// function here

let context = {wasmModule: null, device: null}

let drawCount = 0;
window.exports =
{
    drawObjects: function drawObject(objectID, count)
    {
        if (objectID == 12)
        {
            drawCount += count
            // Requires simModule, device,
            // Get the buffer ptr,
            //

            // Get Current Pass
            // Set Pipeline
            // Bind Buffers
            // Call Draw Function
            // Profit
        }
    }
}

function initSimulation(simWasmModule)
{
    const testCanvasHeight = 1000;
    const testCanvasWidth = 1000;

    const canvasHeight = testCanvasHeight;
    const canvasWidth = testCanvasWidth;

    const simHeight = 3.0;
    const canvasScale = canvasHeight / simHeight;
    //const simWidth = canvasWidth / canvasScale;
    const simWidth = 3.0;

    // NOTE(Doesnt change)
    const tankWidth = 1.0 * simWidth;
    const tankHeight = 1.0 * simHeight;

    const relWaterHeight = 0.8;
    const relWaterWidth = 0.6;
    const density = 1000.0;

    const res = 100;
    const cellSize = tankHeight / res;

    // Particle Radius is with respect to cell size.
    const particleRadius = 0.3 * cellSize;

    const dx = 2.0 * particleRadius;
    const dy = Math.sqrt(3.0) / 2.0 * dx;

    // Number of potential particles on the X axis.
    const numX = Math.floor((relWaterWidth * tankWidth - 2.0 * cellSize - 2.0 * particleRadius) / dx);
    // Number of potential particles on the y axis.
    const numY = Math.floor((relWaterHeight * tankHeight - 2.0 * cellSize - 2.0 * particleRadius) / dy);

    // Maximum possible particles in our simulation.
    const maxParticles = numX * numY;
    console.log(maxParticles);

    // Create WebAssembly + Rust Simulation Handler
    const simHandler = new SimulationHandler(
        numX,
        numY,
        density,
        tankWidth,
        tankHeight,
        cellSize,
        particleRadius,
        maxParticles
    );

    return simHandler
}

function main(device, simModule, shaders)
{
    let wasmMemory = new Uint8Array(simModule.memory.buffer);

    let simHandler = initSimulation(simModule);
    const maxCircleInstances = 16384;
    const maxSquareInstances = 4096 * 3;
    simModule.init_instance_gfx_buffers(maxCircleInstances, maxSquareInstances);

    const canvas = document.querySelector('canvas');
    const presentationFmt = navigator.gpu.getPreferredCanvasFormat();
    const context = canvas.getContext('webgpu');
    context.configure({device, format: presentationFmt} )

    const circleObjs = [];
    const projectionMat = new Float32Array(16);
    const viewMat = new Float32Array(16);

    initObjects(circleObjs, maxCircleInstances);

    const circleBuffers = createCircleBuffers(device, maxCircleInstances);
    // TODO: fill this out
    const squareBuffers = createSquareBuffers(device, maxSquareInstances);

    // Uniform Buffer that will hold projection and view matrix.
    const matrixUniformBuffer = device.createBuffer({
        label: 'Uniform buffer for projection/view matrices',
        size: (16 * 2) * 4,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });

    // HACK: need to swap this out with real buffer
    const instanceValuesF32 = new Float32Array((4+2+2) * maxCircleInstances);
    updateInstanceValues(instanceValuesF32, circleObjs);

    device.queue.writeBuffer(circleBuffers.instanceBuffer, 0, instanceValuesF32);

    let uniformValues = new Float32Array(32);

    const squarePipeline = device.createRenderPipeline(
        {
            label: 'Square Pipeline',
            layout: 'auto',
            vertex:   { module: shaders.squareShader, buffers: squareBuffers.layout},
            fragment: { module: shaders.squareShader, targets: [ {format: presentationFmt} ] }
        }
    );
    const circlePipeline = device.createRenderPipeline(
    {
        label: 'Circle Pipeline',
        layout: 'auto',
        vertex:   { module: shaders.circleShader, buffers: circleBuffers.layout},
        fragment: { module: shaders.circleShader, targets: [ {format: presentationFmt} ]}
    });

    const bindGroup = device.createBindGroup({
        label: 'circle bind group',
        layout: circlePipeline.getBindGroupLayout(0),
        entries: [ { binding: 0, resource: {buffer: matrixUniformBuffer} }]
    });

    const sqrBindGroup = device.createBindGroup({
        label: 'square bind group',
        layout: squarePipeline.getBindGroupLayout(0),
        entries: [ {binding: 0, resource: {buffer: matrixUniformBuffer} }]
    });

    const renderPassDescriptor = {
        label: 'main canvas renderer',
        colorAttachments: [ {clearValue: [0.2,0.2, 0.3, 1], loadOp: 'clear', storeOp: 'store'} ]
    }

    // MENU INPUT
    let switch_1 = document.getElementById("1"),
        switch_3 = document.getElementById("3"), switch_4 = document.getElementById("4"); // refer to switch_X.querySelector('input').checked (returns T/F bool)
    let slider_1 = document.getElementById("slider1");

    let gravitySlider = document.getElementById("slider2"); // refer to slider_X.value
    let gridSwitch = document.getElementById("2");
    const infoElem = document.querySelector('#info');
    let zoom = 3.0;

    // MOUSE INPUT
    let mouse_x = 0, mouse_y = 0, mouse_down = 0; // mouse_down = {0: false; 1: true}
    window.addEventListener('mousedown', function(event){mouse_down = 1;});
    window.addEventListener('mouseup', function(event){mouse_down = 0;});
    window.addEventListener('mousemove', (event) => { mouse_x = event.clientX; mouse_y = event.clientY;} );
    window.addEventListener("wheel", (event) => {zoom += event.deltaY * 0.002; if (zoom < 0.1) { zoom = 0.1;}});

    let then = 0;
    function render(now)
    {
        now *= 0.001;
        const deltaTime = now - then;
        then = now;

        // Calculate Mouse Position
        const aspect = canvas.width/canvas.height;

        // Create projection and view matrices
        const l = (-aspect/2) * zoom;
        const r = -l;
        const t = zoom/2;
        const b = -t;
        ortho(l, r, b, t, 200, -100, projectionMat);
        view(-aspect/2,-zoom/2.5, viewMat);

        let rect = canvas.getBoundingClientRect();

        let mouse_ndc = getMousePosNDC(mouse_x - rect.left, mouse_y - rect.top, canvas.width, canvas.height);
        let mouse_world = mat4MultV2(invertMat4(mulMat4(viewMat, projectionMat)), mouse_ndc);

        // Renderer Start Pass

        // Update Simulation State
        const simStartTime = performance.now();
        simHandler.update(deltaTime, mouse_world[0], mouse_world[1], gravitySlider.value);
        const simTime = performance.now() - simStartTime;


        simHandler.render();

        // Renderer End Pass
        // Get pointer to location of instance buffer in wasm memory
        let circleInstBufPtr = simModule.get_circle_instance_buf_ptr();
        let squareInstBufPtr = simModule.get_square_instance_buf_ptr();

        // Get a F32 array view into the buffer
        let circleInstanceBufferF32 = new Float32Array(simModule.memory.buffer, circleInstBufPtr, maxCircleInstances * (4 + 2 + 2));
        let squareInstanceBufferF32 = new Float32Array(simModule.memory.buffer, squareInstBufPtr, maxSquareInstances * (4 + 2 + 1));

        device.queue.writeBuffer(circleBuffers.instanceBuffer, 0, circleInstanceBufferF32);
        device.queue.writeBuffer(squareBuffers.instanceBuffer, 0, squareInstanceBufferF32);

        for (let i = 0; i < 16; i++)
        {
            uniformValues[i] = projectionMat[i];
            uniformValues[i + 16] = viewMat[i];
        }
        device.queue.writeBuffer(matrixUniformBuffer, 0, uniformValues);

        // Set canvas as texture to render too.
        renderPassDescriptor.colorAttachments[0].view = context.getCurrentTexture().createView();
        const encoder = device.createCommandEncoder({ label: 'epic encoder'});

        const pass = encoder.beginRenderPass(renderPassDescriptor);

        // Draw Circles
        pass.setPipeline(circlePipeline);
        pass.setBindGroup(0, bindGroup);
        pass.setVertexBuffer(0, circleBuffers.vertexBuffer);
        pass.setVertexBuffer(1, circleBuffers.instanceBuffer);
        pass.setIndexBuffer(circleBuffers.indexBuffer, 'uint32');
        pass.drawIndexed(6, maxCircleInstances);

        // Draw Squares
        pass.setPipeline(squarePipeline);
        pass.setBindGroup(0, sqrBindGroup);
        pass.setVertexBuffer(0, squareBuffers.vertexBuffer);
        pass.setVertexBuffer(1, squareBuffers.instanceBuffer);
        pass.setIndexBuffer(squareBuffers.indexBuffer, 'uint32');

        // Draw Grid if option selected.
        if (gridSwitch.querySelector('input').checked) {
            pass.drawIndexed(6,maxSquareInstances);
        }

        pass.end();

        const commandBuffer = encoder.finish();

        // Submit To GPUUUUU
        device.queue.submit([commandBuffer]);
        infoElem.textContent = `fps: ${(1 / deltaTime).toFixed(1)}\nsim: ${simTime.toFixed(1)}ms\n draws: ${drawCount}`;

        requestAnimationFrame(render);
    }
    render();

    const observer = new ResizeObserver(entries =>
    {
        // Iterates over all entries but there should only be one.
        for (const entry of entries)
        {
            const width = entry.devicePixelContentBoxSize?.[0].inlineSize ||
                entry.contentBoxSize[0].inlineSize * devicePixelRatio;
            const height = entry.devicePixelContentBoxSize?.[0].blockSize ||
                entry.contentBoxSize[0].blockSize * devicePixelRatio;

            const canvas = entry.target;
            canvas.width  = Math.max(1, Math.min(width,  device.limits.maxTextureDimension2D));
            canvas.height = Math.max(1, Math.min(height, device.limits.maxTextureDimension2D));
        }
    });

    try   { observer.observe(canvas, {box: 'device-pixel-content-box'}); }
    catch { observer.observe(canvas, {box: 'content-box'}) }


    /*
    CANVAS VALUE | (canvas.width) x (canvas.height)
           MOUSE |(mouse_x) x (mouse_y); mouse_down
    SWITCH VALUE | ${switch_#.querySelector('input').checked}
    SLIDER VALUE | ${slider_#.value}`)
    */

}

const initWasm = async () =>
{
    // Load Wasm module so we can call Rust functions.
    const wasmModule = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");

    // Get the webGPU adapter device. (need it for graphics stuff)
    const device = await initWebGPU();
    // Fetch the circle shader .wgsl file and turn it into a string for GPU compiling.
    const shaders = await compileShaders(device);

    main(device, wasmModule, shaders);
};

initWasm()
