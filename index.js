import wasmInit from "./pkg/Fluid_Simulation_hackED_2025.js";
import {SimWASM} from "./pkg/Fluid_Simulation_hackED_2025.js";
import {ortho, initObjects, updateObjects, updateInstanceValues} from './util.js'
import {testStruct} from "./pkg/Fluid_Simulation_hackED_2025.js";
import {printBalls} from "./util.js";


// render global state

// function here

const num = 0;
window.exports = {
    hello: function () { console.log(num); }
};


async function initWebGPU()
{
    if (!navigator.gpu)
    {
        console.error('This browser does not support WebGPU');
        return;
    }
    const adapter = await navigator.gpu?.requestAdapter();
    if (!adapter)
    {
        console.error('No WebGPU adapters found');
        return;
    }

    let device = await adapter.requestDevice();
    device.lost.then((info) =>
    {
        console.error(`WebGPU device was lost: ${info.message}`)
        device = null;
        if (info.reason != 'destroyed')
        {
            initWebGPU();
        }
    });

    console.log("WebGPU Initialized Successfully");
    return device;
}


function setupScene(){

 /*   var scene = 
	{
		gravity : -9.81,
//		gravity : 0.0,
		dt : 1.0 / 120.0,
		flipRatio : 0.9,
		numPressureIters : 100,
		numParticleIters : 2,
		frameNr : 0,
		overRelaxation : 1.9,
		compensateDrift : true,
		separateParticles : true,
		obstacleX : 0.0,
		obstacleY : 0.0,
		obstacleRadius: 0.15,
		paused: true,
		showObstacle: true,
		obstacleVelX: 0.0,
		obstacleVelY: 0.0,
		showParticles: true,
		showGrid: false,
		fluid: null
	};

    scene.obstacleRadius = 0.15;
    scene.overRelaxation = 1.9;

    scene.dt = 1.0 / 60.0;
    scene.numPressureIters = 50;
    scene.numParticleIters = 2;*/

    var simHeight = 3.0;
    var canvasScale = canvas.height/simHeight;
    var simWidth = canvas.width/canvasScale;

    var tankWidth = 1.0 * simWidth;
    var tankHeight = 1.0 * simHeight;
    var relativeWaterHeight = 0.8;
    var relativeWaterWidth = 0.6;
    var density = 1000.0;

    var res = 100;
    var spacing = tankHeight / res; //h

    //dam break

    //compute number of particles

    var particleRadius = 0.3 * spacing; // r
    var dx = 2.0 * particleRadius;
    var dy = Math.sqrt(3.0) / (2.0 * dx);

    var numX = Math.floor((relativeWaterWidth * tankWidth - 2.0 * spacing - 2.0 * particleRadius) /dx);
    var numY =  Math.floor((relativeWaterHeight * tankHeight - 2.0 * spacing - 2.0 * particleRadius) /dy);
    var maxParticles = numX * numY;


    // create fluid

    var simWasm = new SimWASM(density, tankWidth, tankHeight, spacing, particleRadius, maxParticles); // f

    simWasm.scene.fluid.num_particles= numX * numY
    var p = 0;
    for (let i = 0; i < numX; i++) {
        for (let j = 0; j < numY; j++) {
            simWasm.scene.fluid.particle_pos[p++] = spacing + particleRadius + dx * i + (j % 2 == 0 ? 0.0 : particleRadius);
            simWasm.scene.fluid.particle_pos[p++] = spacing + particleRadius + dy * j;
        }
    }

    //setup grid cells for tank

    var n = simWasm.scene.fluid.f_num_y;

    for (let i = 0; i < simWasm.scene.fluid.f_num_x; i++) {
        for (let j = 0; j < simWasm.scene.fluid.f_num_y; j++) {
            var s = 1.0; // fluid
            if(i == 0 || i == simWasm.scene.fluid.f_num_x-1 || j == 0){
                s = 0.0;
            }
            simWasm.scene.fluid.s[i*n + j] = s;
        }
    }
   // setObstacle(3.0, 2.0, true);
}

function main(device, simModule, circleShaderSrc)
{
    let wasmMemory = new Uint8Array(simModule.memory.buffer);

    // Start the simulation
    simModule.init_simulation();
    const canvas = document.querySelector('canvas');
    const presentationFmt = navigator.gpu.getPreferredCanvasFormat();
    const context = canvas.getContext('webgpu');
    context.configure({device, format: presentationFmt} )

   
    // Circle Render Data ------
    const circleShaderModule = device.createShaderModule( {label: 'Circle Shader Module', code: circleShaderSrc})
    if (!circleShaderModule) { console.error("Failed to create circle shader module!"); }

    const maxObjects = 100;
    const circleObjs = [];
    const projectionMat = new Float32Array(16);

    initObjects(circleObjs, maxObjects);
    console.log(circleObjs);
    // Vertex Data
    const numVertices = 4;
    const numTris = 2;
    const vertexData = new Float32Array([
        //  x     y
        -0.5, -0.5,   0, 0, // BL
        0.5, -0.5,   1, 0,  // BR
        0.5,  0.5,   1, 1,  // TR
        -0.5,  0.5,   0, 1  // TL
    ]);

    // 3--2       2
    // | /       /|
    // |/       / |
    // 0       0--1
    const indexData = new Uint32Array([0, 1, 2,   0, 2, 3]);

    const vertexSize   = 4 * 4; // 4 x 4byte floats
    const instUnitSize = (4 + 2 + 2) * 4; // 8 x 4byte floats

    const vertexLayout = [
        {arrayStride: vertexSize,   stepMode: 'vertex',   attributes: [ {shaderLocation: 0, offset:  0, format: 'float32x2'},
                {shaderLocation: 2, offset:  8, format: 'float32x2'} ]}, // Per Vertex Position
        {arrayStride: instUnitSize, stepMode: 'instance', attributes: [ {shaderLocation: 1, offset:  0, format: 'float32x4'},   // Per Instance Color
                {shaderLocation: 3, offset: 16, format: 'float32x2'},  // Per Instance Offset
                {shaderLocation: 4, offset: 24, format: 'float32x2'}]}  // Per Instance Scale
    ];

    const vertexBufferSize = vertexData.byteLength;
    const instBufferSize   = instUnitSize * 100;
    const indexBufferSize  = indexData.byteLength;

    // create buffers
    const vertexBuf = device.createBuffer({
        label: 'Vertex buffer for objects',
        size: vertexBufferSize,
        usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
    });

    const instBuf = device.createBuffer({
        label: 'Instance data buffer for objects',
        size: instBufferSize,
        usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
    });

    const indexBuf = device.createBuffer({
        label: 'Index buffer for objects',
        size: indexBufferSize,
        usage: GPUBufferUsage.INDEX | GPUBufferUsage.COPY_DST
    });

    const uniformBuf = device.createBuffer({
        label: 'Uniform buffer for objects',
        size: 16 * 4,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });

    // write to buffer
    device.queue.writeBuffer(vertexBuf, 0, vertexData);
    device.queue.writeBuffer(indexBuf, 0, indexData);

    // HACK: need to swap this out with real buffer
    const instanceValuesF32 = new Float32Array(instBufferSize / 4);
    updateInstanceValues(instanceValuesF32, circleObjs);
    /*
    for (let i = 0; i < maxObjects; i++)
    {
        const strideF32 = i * 8; // Stride
        instanceValuesF32.set([rand(), rand(), rand(), 1], strideF32 + 0);
        instanceValuesF32.set([(rand()-0.5) * 2, (rand()-0.5) * 2], strideF32 + 4)
        instanceValuesF32.set([0.1,0.1], strideF32 + 6);
    }
     */
    device.queue.writeBuffer(instBuf, 0, instanceValuesF32);

    const circlePipeline = device.createRenderPipeline(
        {
            label: 'Circle Pipeline',
            layout: 'auto',
            vertex:   { module: circleShaderModule, buffers: vertexLayout},
            fragment: { module: circleShaderModule, targets: [ {format: presentationFmt} ]}
        });

    const renderPassDescriptor = {
        label: 'main canvas renderer',
        colorAttachments: [ {clearValue: [0.2,0.2, 0.3, 1], loadOp: 'clear', storeOp: 'store'} ]
    }

    const bindGroup = device.createBindGroup({
        label: 'circle bind group',
        layout: circlePipeline.getBindGroupLayout(0),
        entries: [ {binding: 0, resource: {buffer: uniformBuf} }]
    });
    function render()
    {
        // Update Simulation State
        //fluidSim.update(0.001);
       // simModule.update(0.01, testS);
        // Get pointer to location of instance buffer in wasm memory
        let instanceBufferPtr = simModule.get_instance_buffer_ptr();
        // Get a F32 array view into the buffer
        let instanceBuffer = new Float32Array(simModule.memory.buffer, instanceBufferPtr, 800);

        //updateObjects(circleObjs);
        //updateInstanceValues(instanceValuesF32, circleObjs);
        device.queue.writeBuffer(instBuf, 0, instanceBuffer);
        const aspect = canvas.width/canvas.height;
        const zoom = 3;
        const l = (-aspect/2) * zoom;
        const r = -l;
        const t = zoom/2;
        const b = -t;
        ortho(l, r, b, t, 200, -100, projectionMat);
        device.queue.writeBuffer(uniformBuf, 0, projectionMat);
        //
        // Set canvas as texture to render too.
        renderPassDescriptor.colorAttachments[0].view = context.getCurrentTexture().createView();
        const encoder = device.createCommandEncoder({ label: 'epic encoder'});

        const pass = encoder.beginRenderPass(renderPassDescriptor);

        pass.setPipeline(circlePipeline);
        pass.setBindGroup(0, bindGroup);
        pass.setVertexBuffer(0, vertexBuf);
        pass.setVertexBuffer(1, instBuf);
        pass.setIndexBuffer(indexBuf, 'uint32');
        pass.drawIndexed(6, maxObjects);

        pass.end();
        const commandBuffer = encoder.finish();

        // Submit To GPUUUUU
        device.queue.submit([commandBuffer]);

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
}

// returning mouse x, y positions
function mouse_position_x()
{
    window.addEventListener('mousemove', function (e)
    {
        return e.x;
    });
}

function mouse_position_y()
{
    window.addEventListener('mousemove', function (e)
    {
        return e.x;
    });
}

const initWasm = async () =>
{
    // Load Wasm module so we can call Rust functions.
    const wasmModule = await wasmInit("./pkg/Fluid_Simulation_hackED_2025_bg.wasm");
    // Get the webGPU adapter device. (need it for graphics stuff)
    const device = await initWebGPU();
    // Fetch the circle shader .wgsl file and turn it into a string for GPU compiling.
    const circleShaderSrc = await fetch("./shaders/circle.wgsl").then(x => x.text())

    main(device, wasmModule, circleShaderSrc);
};

initWasm()
