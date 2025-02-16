console.log("Hello from WebGPU js file! :)");
const rand = (min, max) => {

    if (min === undefined)
    {
        min = 0;
        max = 1;
    }
    else if (max === undefined)
    {
        max = min;
        min = 0;
    }

    return min + Math.random() * (max - min);
};

function ortho(left, right, bottom, top, near, far, mat)
{
    mat[0] = 2 / (right - left);
    mat[1] = 0;
    mat[2] = 0;
    mat[3] = 0;

    mat[4] = 0;
    mat[5] = 2 / (top - bottom);
    mat[6] = 0;
    mat[7] = 0;

    mat[8] = 0;
    mat[9] = 0;
    mat[10] = -2 / (far - near);
    mat[11] = 0;

    mat[12] = -(right + left) / (right - left);
    mat[13] = -(top + bottom) / (top - bottom);
    mat[14] = -(far + near) / (far - near);
    mat[15] = 1;
}

function initObjects(objArray, numObjects)
{
    for (let i = 0; i < numObjects; i++)
    {
        objArray.push(
            {
                color: new Float32Array([rand(), rand(), rand(), 1.0]),
                position: new Float32Array([-1 + 2*rand(),-1 + 2*rand()]),
                scale: new Float32Array([0.1, 0.1]),
                velocity: [rand(0.0, 0.0), rand(-0.1, 0.1)]
            }
        )
    }
}
function updateObjects(objArray)
{
    for (const obj of objArray)
    {
        obj.position[0] += obj.velocity[0];
    }
}

function updateInstanceValues(instanceValuesF32, objects)
{
    for (let i = 0; i < objects.length; i++)
    {
        const strideF32 = i * 8; // Stride
        const object = objects[i];
        instanceValuesF32.set([object.color[0], object.color[1], object.color[2], object.color[3]], strideF32 + 0);
        instanceValuesF32.set([object.position[0], object.position[1]], strideF32 + 4);
        instanceValuesF32.set([object.scale[0], object.scale[1]], strideF32 + 6);
    }
}

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
    const circleShaderSrc = await fetch("./shaders/circle.wgsl").then(x => x.text())

    main(device, circleShaderSrc);
}

function main(device, circleShaderSrc)
{
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
    vertexBuf = device.createBuffer({
        label: 'Vertex buffer for objects',
        size: vertexBufferSize,
        usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
    });

    instBuf = device.createBuffer({
        label: 'Instance data buffer for objects',
        size: instBufferSize,
        usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
    });

    indexBuf = device.createBuffer({
        label: 'Index buffer for objects',
        size: indexBufferSize,
        usage: GPUBufferUsage.INDEX | GPUBufferUsage.COPY_DST
    });

    uniformBuf = device.createBuffer({
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

        updateObjects(circleObjs);
        updateInstanceValues(instanceValuesF32, circleObjs);
        device.queue.writeBuffer(instBuf, 0, instanceValuesF32);
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

initWebGPU();