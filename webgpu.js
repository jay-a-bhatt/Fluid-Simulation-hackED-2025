export async function initWebGPU()
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

export async function compileShaders(device)
{
    const circleShaderSrc = await fetch("./shaders/circle.wgsl").then(x => x.text())
    const squareShaderSrc = await fetch("./shaders/square.wgsl").then(x => x.text())

    const circleShaderModule = device.createShaderModule( {label: 'Circle Shader Module', code: circleShaderSrc})
    if (!circleShaderModule) { console.error("Failed to create circle shader module!"); }


    const squareShaderModule = device.createShaderModule( {label: 'Square Shader Module', code: squareShaderSrc})
    if (!squareShaderModule) { console.error("Failed to create square shader module!"); }

    return {circleShader: circleShaderModule, squareShader: squareShaderModule}
}

const uvQuadVertexData = new Float32Array([
    // x     y
    -0.5, -0.5,   0, 0,  // BL
     0.5, -0.5,   1, 0,  // BR
     0.5,  0.5,   1, 1,  // TR
    -0.5,  0.5,   0, 1   // TL
]);

// 3--2       2
// | /       /|
// |/       / |
// 0       0--1
const quadIndexData = new Uint32Array([0, 1, 2,   0, 2, 3]);

export function createCircleBuffers(device, maxCircleInstances)
{

    const vertexSize   = 4 * 4;
    const instUnitSize = (4 + 2 + 2) * 4;

    const vertexLayout = [
        {arrayStride: vertexSize,   stepMode: 'vertex',   attributes: [
                {shaderLocation: 0, offset:  0, format: 'float32x2'},    // Per Vertex Position
                {shaderLocation: 1, offset:  8, format: 'float32x2'} ]}, // Per Vertex UV Coordinates
        {arrayStride: instUnitSize, stepMode: 'instance', attributes: [
                {shaderLocation: 2, offset:  0, format: 'float32x4'},    // Per Instance Color
                {shaderLocation: 3, offset: 16, format: 'float32x2'},    // Per Instance Offset
                {shaderLocation: 4, offset: 24, format: 'float32x2'}]}    // Per Instance Scale
    ];

    const vertexBufferSize = uvQuadVertexData.byteLength;
    const instBufferSize   = instUnitSize * maxCircleInstances;
    const indexBufferSize  = quadIndexData.byteLength;

    // Create Buffers
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

    // write to buffer
    device.queue.writeBuffer(vertexBuf, 0, uvQuadVertexData);
    device.queue.writeBuffer(indexBuf, 0, quadIndexData);

    return {vertexBuffer: vertexBuf, indexBuffer: indexBuf, instanceBuffer: instBuf, layout: vertexLayout};
}

export function createSquareBuffers(device, maxSquareInstances)
{
    const vertexSize   = 4 * 4;
    const instanceUnitSize = (4 + 2 + 1) * 4; // 4 for color, 2 for pos, 1 for scale

    const vertexLayout = [
        {arrayStride: vertexSize,   stepMode: 'vertex',   attributes: [
                {shaderLocation: 0, offset:  0, format: 'float32x2'},    // Per Vertex Position
                {shaderLocation: 1, offset:  8, format: 'float32x2'} ]}, // Per Vertex UV Coordinates
        {arrayStride: instanceUnitSize, stepMode: 'instance', attributes: [
                {shaderLocation: 2, offset:  0, format: 'float32x4'},    // Per Instance Color
                {shaderLocation: 3, offset: 16, format: 'float32x2'},    // Per Instance Offset
                {shaderLocation: 4, offset: 24, format: 'float32'}]}     // Per Instance Scale
    ];

    const vertexBufferSize = uvQuadVertexData.byteLength;
    const instBufferSize   = instanceUnitSize * maxSquareInstances;
    const indexBufferSize  = quadIndexData.byteLength;

    // Create Buffers
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

    // Write static vertex and index values to buffers
    device.queue.writeBuffer(vertexBuf, 0, uvQuadVertexData);
    device.queue.writeBuffer(indexBuf, 0, quadIndexData);

    return {vertexBuffer: vertexBuf, indexBuffer: indexBuf, instanceBuffer: instBuf, layout: vertexLayout}
}

