struct vsOutput
{
    @builtin(position) position: vec4f,
    @location(0) color: vec4f,
    @location(1) texCoord: vec2f
};

struct Uniforms
{
    proj: mat4x4f
}

struct vert
{
    @location(0) vertexPos: vec2f,    // Position in model space.
    @location(1) color: vec4f,        // Instance: Color
    @location(2) texCoord: vec2f,     // UV Coords:
    @location(3) worldPos: vec2f,     // Position of object in world
    @location(4) scale: vec2f         // Local scale of object
};

@group(0) @binding(0) var<uniform> uniforms: Uniforms;

@vertex
fn vs(vertex: vert) -> vsOutput
{
    var model = mat4x4f(
        vertex.scale.x, 0.0, 0.0, 0,
        0.0, vertex.scale.y, 0.0, 0,
        0.0, 0.0, 1.0, 0.0,
        vertex.worldPos.x, vertex.worldPos.y, 0.0, 1.0
    );
    var proj = uniforms.proj;
    var vsOut: vsOutput;
    var mp = proj * model;
    vsOut.position = mp * vec4f(vertex.vertexPos, 0.0, 1.0);
    vsOut.color = vertex.color;
    vsOut.texCoord = vertex.texCoord;
    return vsOut;
}

@fragment
fn fs(fsInput : vsOutput) -> @location(0) vec4f
{
    var uv = fsInput.texCoord;
    uv -= 0.5;
    uv *= 2;

    var distance = 1.0 - length(uv);
    if (distance > 0.0)
    {
        distance = 1.0;
    }
    var color = vec3f(distance) * fsInput.color.xyz;
    return vec4f(color, 0.0);
}