struct vsOutput
{
    @builtin(position) position: vec4f,
    @location(0) color: vec4f,
    @location(1) texCoord: vec2f
};

struct Uniforms
{
    proj: mat4x4f,
    view: mat4x4f
}

struct vert
{
    @location(0) vertexPos: vec2f,    // Position in model space.
    @location(1) texCoord: vec2f,     // UV Coordinates
    @location(2) color: vec4f,        // Per Instance: Color
    @location(3) worldPos: vec2f,     // Per Instance: Position of object in world
    @location(4) scale: vec2f         // Per Instance: Local scale of object
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
    var view = uniforms.view;
    var mvp = proj * view * model;

    var vsOut: vsOutput;

    vsOut.position = mvp * vec4f(vertex.vertexPos, 0.0, 1.0);
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
    else
    {
        discard;
    }
    var color = vec3f(distance) * fsInput.color.xyz;
    return vec4f(color, 0.0);
}