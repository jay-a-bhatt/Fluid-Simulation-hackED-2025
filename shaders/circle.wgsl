struct vsOutput
{
    @builtin(position) position: vec4f,
    @location(0) color: vec4f,
    @location(1) texCoord: vec2f
};

struct vert
{
    @location(0) vertexPos: vec2f,    // Position in model space.
    @location(1) color: vec4f,        // Instance: Color
    @location(2) texCoord: vec2f,     // UV Coords:
    @location(3) worldPos: vec2f,     // Position of object in world
    @location(4) scale: vec2f         // Local scale of object
};

@vertex
fn vs(vertex: vert) -> vsOutput
{
    var vsOut: vsOutput;
    vsOut.position = vec4f((vertex.vertexPos * vertex.scale) + vertex.worldPos, 0.0, 1.0);
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
    return vec4f(color, 1.0);
}