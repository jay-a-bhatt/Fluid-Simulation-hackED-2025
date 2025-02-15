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
    @location(2) texCoord: vec2f      // UV Coords:
};

@vertex
fn vs(vertex: vert) -> vsOutput
{
    var vsOut: vsOutput;
    vsOut.position = vec4f(vertex.vertexPos, 0.0, 1.0);
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
    return vec4f(uv, 0.0, 1.0);
}
