struct vsOutput
{
    @builtin(position) position: vec4f,
    @location(0) color: vec4f
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
    vsOut.color = vec4f(vertex.texCoord, 1, 1);
    return vsOut;
}

@fragment
fn fs(fsInput : vsOutput) -> @location(0) vec4f
{
    return fsInput.color;
}
