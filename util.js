export const rand = (min, max) => {

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

export function ortho(left, right, bottom, top, near, far, mat)
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

export function view(x, y, mat)
{
    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;
    mat[3] = 0;

    mat[4] = 0;
    mat[5] = 1;
    mat[6] = 0;
    mat[7] = 0;

    mat[8] = 0;
    mat[9] = 0;
    mat[10] = 1;
    mat[11] = 0;

    mat[12] = x;
    mat[13] = y;
    mat[14] = 0;
    mat[15] = 1;
}
export function mat4Identity(mat)
{
    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;
    mat[3] = 0;

    mat[4] = 0;
    mat[5] = 1;
    mat[6] = 0;
    mat[7] = 0;

    mat[8] = 0;
    mat[9] = 0;
    mat[10] = 1;
    mat[11] = 0;

    mat[12] = 0;
    mat[13] = 0;
    mat[14] = 0;
    mat[15] = 1;
}

export function getMousePosNDC(mouse_x, mouse_y, w, h)
{
    return [((mouse_x / w) - 0.5) * 2, ((mouse_y / h) - 0.5) * -2];
}

export function mat4MultV2(m, v)
{
    return [
        m[0] * v[0] + m[4] * v[1] + m[12],
        m[1] * v[0] + m[5] * v[1] + m[13]
    ]
}
export function mulMat4(a, b) {
    return [
        a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12],
        a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13],
        a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14],
        a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15],

        a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12],
        a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13],
        a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14],
        a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15],

        a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12],
        a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13],
        a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14],
        a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15],

        a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12],
        a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13],
        a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14],
        a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15]
    ];
}

export function initObjects(objArray, numObjects)
{
    for (let i = 0; i < numObjects; i++)
    {
        objArray.push(
            {
                color: new Float32Array([1.0, 0.0, 0.0, 1.0]),
                position: new Float32Array([-1 + rand(),-1 + 2*rand()]),
                scale: new Float32Array([0.1, 0.1]),
                velocity: [rand(0.0, 0.0), rand(-0.1, 0.1)]
            }
        )
    }
}
export function updateObjects(objArray)
{
    for (const obj of objArray)
    {
        obj.position[0] += obj.velocity[0];
    }
}

export function updateInstanceValues(instanceValuesF32, objects)
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

export function invertMat4(m) {
    const det = (
        m[0] * m[5] * m[10] * m[15] + m[0] * m[6] * m[11] * m[13] + m[0] * m[7] * m[9] * m[14]
        - m[0] * m[7] * m[10] * m[13] - m[0] * m[6] * m[9] * m[15] - m[0] * m[5] * m[11] * m[14]
        - m[1] * m[4] * m[10] * m[15] - m[1] * m[6] * m[11] * m[12] - m[1] * m[7] * m[8] * m[14]
        + m[1] * m[7] * m[10] * m[12] + m[1] * m[6] * m[8] * m[15] + m[1] * m[4] * m[11] * m[14]
        + m[2] * m[4] * m[9] * m[15] + m[2] * m[5] * m[11] * m[12] + m[2] * m[7] * m[8] * m[13]
        - m[2] * m[7] * m[9] * m[12] - m[2] * m[5] * m[8] * m[15] - m[2] * m[4] * m[11] * m[13]
        - m[3] * m[4] * m[9] * m[14] - m[3] * m[5] * m[10] * m[12] - m[3] * m[6] * m[8] * m[13]
        + m[3] * m[6] * m[9] * m[12] + m[3] * m[5] * m[8] * m[14] + m[3] * m[4] * m[10] * m[13]
    );

    if (det === 0) {
        throw new Error("Matrix is not invertible");
    }

    const invDet = 1 / det;

    const inv = [
        (m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10]) * invDet,
        (-m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10]) * invDet,
        (m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6]) * invDet,
        (-m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6]) * invDet,
        (-m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10]) * invDet,
        (m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10]) * invDet,
        (-m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6]) * invDet,
        (m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6]) * invDet,
        (m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9]) * invDet,
        (-m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9]) * invDet,
        (m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5]) * invDet,
        (-m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5]) * invDet,
        (-m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9]) * invDet,
        (m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9]) * invDet,
        (-m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5]) * invDet,
        (m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5]) * invDet
    ];

    return inv;
}