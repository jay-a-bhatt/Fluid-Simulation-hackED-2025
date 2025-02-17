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