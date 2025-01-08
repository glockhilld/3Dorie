import numpy as np

# TODO: considering rotation center and expend to 3D cases
# TODO: add the siddon projection
def parallel_ray(
    nray,
    angles, 
    ray_width, 
    nx, ny, 
    center,
    siddon=False,
    boolean_wieght=False
    ):
    '''parallel ray projection geometry'''
    nangle = angles.size
    axis_x = np.linspace(-1/2*(nx-1), 1/2*(nx-1), nx)
    axis_y = np.linspace(-1/2*(ny-1), 1/2*(ny-1), ny)
    
    # create the projection matrix
    detector = np.linspace((nray-1)*.5, (nray-1)*.5, nray)

    row = []
    col = []
    vals = []
    for angle in angles:
        cos = np.cos(angle)
        sin = np.sin(angle)
        v_det = np.array([cos, sin])
        v_ray = np.array([-sin, cos])
        _axis_x = axis_x * cos
        _axis_y = axis_y * sin
        
        for i in range(nray):
            _distance = np.linalg.norm((_axis_x[i], _axis_y[i]))
            for j in range(nx):
                for k in range(nx):
                    pos_pixel = np.array([axis_x[j], axis_y[k]])
                    # get distance to ray vector 
                    # image pixel vector projection on the rotated detecotr vector
                    _distance_pixel = np.abs(np.dot(pos_pixel, v_det))
                    _dist = np.linalg.norm(_distance_pixel - _distance)
                    if _dist <= 1/2*ray_width:
                        row.append(j)
                        col.append(k)
                        vals.append(_dist)
    if boolean_wieght:
        vals = np.where(np.array(vals) >0, 1, 0)

    return row, col, vals



