import numpy as np
cimport numpy as np

np.import_array()
ctypedef np.double_t DTYPE_t

def fdtd(double dt, double dx, int steps, int points_x, int points_y, int points_z,
         np.ndarray [DTYPE_t, ndim=4] ex, np.ndarray [DTYPE_t, ndim=4] ey, np.ndarray [DTYPE_t, ndim=4] ez,
         np.ndarray [DTYPE_t, ndim=4] hx, np.ndarray [DTYPE_t, ndim=4] hy, np.ndarray [DTYPE_t, ndim=4] hz,         
         np.ndarray [DTYPE_t, ndim=3] mu, np.ndarray [DTYPE_t, ndim=3] ep, np.ndarray [np.uint8_t, ndim=3] co,
         np.ndarray [DTYPE_t, ndim=1] s , np.ndarray [np.int64_t, ndim=1]  sx, np.ndarray [np.int64_t, ndim=1]  sy, np.ndarray [np.int64_t, ndim=1] sz):

    cdef int t, x, y, z
    cdef DTYPE_t Sc = dt/dx
    cdef np.ndarray[DTYPE_t, ndim=3] ch = np.zeros([points_x, points_y, points_z])
    cdef np.ndarray[DTYPE_t, ndim=3] ce = np.zeros([points_x, points_y, points_z])
    cdef np.ndarray[DTYPE_t, ndim=3] cd = np.zeros([points_x, points_y, points_z])
    cdef np.ndarray[DTYPE_t, ndim=3] ca = np.zeros([points_x, points_y, points_z])

    for x in range(0, points_x):
        for y in range(0, points_y):
            for z in range(0, points_z):
                ch[x,y,z] = dt / (mu[x, y, z] * dx)                    
                ca[x,y,z] = (Sc / np.sqrt(ep[x, y, z] * mu[x, y, z]) - 1) / (Sc / np.sqrt(ep[x, y, z] * mu[x, y, z]) + 1)
                if (co[x,y,z] == 1):
                    ce[x,y,z] = 0
                    cd[x,y,z] = 0
                else:
                    ce[x,y,z] = dt / (ep[x, y, z] * dx)
                    cd[x,y,z] = 1

    for t in range(1, steps):
        for x in range(0, points_x):
            for y in range(0, points_y - 1):
                for z in range(0, points_z - 1):
                    hx[t, x, y, z] = hx[t - 1, x, y, z] + ch[x,y,z] * ((ey[t - 1, x, y, z + 1] - ey[t - 1, x, y, z]) - (ez[t - 1, x, y + 1, z] - ez[t - 1, x, y, z]))
        hx[t, :, points_y - 1, :] = hx[t - 1, :, points_y - 2, :] + ca[x,y,z] * (hx[t, :, points_y - 2, :] - hx[t - 1, :, points_y - 1, :])
        hx[t, :, :, points_z - 1] = hx[t - 1, :, :, points_z - 2] + ca[x,y,z] * (hx[t, :, :, points_z - 2] - hx[t - 1, :, :, points_z - 1])

        for x in range(0, points_x - 1):
            for y in range(0, points_y):
                for z in range(0, points_z - 1):
                    hy[t, x, y, z] = hy[t - 1, x, y, z] + ch[x,y,z] * ((ez[t - 1, x + 1, y, z] - ez[t - 1, x, y, z]) - (ex[t - 1, x, y, z + 1] - ex[t - 1, x, y, z]))
        hy[t, points_x - 1, :, :] = hy[t - 1, points_x - 2, :, :] + ca[x,y,z] * (hy[t, points_x - 2, :, :] - hy[t - 1, points_x - 1, :, :])
        hy[t, :, :, points_z - 1] = hy[t - 1, :, :, points_z - 2] + ca[x,y,z] * (hy[t, :, :, points_z - 2] - hy[t - 1, :, :, points_z - 1])

        for x in range(0, points_x - 1):
            for y in range(0, points_y - 1):
                for z in range(0, points_z):
                    hz[t, x, y, z] = hz[t - 1, x, y, z] + ch[x,y,z] * ((ex[t - 1, x, y + 1, z] - ex[t - 1, x, y, z]) - (ey[t - 1, x + 1, y, z] - ey[t - 1, x, y, z]))
        hz[t, points_x - 1, :, :] = hz[t - 1, points_x - 2, :, :] + ca[x,y,z] * (hz[t, points_x - 2, :, :] - hz[t - 1, points_x - 1, :, :])
        hz[t, :, points_y - 1, :] = hz[t - 1, :, points_y - 2, :] + ca[x,y,z] * (hz[t, :, points_y - 2, :] - hz[t - 1, :, points_y - 1, :])

        for x in range(0, points_x):
            for y in range(1, points_y):
                for z in range(1, points_z):
                    ex[t, x, y, z] = cd[x,y,z] * ex[t - 1, x, y, z] + ce[x,y,z] * ((hz[t, x, y, z] - hz[t, x, y - 1, z]) - (hy[t, x, y, z] - hy[t, x, y, z - 1]))
        ex[t, :, 0, :] = ex[t - 1, :, 1, :] + ca[x,y,z] * (ex[t, :, 1, :] - ex[t - 1, :, 0, :])
        ex[t, :, :, 0] = ex[t - 1, :, :, 1] + ca[x,y,z] * (ex[t, :, :, 1] - ex[t - 1, :, :, 0])

        for x in range(1, points_x):
            for y in range(0, points_y):
                for z in range(1, points_z):
                    ey[t, x, y, z] = cd[x,y,z] * ey[t - 1, x, y, z] + ce[x,y,z] * ((hx[t, x, y, z] - hx[t, x, y, z - 1]) - (hz[t, x, y, z] - hz[t, x - 1, y, z]))
        ey[t, 0, :, :] = ey[t - 1, 1, :, :] + ca[x,y,z] * (ey[t, 1, :, :] - ey[t - 1, 0, :, :])
        ey[t, :, :, 0] = ey[t - 1, :, :, 1] + ca[x,y,z] * (ey[t, :, :, 1] - ey[t - 1, :, :, 0])

        for x in range(1, points_x):
            for y in range(1, points_y):
                for z in range(0, points_z):
                    ez[t, x, y, z] = cd[x,y,z] * ez[t - 1, x, y, z] + ce[x,y,z] * ((hy[t, x, y, z] - hy[t, x - 1, y, z]) - (hx[t, x, y, z] - hx[t, x, y - 1, z]))
        ez[t, 0, :, :] = ez[t - 1, 1, :, :] + ca[x,y,z] * (ez[t, 1, :, :] - ez[t - 1, 0, :, :])
        ez[t, :, 0, :] = ez[t - 1, :, 1, :] + ca[x,y,z] * (ez[t, :, 1, :] - ez[t - 1, :, 0, :])

        ez[t, sx, sy, sz] = s[t]