import numpy as np
cimport numpy as np

np.import_array()
DTYPE = np.double
ctypedef np.double_t DTYPE_t

def fdtd(double dt, double dx, int steps, int points_x, int points_y, int points_z,
         np.ndarray [DTYPE_t, ndim=4] ex, np.ndarray [DTYPE_t, ndim=4] ey, np.ndarray [DTYPE_t, ndim=4] ez,
         np.ndarray [DTYPE_t, ndim=4] hx, np.ndarray [DTYPE_t, ndim=4] hy, np.ndarray [DTYPE_t, ndim=4] hz,         
         np.ndarray [DTYPE_t, ndim=3] mu, np.ndarray [DTYPE_t, ndim=3] ep, np.ndarray [DTYPE_t, ndim=3] co):

    cdef int t, x, y, z

    for t in range(1, steps):
        for x in range(0, points_x):
            for y in range(0, points_y - 1):
                for z in range(0, points_z - 1):
                    hx[t, x, y, z] = hx[t - 1, x, y, z] + dt / (mu[x, y, z] * dx) * ((ey[t - 1, x, y, z + 1] - ey[t - 1, x, y, z]) - (ez[t - 1, x, y + 1, z] - ez[t - 1, x, y, z]))

        for x in range(0, points_x - 1):
            for y in range(0, points_y):
                for z in range(0, points_z - 1):
                    hy[t, x, y, z] = hy[t - 1, x, y, z] + dt / (mu[x, y, z] * dx) * ((ez[t - 1, x + 1, y, z] - ez[t - 1, x, y, z]) - (ex[t - 1, x, y, z + 1] - ex[t - 1, x, y, z]))

        for x in range(0, points_x - 1):
            for y in range(0, points_y - 1):
                for z in range(0, points_z):
                    hz[t, x, y, z] = hz[t - 1, x, y, z] + dt / (mu[x, y, z] * dx) * ((ex[t - 1, x, y + 1, z] - ex[t - 1, x, y, z]) - (ey[t - 1, x + 1, y, z] - ey[t - 1, x, y, z]))

        for x in range(0, points_x):
            for y in range(1, points_y):
                for z in range(1, points_z):
                    ex[t, x, y, z] = (1 - co[x, y, z] * dt / ep[x,y,z] * 0.5) / (1 + co[x, y, z] * dt / ep[x,y,z] * 0.5) * ex[t - 1, x, y, z] + 1 / (1 + co[x, y, z] * dt / ep[x, y, z] * 0.5) * dt / (ep[x, y, z] * dx) * ((hz[t, x, y, z] - hz[t, x, y - 1, z]) - (hy[t, x, y, z] - hy[t, x, y, z - 1]))

        for x in range(1, points_x):
            for y in range(0, points_y):
                for z in range(1, points_z):
                    ey[t, x, y, z] = (1 - co[x, y, z] * dt / ep[x,y,z] * 0.5) / (1 + co[x, y, z] * dt / ep[x,y,z] * 0.5) * ey[t - 1, x, y, z] + 1 / (1 + co[x, y, z] * dt / ep[x, y, z] * 0.5) * dt / (ep[x, y, z] * dx) * ((hx[t, x, y, z] - hx[t, x, y, z - 1]) - (hz[t, x, y, z] - hz[t, x - 1, y, z]))

        for x in range(1, points_x):
            for y in range(1, points_y):
                for z in range(0, points_z):
                    ez[t, x, y, z] = (1 - co[x, y, z] * dt / ep[x,y,z] * 0.5) / (1 + co[x, y, z] * dt / ep[x,y,z] * 0.5) * ez[t - 1, x, y, z] + 1 / (1 + co[x, y, z] * dt / ep[x, y, z] * 0.5) * dt / (ep[x, y, z] * dx) * ((hy[t, x, y, z] - hy[t, x - 1, y, z]) - (hx[t, x, y, z] - hx[t, x, y - 1, z]))

        ez[t, 50, 50, 50] += np.sin(t * dt)