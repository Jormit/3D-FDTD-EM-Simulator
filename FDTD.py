import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from FDTD import fdtd

# Settings
dt = 0.1
dx = 0.5
steps = 200;
points_x = 100;
points_y = 100;
points_z = 100;

# Fields
ex = np.zeros((steps, points_x, points_y, points_z));
ey = np.zeros((steps, points_x, points_y, points_z));
ez = np.zeros((steps, points_x, points_y, points_z));

hx = np.zeros((steps, points_x, points_y, points_z));
hy = np.zeros((steps, points_x, points_y, points_z));
hz = np.zeros((steps, points_x, points_y, points_z));


# Material Properties
mu = np.ones((points_x, points_y, points_z)) # Permeability
ep = np.ones((points_x, points_y, points_z)) # Permittivity
co = np.zeros((points_x, points_y, points_z)) # Conductivity

# Pass to C code for main processing loop
fdtd(dt, dx, steps, points_x, points_y, points_z, ex, ey, ez, hx, hy, hz, mu, ep, co)

plt.imshow(ez[steps-1, 45, :, :], animated=True, cmap='jet')
plt.show()