import numpy as np
from mayavi import mlab
from FDTD import fdtd

# Settings
dt = 0.1
dx = 0.5
steps = 200;
points_x = 50;
points_y = 50;
points_z = 50;

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

# Plot using mayavi
mlab.figure(bgcolor=(0,0,0))
sx = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="x_axes")
sy = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="y_axes")
sz = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="z_axes")
mlab.colorbar()

mlab.figure(bgcolor=(0,0,0))
s = mlab.imshow(np.abs(ez[0,:,:,20]), vmax=0.05)
mlab.colorbar()
@mlab.animate(delay=30)
def anim():
    i = 0
    while True:
        i%=steps
        s.mlab_source.scalars = np.abs(ez[i,:,:,20])
        i+=1
        yield
anim()
mlab.show()
