import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
from FDTD import fdtd

# Settings
f = 1
dr = 6 / f
dt = 0.03 / f
dx = 2 * dt
steps = 1024
points_x = 100
points_y = 50
points_z = 50

t = np.linspace(0, (steps - 1) * dt, steps)
# s = (1-2*(np.pi*f*(t-dr))**2)*np.exp(-(np.pi*f*(t-dr))**2) # Source waveform
s = np.sin(2 * np.pi * t)

# Fields
ex = np.zeros((steps, points_x, points_y, points_z))
ey = np.zeros((steps, points_x, points_y, points_z))
ez = np.zeros((steps, points_x, points_y, points_z))

hx = np.zeros((steps, points_x, points_y, points_z))
hy = np.zeros((steps, points_x, points_y, points_z))
hz = np.zeros((steps, points_x, points_y, points_z))

# Material Properties
mu = np.ones((points_x, points_y, points_z))  # Permeability
ep = np.ones((points_x, points_y, points_z))  # Permittivity
co = np.zeros((points_x, points_y, points_z), dtype=bool)  # Conductivity

# 2 Parallel Conductors
co[:, 22:28, 28:30] = 1
co[:, 22:28, 20:22] = 1

# Source locations
source_pos = np.zeros((points_x, points_y, points_z))
source_pos[0, 22:28, 22:28] = 1
sx, sy, sz = np.where(source_pos == 1)

# Pass to C code for main processing loop
fdtd(
    dt,
    dx,
    steps,
    points_x,
    points_y,
    points_z,
    ex,
    ey,
    ez,
    hx,
    hy,
    hz,
    mu,
    ep,
    co,
    s,
    sx,
    sy,
    sz,
)

# Problem Plot
e_space = ep > 1
c_space = co > 0
s_space = source_pos

mlab.figure(bgcolor=(0, 0, 0))
xx, yy, zz = np.where(e_space == 1)
mlab.points3d(xx, yy, zz, mode="cube", color=(1, 0, 0), scale_factor=1)

xx, yy, zz = np.where(c_space == 1)
mlab.points3d(xx, yy, zz, mode="cube", color=(0.72, 0.45, 0.2), scale_factor=1)

xx, yy, zz = np.where(s_space == 1)
mlab.points3d(xx, yy, zz, mode="cube", color=(0, 0, 1), scale_factor=1)

# Plot 3d Fields
ind = int(steps / 2)
mlab.figure(bgcolor=(0, 0, 0))
sx = mlab.volume_slice(np.abs(ez[ind]), vmax=0.5, plane_orientation="x_axes")
sy = mlab.volume_slice(np.abs(ez[ind]), vmax=0.5, plane_orientation="y_axes")
sz = mlab.volume_slice(np.abs(ez[ind]), vmax=0.5, plane_orientation="z_axes")
mlab.colorbar()
mlab.outline()

# Full Vector Fields
mlab.figure(bgcolor=(0, 0, 0))
mlab.quiver3d(ex[ind], ey[ind], ez[ind])
mlab.title("E Field")
mlab.colorbar()
mlab.outline()

mlab.figure(bgcolor=(0, 0, 0))
mlab.quiver3d(hx[ind], hy[ind], hz[ind])
mlab.title("H Field")
mlab.colorbar()
mlab.outline()

# Masked Vector Fields
mask = np.zeros((points_x, points_y, points_z))
mask[:, 22:28, 22:28] = 1

mlab.figure(bgcolor=(0, 0, 0))
mlab.quiver3d(mask * ex[ind], mask * ey[ind], mask * ez[ind])
mlab.title("Masked E Field")
mlab.colorbar()
mlab.outline()

# Masked Vector Fields
mask = np.zeros((points_x, points_y, points_z))
mask[10:, 18:31, 18:32] = 1

mlab.figure(bgcolor=(0, 0, 0))
mlab.quiver3d(mask * hx[ind], mask * hy[ind], mask * hz[ind])
mlab.title("Masked H Field")
mlab.colorbar()
mlab.outline()
mlab.show()

"""
# Animate Slice
mlab.figure(bgcolor=(0,0,0))
s = mlab.imshow(ez[0,:,:,25], vmax=0.05, vmin=-0.05)
mlab.colorbar()
@mlab.animate(delay=30)
def anim():
    i = 0
    while True:
        i%=steps
        s.mlab_source.scalars = (ez[i,:,:,25])
        i+=1
        yield
anim()
mlab.show()
"""
