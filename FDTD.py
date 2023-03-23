import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
from FDTD import fdtd

# Settings
f = 300e6

dr = 6/f
dt = 0.03/f
dx = 2*dt
steps = 512;
points_x = 50;
points_y = 50;
points_z = 50;

t = np.linspace(0,(steps-1)*dt,steps)
s = (1-2*(np.pi*f*(t-dr))**2)*np.exp(-(np.pi*f*(t-dr))**2) # Source waveform

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

# 2 Parallel Conductors
ep[1:48, 20:30, 28:30] = 1e9
ep[1:48, 20:30, 20:22] = 1e9
co[1:48, 20:30, 20:22] = 1e9
co[1:48, 20:30, 20:22] = 1e9

# Source locations
source_pos = np.zeros((points_x, points_y, points_z))
source_pos[1, 20:30, 22:28] = 1
sx, sy, sz = np.where(source_pos == 1)

# Pass to C code for main processing loop
fdtd(dt, dx, steps, points_x, points_y, points_z, ex, ey, ez, hx, hy, hz, mu, ep, co, s, sx, sy, sz)

v1 = s
v2 = ez[:,48,25,25]

plt.figure()
plt.plot(s)

plt.figure()
plt.plot(np.abs(np.fft.fft(v1)))
plt.plot(np.abs(np.fft.fft(v2)))

plt.figure()
plt.plot(v1)
plt.plot(v2)

H = (np.fft.fft(v2)/np.fft.fft(v1))[0:int(steps/2)]
plt.figure()
plt.plot(abs(H))
plt.title("Amplitude Response")
plt.ylim([0,10])

plt.figure()
plt.plot(np.angle(H))
plt.title("Phase Response")
plt.show()

# Plot 3d Fields

mlab.figure(bgcolor=(0,0,0))
sx = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="x_axes")
sy = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="y_axes")
sz = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="z_axes")
mlab.colorbar()

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
