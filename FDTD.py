import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
from FDTD import fdtd

# Settings
f = 300e6

dr = 6/f
dt = 0.03/f
dx = 2*dt
steps = 1024;
points_x = 50;
points_y = 50;
points_z = 50;

t = np.linspace(0,(steps-1)*dt,steps)
s = (1-2*(np.pi*f*(t-dr))**2)*np.exp(-(np.pi*f*(t-dr))**2)

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
fdtd(dt, dx, steps, points_x, points_y, points_z, ex, ey, ez, hx, hy, hz, mu, ep, co, s)

plt.figure()
plt.plot(s)

plt.figure()
plt.plot(np.abs(np.fft.fft(s)))
plt.plot(np.abs(np.fft.fft(ez[:,15,25,25])))

plt.figure()
plt.plot(s)
plt.plot(ez[:,15,25,25])

H = (np.fft.fft(ez[:,15,25,25])/np.fft.fft(s))[0:int(steps/2)]
plt.figure()
plt.plot(abs(H))
plt.title("Amplitude Response")
plt.ylim([0,1])

plt.figure()
plt.plot(np.angle(H))
plt.title("Phase Response")
plt.show()

# Plot 3d Fields
'''
mlab.figure(bgcolor=(0,0,0))
sx = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="x_axes")
sy = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="y_axes")
sz = mlab.volume_slice(np.abs(ez[steps-1]), vmax=0.2, plane_orientation="z_axes")
mlab.colorbar()

# Animate Slice
mlab.figure(bgcolor=(0,0,0))
s = mlab.imshow(np.abs(ez[0,:,:,20]), vmax=0.05)
mlab.colorbar()
@mlab.animate(delay=30)
def anim():
    i = 0
    while True:
        i%=steps
        s.mlab_source.scalars = (ez[i,:,:,20])
        i+=1
        yield
anim()
mlab.show()
'''
