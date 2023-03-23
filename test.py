import mayavi.mlab
import numpy

data = (50, 50, 50)
data = numpy.zeros(data)
data[1:48, 20:30, 28:30] = 1
data[1:48, 20:30, 20:22] = 1

data[1, 20:30, 22:28] = 2

data[1:48, 20:30, 22:28] = 3

xx, yy, zz = numpy.where(data == 1)

mayavi.mlab.points3d(xx, yy, zz,
                     mode="cube",
                     color=(0, 1, 0),
                     scale_factor=1)

mayavi.mlab.outline()
mayavi.mlab.axes()

xx, yy, zz = numpy.where(data == 2)

mayavi.mlab.points3d(xx, yy, zz,
                     mode="cube",
                     color=(1, 0, 0),
                     scale_factor=1)

xx, yy, zz = numpy.where(data == 3)

mayavi.mlab.points3d(xx, yy, zz,
                     mode="cube",
                     color=(0, 0, 1),
                     scale_factor=1)

mayavi.mlab.show()