from setuptools import setup
from Cython.Build import cythonize

import numpy

setup(
    name='FDTD',
    ext_modules=cythonize("FDTD.pyx"),
    zip_safe=False,
    include_dirs=[numpy.get_include()]
)