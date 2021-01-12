import distutils.core
import Cython.Build
import numpy
distutils.core.setup(
    ext_modules = Cython.Build.cythonize("test_0.pyx"),
    include_dirs=[numpy.get_include()])