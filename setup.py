from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = "PyNBody",
        ext_modules = cythonize(('PyNBody.pyx'))
        )
