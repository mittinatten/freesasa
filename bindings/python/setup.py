from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("*", ["*.pyx"],
              include_dirs = ["../../src"],
              language='c',
              extra_objects = ["../../src/libfreesasa.a"],
              extra_compile_args = ["-w"] 
              )
]

setup(
    ext_modules = cythonize(extensions)
)
