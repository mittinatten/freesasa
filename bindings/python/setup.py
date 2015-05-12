from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("*", ["*.pyx"],
              include_dirs = ["../../src/"],
              libraries = ["libfreesasa"],
              library_dirs = ["../../src/"]
          )
]

setup(
    ext_modules = cythonize(extensions)
)
