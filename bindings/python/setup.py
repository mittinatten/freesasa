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
    name='FreeSASA',
    description='Calculate solvent accessible surface areas of proteins',
    version='0.5.0', # how to make this use the version of the rest of the library
    author='Simon Mitternacht',
    url='http://mittinatten.github.io/freesasa/',
    license='GPL-v3',
    ext_modules = cythonize(extensions)
)
