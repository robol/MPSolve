from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("mpsolve",
                             sources=["mpsolve.pyx"],
                             libraries=["mps", "gmp"],
                             include_dirs=[os.environ['SAGE_LOCAL']+'/include/mps'])]
)
