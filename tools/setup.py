from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import sys

include_gsl_dir = "/usr/local/opt/gsl/include/"
lib_gsl_dir = "/usr/local/opt/"

name_mdl = "flar"

extensions= [
    Extension(name_mdl,
        sources=["Wrapper.pyx", "splinecoeffs.c", "struct.c", "Faddeeva.c", "fresnel.c", "likelihoodKurz.c"],
        include_dirs=[numpy.get_include(), include_gsl_dir],
        library_dirs=[lib_gsl_dir],
        libraries=["gsl", "gslcblas"])
]

from Cython.Build import cythonize
extensions = cythonize(extensions)

NUMPY_DEP = 'numpy>=1.11'

SETUP_REQUIRES = [NUMPY_DEP]

setup(name=name_mdl,
    ext_modules=extensions,
    install_requires=SETUP_REQUIRES,
    py_modules = ['pyFDresponse']
    )


#    cmdclass = {'build_ext': build_ext})
