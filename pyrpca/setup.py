"""Build script for the rpca extension.

The extension links the package's thin C API (rpca_api.c) against the existing
numerical core in the parent repository (src/linear_algebra.c), so there is a
single source of truth for the algorithms. Productionising this would vendor a
copy of the core into the package; for the proof of concept we reference it in
place.
"""
import os

import numpy as np
from setuptools import Extension, setup
from Cython.Build import cythonize

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)  # parent RPCA repository

CORE_SRC = os.path.join(REPO, "src", "linear_algebra.c")
CORE_INC = os.path.join(REPO, "include")
PKG_CSRC = os.path.join(HERE, "src", "rpca", "_csrc")

# OpenBLAS bundles both CBLAS and the LAPACK routine (dsyevr_) the core needs;
# liblapack is listed as a fallback provider of the Fortran symbols.
blas_libs = ["openblas", "lapack", "m"]

extensions = [
    Extension(
        name="rpca._sdiag",
        sources=[
            os.path.join(HERE, "src", "rpca", "_sdiag.pyx"),
            os.path.join(PKG_CSRC, "rpca_api.c"),
            CORE_SRC,
        ],
        include_dirs=[CORE_INC, PKG_CSRC, np.get_include()],
        libraries=blas_libs,
        extra_compile_args=["-O3"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]

setup(
    name="rpca",
    version="0.0.1",
    package_dir={"": "src"},
    packages=["rpca"],
    ext_modules=cythonize(extensions, language_level=3),
    zip_safe=False,
)
