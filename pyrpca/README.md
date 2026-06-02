# rpca (Python package) — proof of concept

A Python/NumPy package that reuses the **original C numerical core** of this
project for speed, with **no dependency on GROMACS or MDAnalysis**.

This is the first slice of a larger plan (see *Roadmap* below). It wires up the
highest-value, lowest-risk piece — the simultaneous diagonalization + relative
PCA Kullback–Leibler ranking — and proves the whole approach end to end.

## Why this works

The numerical core (`src/linear_algebra.c`) is pure C + LAPACK/BLAS: it does not
include a single GROMACS header. The GROMACS coupling in this project lives only
in (a) the four `main()` drivers (`gmx_*.c`, `SIMDAIG.c`) and (b) typedefs/PBC
helpers pulled in through some headers — none of which the linear-algebra core
touches. So the core can be compiled and called directly from Python.

```
Python (numpy) ──▶ Cython (_sdiag.pyx) ──▶ rpca_api.c ──▶ Simultaneous_Diagonalization()
                                                            (src/linear_algebra.c, LAPACK)
```

The built extension links **only** `libopenblas`, `libm`, `libc`, `libgfortran`
— verified with `ldd`. No GROMACS, no MDAnalysis.

## Layout

```
pyrpca/
├── pyproject.toml          build metadata (setuptools + Cython + numpy)
├── setup.py                compiles the extension against ../src/linear_algebra.c
├── src/rpca/
│   ├── __init__.py
│   ├── _sdiag.pyx          Cython binding (numpy <-> C, zero-copy, nogil)
│   └── _csrc/
│       ├── rpca_api.h      thin buffer-in/buffer-out C API
│       └── rpca_api.c      calls the core + KL scoring/ordering glue
└── tests/test_sdiag.py     validates G^T A G = I, G^T B G = diag, KL, vs SciPy
```

## Build & test

Requires a C compiler, LAPACK/BLAS dev libraries, and Python build deps:

```bash
# system (Debian/Ubuntu)
sudo apt-get install -y libopenblas-dev liblapack-dev
# python
pip install numpy cython pytest scipy

cd pyrpca
python setup.py build_ext --inplace
PYTHONPATH=src python -m pytest tests/ -v
```

## Usage

```python
import numpy as np
from rpca import sdiag

res = sdiag(cov_a, cov_b, mean_a, mean_b)   # symmetric (n, n) covariances + means
res["geigval"]   # generalized eigenvalues (var_B / var_A per mode)
res["gevec"]     # (n, rank) eigenvectors; column k is mode k
res["kl"]        # per-mode KL divergence, descending
res["kl_m"]      # mean-shift contribution to KL
res["acc_kl"]    # cumulative % of total KL
res["rank"]      # effective rank
```

## Roadmap

1. **[done]** De-GROMACS the core for the sdiag path; Cython binding; tests.
2. Stream covariance + GPA over `.xtc` directly in C (reuse `covariance.c`,
   `fitting.c` behind a small `compat.h` that replaces the GROMACS typedefs and
   makes PBC optional).
3. Restore covariance-weighted Mahalanobis fitting (`CWfitting.c` + the BFGS
   minimizer) — a step the earlier numba port dropped.
4. Pure-Python `io` (gro/pdb parsers, PDB B-factor writer) + thin `xtc` wrapper
   over the bundled `xdrfile`, plus a minimal atom-selection layer.
5. High-level `RPCA` pipeline class + `argparse` CLI mirroring the old tools.
6. Switch packaging to `scikit-build-core`/CMake (or `meson-python`) and vendor
   the C core into the package for distributable wheels.

## Notes

- For the PoC the build references `../src/linear_algebra.c` so there is a single
  source of truth. Productionising (step 6) vendors a copy into the package.
- The KL formula matches the published RPCA definition (Ahmad et al., *JCTC*
  2019) and the existing Python reference implementation.
