"""rpca - Relative Principal Components Analysis backed by the original C core.

This package reuses the proven C numerical routines (simultaneous
diagonalization, whitening, LAPACK eigendecomposition) from the parent project
and exposes them to Python/NumPy without any GROMACS or MDAnalysis dependency.

Proof-of-concept stage: only the simultaneous-diagonalization step is wired up.
"""

from ._sdiag import sdiag

__all__ = ["sdiag"]
__version__ = "0.0.1"
