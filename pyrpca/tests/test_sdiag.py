"""Tests for rpca.sdiag.

The defining property of simultaneous diagonalization gives a self-contained
ground truth (no second library required): for the generalized eigenbasis G,

    G^T A G = I        and        G^T B G = diag(geigval).

We also check the KL bookkeeping and the descending ordering, and - when SciPy
is available - cross-check the eigenvalues against scipy.linalg.eigh.
"""
import numpy as np
import pytest

from rpca import sdiag


def _random_spd(n, rng, scale=1.0):
    """A random symmetric positive-definite n x n matrix."""
    M = rng.standard_normal((n, n))
    return scale * (M @ M.T) + n * np.eye(n)


def test_simultaneous_diagonalization_property():
    rng = np.random.default_rng(0)
    n = 15
    A = _random_spd(n, rng)
    B = _random_spd(n, rng, scale=2.0)

    res = sdiag(A, B)
    G = res["gevec"]
    r = res["rank"]

    assert G.shape == (n, r)
    assert r == n  # both inputs full rank

    GtAG = G.T @ A @ G
    GtBG = G.T @ B @ G

    np.testing.assert_allclose(GtAG, np.eye(r), atol=1e-8)
    # off-diagonal of G^T B G must vanish; diagonal equals the eigenvalues
    np.testing.assert_allclose(GtBG - np.diag(np.diag(GtBG)), 0.0, atol=1e-8)
    np.testing.assert_allclose(np.diag(GtBG), res["geigval"], atol=1e-8)


def test_kl_formula_and_ordering():
    rng = np.random.default_rng(1)
    n = 12
    A = _random_spd(n, rng)
    B = _random_spd(n, rng, scale=1.5)
    mean_a = rng.standard_normal(n)
    mean_b = rng.standard_normal(n)

    res = sdiag(A, B, mean_a, mean_b)
    G, gv = res["gevec"], res["geigval"]

    # KL must be sorted descending
    assert np.all(np.diff(res["kl"]) <= 1e-12)

    # reproduce the published per-mode KL from the returned pieces
    proj = G.T @ (mean_b - mean_a)
    kl_m_expected = 0.5 * proj ** 2
    kl_expected = 0.5 * (gv - np.log(gv) - 1.0) + kl_m_expected

    np.testing.assert_allclose(res["kl_m"], kl_m_expected, atol=1e-9)
    np.testing.assert_allclose(res["kl"], kl_expected, atol=1e-9)

    # accumulated KL ends at 100%
    np.testing.assert_allclose(res["acc_kl"][-1], 100.0, atol=1e-6)
    assert res["sum_kl"] == pytest.approx(res["kl"].sum())


def test_eigenvalues_match_scipy():
    scipy_linalg = pytest.importorskip("scipy.linalg")
    rng = np.random.default_rng(2)
    n = 10
    A = _random_spd(n, rng)
    B = _random_spd(n, rng, scale=3.0)

    res = sdiag(A, B)
    # generalized problem B v = lambda A v -> eigenvalues are var_B/var_A
    ref = np.sort(scipy_linalg.eigh(B, A, eigvals_only=True))
    np.testing.assert_allclose(np.sort(res["geigval"]), ref, rtol=1e-6, atol=1e-8)
