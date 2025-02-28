# RPCA-Python

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/downloads/)

A high-performance Python implementation of Relative Principal Components Analysis (RPCA) for analyzing conformational changes in biomolecular systems. This package provides tools to identify and characterize the most significant differences between two molecular states (e.g., bound/unbound, wild-type/mutant).

## Features

- **Simultaneous diagonalization** of covariance matrices from two molecular states
- **Kullback-Leibler divergence** analysis to quantify conformational differences
- **Projection** of trajectories onto generalized eigenvectors
- **Residue contribution analysis** to identify conformational hotspots
- **Highly optimized performance** using parallelization, vectorization, and JIT compilation
- **Memory-efficient operation** for large molecular systems
- **User-friendly interface** with both command-line and Python API options

## Installation

### Using pip

```bash
pip install rpca-python
```

### From source

```bash
git clone https://github.com/username/rpca-python.git
cd rpca-python
pip install -e .
```

## Dependencies

- Python 3.7+
- NumPy
- SciPy
- MDAnalysis
- Matplotlib
- Numba

## Quick Start

### Command Line Interface

```bash
python -m rpca -fa trajectory_a.xtc -fb trajectory_b.xtc -sa structure_a.pdb -sb structure_b.pdb
```

### Python API

```python
from rpca import RPCAAnalysis

# Create RPCA object
rpca = RPCAAnalysis()

# Read trajectories
universe_a = rpca.read_trajectory("trajectory_a.xtc", "structure_a.pdb")
universe_b = rpca.read_trajectory("trajectory_b.xtc", "structure_b.pdb")

# Perform analysis
avg_coords_a = rpca.perform_gpa(universe_a, selection="protein and name CA")
avg_coords_b = rpca.perform_gpa(universe_b, selection="protein and name CA")
cov_a, mean_a = rpca.compute_covariance_matrix(universe_a, selection="protein and name CA")
cov_b, mean_b = rpca.compute_covariance_matrix(universe_b, selection="protein and name CA")
rpca_results = rpca.simultaneous_diagonalization(cov_a, cov_b, mean_a, mean_b)

# Visualize results
rpca.plot_kl_divergence(rpca_results['kl'], rpca_results['kl_m'], 
                       rpca_results['acc_kl'], "kl_divergence.png")
```

## Performance Optimization

RPCA-Python is optimized for high performance with large molecular systems:

- **Parallel processing** utilizing multiple CPU cores
- **Numba JIT compilation** for computational hotspots
- **Vectorized operations** leveraging NumPy's optimized array handling
- **Batched processing** for memory efficiency with large trajectories
- **Single precision option** for faster computation with minimal accuracy loss

Example with performance options:

```bash
python -m rpca -fa traj_a.xtc -fb traj_b.xtc -sa struct_a.pdb -sb struct_b.pdb \
  -sel "protein and name CA" \
  -n_jobs 8 -batch 1000 -fp32
```

## Documentation

Complete documentation is available at [https://rpca-python.readthedocs.io/](https://rpca-python.readthedocs.io/)

## Examples

See the `examples` directory for detailed examples of different analysis scenarios.

## Citation

If you use RPCA-Python in your research, please cite:

```
Ahmad, M., Helms, V., Kalinina, O. V. & Lengauer, T. Relative Principal Components Analysis:
Application to Analyzing Biomolecular Conformational Changes. J. Chem. Theory Comput. 15, 2166–2178 (2019).
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The original C implementation authors for the RPCA algorithm
- MDAnalysis team for the trajectory analysis framework
- NumPy, SciPy, and Numba developers for high-performance scientific computing tools
