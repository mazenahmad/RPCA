# Relative Principal Components Analysis (RPCA) - Python Implementation

## Overview

This is a Python implementation of the RPCA (Relative Principal Components Analysis) tools originally written in C. The tool performs analysis of conformational changes upon biomolecular interactions.

## Installation

### Prerequisites

The Python implementation requires the following packages:

- Python 3.6 or higher
- NumPy
- SciPy
- MDAnalysis
- Matplotlib

### Installation Steps

1. Create a virtual environment (recommended):

```bash
python -m venv rpca_env
source rpca_env/bin/activate  # On Windows: rpca_env\Scripts\activate
```

2. Install the required packages:

```bash
pip install numpy scipy MDAnalysis matplotlib
```

3. Download the RPCA Python implementation:

```bash
git clone https://github.com/your_username/rpca-python.git  # Replace with actual repository
cd rpca-python
```

## Usage

The RPCA tool can be used as a standalone command-line tool or imported as a module in your Python scripts.

### Command-line Usage

Basic usage:

```bash
python rpca.py -fa trajectory_a.xtc -fb trajectory_b.xtc -sa structure_a.pdb -sb structure_b.pdb
```

Advanced options:

```bash
python rpca.py -fa trajectory_a.xtc -fb trajectory_b.xtc -sa structure_a.pdb -sb structure_b.pdb \
  -sel "protein and name CA" \
  -geig geigen.npy -dkl d_kl.png \
  -proj_a projection_a.npy -proj_b projection_b.npy \
  -res residue_interactions.npy -respdb contributions.pdb \
  -bt 0 -et 1000 -first 0 -last 10 -verbose 1
```

### Usage as a Python Module

```python
from rpca import RPCAAnalysis

# Create RPCA object
rpca = RPCAAnalysis()

# Read trajectories
universe_a = rpca.read_trajectory("trajectory_a.xtc", "structure_a.pdb")
universe_b = rpca.read_trajectory("trajectory_b.xtc", "structure_b.pdb")

# Perform GPA to find average structures
avg_coords_a = rpca.perform_gpa(universe_a, selection="protein and name CA")
avg_coords_b = rpca.perform_gpa(universe_b, selection="protein and name CA")

# Compute covariance matrices
cov_a, mean_a = rpca.compute_covariance_matrix(universe_a, selection="protein and name CA", 
                                             average_coords=avg_coords_a)
cov_b, mean_b = rpca.compute_covariance_matrix(universe_b, selection="protein and name CA", 
                                             average_coords=avg_coords_b)

# Perform simultaneous diagonalization
rpca_results = rpca.simultaneous_diagonalization(cov_a, cov_b, mean_a, mean_b)

# Plot KL divergence
rpca.plot_kl_divergence(rpca_results['kl'], rpca_results['kl_m'], 
                      rpca_results['acc_kl'], "kl_divergence.png")

# Project trajectories
proj_a = rpca.project_trajectory(universe_a, "protein and name CA", 
                               rpca_results['gevec'], mean_b,
                               first_vec=0, last_vec=10)

proj_b = rpca.project_trajectory(universe_b, "protein and name CA",
                               rpca_results['gevec'], mean_b,
                               first_vec=0, last_vec=10)

# Plot projections
rpca.plot_projections(proj_a, proj_b, first_vec=0, filename="projections.png")
```

## Output Files

- **geigen.npy**: Contains the generalized eigenpairs and related data
- **d_kl.png**: Plot of KL divergence components
- **projection_a.npy/projection_b.npy**: Projections of trajectories onto generalized eigenvectors
- **residue_interactions.npy**: Residue-residue interaction matrix
- **contributions.pdb**: PDB file with residue contributions as B-factors

## References

- Ahmad, M., Helms, V., Kalinina, O. V. & Lengauer, T. Relative Principal Components Analysis: Application to Analyzing Biomolecular Conformational Changes. J. Chem. Theory Comput. 15, 2166–2178 (2019).
- Ahmad, M., Helms, V., Kalinina, O. V. & Lengauer, T. Elucidating the energetic contributions to the binding free energy. J. Chem. Phys. 146, 014105 (2017).
- Ahmad, M., Helms, V., Kalinina, O. V. & Lengauer, T. The Role of Conformational Changes in Molecular Recognition. J. Phys. Chem. B 120, 2138–2144 (2016).
