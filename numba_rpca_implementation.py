"""
Relative Principal Components Analysis (RPCA) for Biomolecular Conformational Changes

This is a Python implementation of the RPCA tools originally written in C.
It performs RPCA analysis of conformational changes upon biomolecular interactions.

Based on research by:
Ahmad, M., Helms, V., Kalinina, O. V. & Lengauer, T. Relative Principal Components Analysis:
Application to Analyzing Biomolecular Conformational Changes. J. Chem. Theory Comput. 15, 2166–2178 (2019).
"""

import numpy as np
import os
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisBase
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import time
import multiprocessing
from numba import jit, prange, float64, int64
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


class RPCAAnalysis:
    """Main class for Relative Principal Components Analysis"""
    
    def __init__(self):
        """Initialize RPCA analysis object"""
        self.precision = 1000.0  # Default precision for coordinates
        self.verbose = 1         # Default verbosity level
        
    def read_trajectory(self, trajectory_file, topology_file=None, start_time=0, end_time=-1):
        """Read molecular dynamics trajectory
        
        Args:
            trajectory_file: Path to trajectory file (xtc, dcd, etc.)
            topology_file: Path to topology file (pdb, gro, etc.)
            start_time: Start time (ps) for analysis
            end_time: End time (ps) for analysis, -1 means until the end
            
        Returns:
            MDAnalysis Universe object
        """
        try:
            if topology_file:
                universe = mda.Universe(topology_file, trajectory_file)
            else:
                universe = mda.Universe(trajectory_file)
                
            # Set trajectory start/end if specified
            if start_time > 0 or end_time > 0:
                # Find frames corresponding to time range
                times = universe.trajectory.time
                start_frame = 0
                end_frame = len(universe.trajectory) - 1
                
                if start_time > 0:
                    for i, t in enumerate(times):
                        if t >= start_time:
                            start_frame = i
                            break
                            
                if end_time > 0:
                    for i, t in enumerate(times):
                        if t > end_time:
                            end_frame = i - 1
                            break
                
                print(f"Using trajectory frames from {start_frame} to {end_frame}")
                universe.trajectory[start_frame:end_frame]
            
            return universe
            
        except Exception as e:
            print(f"Error reading trajectory: {e}")
            return None
    
    @jit(nopython=True, parallel=True)
    def _accumulate_coords(self, mobile_coords_array, reference_coords, n_frames, n_atoms):
        """Numba-accelerated function to accumulate coordinates for GPA"""
        avg_coords = np.zeros((n_atoms, 3), dtype=np.float64)
        
        for i in prange(n_frames):
            mobile_coords = mobile_coords_array[i]
            
            # Center coordinates
            mobile_mean = np.mean(mobile_coords, axis=0)
            ref_mean = np.mean(reference_coords, axis=0)
            
            mobile_centered = mobile_coords - mobile_mean
            ref_centered = reference_coords - ref_mean
            
            # Calculate correlation matrix
            correlation_matrix = np.dot(mobile_centered.T, ref_centered)
            
            # SVD
            u, s, vh = np.linalg.svd(correlation_matrix)
            
            # Ensure proper rotation (not reflection)
            d = np.sign(np.linalg.det(np.dot(vh.T, u.T)))
            
            # Construct rotation matrix
            if d < 0:
                vh[-1] = -vh[-1]  # Flip last row of vh
            
            rotation = np.dot(u, vh)
            
            # Apply rotation
            rotated_coords = np.dot(mobile_centered, rotation)
            
            # Accumulate for average
            avg_coords += rotated_coords
            
        return avg_coords / n_frames
    
    def perform_gpa(self, universe, selection='protein', max_iterations=10, 
                    convergence=0.00001, ref_frame=0, n_jobs=None):
        """Perform Generalized Procrustes Analysis to find average structure
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            max_iterations: Maximum number of GPA iterations
            convergence: RMSD convergence criterion
            ref_frame: Initial reference frame number
            n_jobs: Number of parallel jobs to run (None uses all available cores)
            
        Returns:
            Average coordinates after GPA
        """
        print("Performing Generalized Procrustes Analysis (GPA)...")
        
        if n_jobs is None:
            n_jobs = multiprocessing.cpu_count()
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        
        # Initialize with the first frame as reference
        universe.trajectory[ref_frame]
        reference_coords = atoms.positions.copy()
        
        # Pre-load all coordinates into memory for faster processing
        # This trades memory for speed
        print("Pre-loading trajectory coordinates...")
        coords_array = []
        for ts in universe.trajectory:
            coords_array.append(atoms.positions.copy())
        
        coords_array = np.array(coords_array)
        n_frames = len(coords_array)
        print(f"Loaded {n_frames} frames")
        
        # Iterative GPA using optimized parallel processing
        for iteration in range(max_iterations):
            print(f"GPA Iteration {iteration + 1}/{max_iterations}")
            
            # Process frames in parallel chunks
            chunk_size = max(1, n_frames // n_jobs)
            chunks = [coords_array[i:i+chunk_size] for i in range(0, n_frames, chunk_size)]
            
            avg_coords = np.zeros_like(reference_coords)
            
            with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                futures = []
                for chunk in chunks:
                    futures.append(
                        executor.submit(
                            self._accumulate_coords_wrapper,
                            chunk, reference_coords, len(chunk), n_atoms
                        )
                    )
                
                # Collect results
                chunk_results = [future.result() for future in futures]
                
                # Weight by chunk size and combine
                for i, result in enumerate(chunk_results):
                    weight = len(chunks[i]) / n_frames
                    avg_coords += result * weight
            
            # Calculate RMSD between old and new reference
            rmsd = np.sqrt(np.mean(np.sum((reference_coords - avg_coords)**2, axis=1)))
            print(f"  RMSD between iterations: {rmsd:.6f}")
            
            # Check convergence
            if rmsd < convergence and iteration > 0:
                print(f"GPA converged after {iteration+1} iterations")
                break
                
            # Update reference for next iteration
            reference_coords = avg_coords.copy()
            
        return avg_coords
    
    def _accumulate_coords_wrapper(self, coords_chunk, reference_coords, chunk_size, n_atoms):
        """Wrapper for the numba-accelerated function to handle Python objects"""
        # Convert to numpy arrays of the right shape and type for numba
        coords_array_np = np.array(coords_chunk, dtype=np.float64)
        reference_coords_np = np.array(reference_coords, dtype=np.float64)
        
        return self._accumulate_coords(coords_array_np, reference_coords_np, chunk_size, n_atoms)
    
    @staticmethod
    @jit(nopython=True)
    def _calculate_rotation_matrix(mobile, reference):
        """Calculate optimal rotation matrix using SVD (Numba-accelerated)
        
        Args:
            mobile: Mobile coordinates (centered)
            reference: Reference coordinates (centered)
            
        Returns:
            Rotation matrix and RMSD
        """
        # Compute correlation matrix
        correlation_matrix = np.dot(mobile.T, reference)
        
        # Singular value decomposition
        U, S, Vt = np.linalg.svd(correlation_matrix)
        
        # Ensure proper rotation (not reflection)
        d = np.sign(np.linalg.det(np.dot(Vt.T, U.T)))
        
        # Construct rotation matrix
        if d < 0:
            S = np.diag([1, 1, d])
            rotation = np.dot(U, np.dot(S, Vt))
        else:
            rotation = np.dot(U, Vt)
        
        # Calculate RMSD
        rotated = np.dot(mobile, rotation)
        rmsd = np.sqrt(np.mean(np.sum((rotated - reference)**2, axis=1)))
        
        return rotation, rmsd
    
    class FastCovarianceAnalysis(AnalysisBase):
        """Optimized covariance calculation using MDAnalysis analysis framework"""
        def __init__(self, atomgroup, reference_coords=None, **kwargs):
            super(FastCovarianceAnalysis, self).__init__(atomgroup.universe.trajectory, **kwargs)
            self.atomgroup = atomgroup
            self.n_atoms = len(atomgroup)
            self.n_dims = self.n_atoms * 3
            self.reference_coords = reference_coords
            self.mean_coords = None
            self.covariance_matrix = None
            
        def _prepare(self):
            # Initialize arrays
            self.covariance_matrix = np.zeros((self.n_dims, self.n_dims), dtype=np.float64)
            
            if self.reference_coords is None:
                # Need to compute mean coordinates first
                self.mean_coords = np.zeros((self.n_atoms, 3), dtype=np.float64)
                self.collect_coordinates = True
                self.coordinates = []
            else:
                # Use provided reference coordinates
                self.mean_coords = self.reference_coords
                self.collect_coordinates = False
            
        def _single_frame(self):
            if self.collect_coordinates:
                self.coordinates.append(self.atomgroup.positions.copy())
            else:
                # Direct covariance accumulation
                coords = self.atomgroup.positions.flatten()
                deviation = coords - self.mean_coords.flatten()
                self.covariance_matrix += np.outer(deviation, deviation)
        
        def _conclude(self):
            if self.collect_coordinates:
                # Calculate mean if not provided
                self.coordinates = np.array(self.coordinates)
                self.mean_coords = np.mean(self.coordinates, axis=0)
                
                # Calculate covariance using vectorized operations
                flattened_coords = self.coordinates.reshape(len(self.coordinates), -1)
                mean_flattened = self.mean_coords.flatten()
                
                # Use np.einsum for efficient matrix multiplication
                # This computes the sum of outer products efficiently
                self.covariance_matrix = np.einsum('ij,ik->jk', 
                                                 flattened_coords - mean_flattened, 
                                                 flattened_coords - mean_flattened) / len(self.coordinates)
            else:
                # Normalize the accumulated covariance
                self.covariance_matrix /= self.n_frames
    
    def compute_covariance_matrix(self, universe, selection='protein', average_coords=None, n_jobs=None):
        """Compute covariance matrix from trajectory (optimized)
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            average_coords: Pre-computed average coordinates (optional)
            n_jobs: Number of parallel jobs for calculations
            
        Returns:
            Covariance matrix and mean coordinates
        """
        print("Computing covariance matrix...")
        
        # Set number of threads for BLAS operations
        if n_jobs is None:
            n_jobs = multiprocessing.cpu_count()
            
        os.environ['OMP_NUM_THREADS'] = str(n_jobs)
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        
        # Initialize and run the optimized analysis
        cov_analysis = self.FastCovarianceAnalysis(atoms, reference_coords=average_coords, n_jobs=n_jobs)
        cov_analysis.run()
        
        return cov_analysis.covariance_matrix, cov_analysis.mean_coords.flatten()
    
    @staticmethod
    @jit(nopython=True)
    def _compute_kl_divergences(gevec, geigval, mean_diff, rank):
        """Compute KL divergences with Numba acceleration"""
        kl = np.zeros(rank, dtype=np.float64)
        kl_m = np.zeros(rank, dtype=np.float64)
        
        for i in range(rank):
            # Compute mean contribution
            mean_proj = np.dot(gevec[:, i], mean_diff)
            kl_m[i] = 0.5 * mean_proj**2
            
            # Total KL divergence
            kl[i] = 0.5 * (-np.log(geigval[i]) + geigval[i] - 1) + kl_m[i]
        
        return kl, kl_m
    
    def simultaneous_diagonalization(self, cov_a, cov_b, mean_a, mean_b):
        """Perform simultaneous diagonalization of two covariance matrices (optimized)
        
        Args:
            cov_a: Covariance matrix of state A
            cov_b: Covariance matrix of state B
            mean_a: Mean coordinates of state A
            mean_b: Mean coordinates of state B
            
        Returns:
            Dictionary containing:
                geigval: Generalized eigenvalues
                gevec: Generalized eigenvectors
                kl: Kullback-Leibler divergences
                kl_m: KL divergences due to mean shifts
                rank: Rank of the decomposition
        """
        print("Performing simultaneous diagonalization...")
        
        start_time = time.time()
        
        # Dimension of the problem
        n_dims = cov_a.shape[0]
        
        # Use scipy's optimized eigenvalue solver for symmetric matrices
        # This is faster than numpy.linalg.eigh for large matrices
        evals_a, evecs_a = linalg.eigh(cov_a, driver='evd')
        
        # Filter out near-zero eigenvalues
        tol = np.sqrt(np.finfo(float).eps) * np.max(evals_a)
        mask = evals_a > tol
        rank = np.sum(mask)
        
        print(f"Matrix A has rank {rank} out of {n_dims} dimensions")
        
        # Use direct indexing and broadcasting for better performance
        evals_masked = evals_a[mask]
        evecs_masked = evecs_a[:, mask]
        
        # Compute whitening transformation using broadcasting
        W = evecs_masked * (1.0 / np.sqrt(evals_masked))
        
        # Transform B with whitening transformation efficiently
        # Use matrix multiplication optimizations
        WtB = np.dot(W.T, cov_b)
        transformed_b = np.dot(WtB, W)
        
        # Make sure transformed_b is symmetric for better eigenvalue computation
        transformed_b = 0.5 * (transformed_b + transformed_b.T)
        
        # Eigendecomposition of transformed B
        geigval, gevec_temp = linalg.eigh(transformed_b)
        
        # Full eigenvectors in original space
        gevec = np.dot(W, gevec_temp)
        
        # Mean difference vector
        mean_diff = mean_b - mean_a
        
        # Compute KL divergences using Numba acceleration
        kl, kl_m = self._compute_kl_divergences(gevec, geigval, mean_diff, rank)
        
        # Sort by KL divergence efficiently
        sort_idx = np.argsort(-kl)
        
        # Use vectorized operations for sorting
        kl = kl[sort_idx]
        kl_m = kl_m[sort_idx]
        geigval = geigval[sort_idx]
        gevec = gevec[:, sort_idx]
        
        # Accumulated KL divergence (percentage)
        kl_sum = np.sum(kl)
        acc_kl = np.cumsum(kl) / kl_sum * 100
        
        elapsed_time = time.time() - start_time
        print(f"Simultaneous diagonalization completed in {elapsed_time:.2f} seconds")
        
        return {
            'geigval': geigval,
            'gevec': gevec,
            'kl': kl,
            'kl_m': kl_m,
            'acc_kl': acc_kl,
            'rank': rank,
            'sum_kl': kl_sum,
            'sum_kl_m': np.sum(kl_m)
        }
    
    @staticmethod
    @jit(nopython=True, parallel=True)
    def _project_coordinates(coords_array, mean_coords, gevec, first_vec, n_vecs):
        """Optimized projection of coordinates onto eigenvectors"""
        n_frames = coords_array.shape[0]
        projections = np.zeros((n_frames, n_vecs), dtype=np.float64)
        
        for i in prange(n_frames):
            coords = coords_array[i].flatten()
            deviation = coords - mean_coords
            
            for j in range(n_vecs):
                vec_idx = first_vec + j
                projections[i, j] = np.dot(deviation, gevec[:, vec_idx])
                
        return projections
    
    def project_trajectory(self, universe, selection, gevec, mean_coords, first_vec=0, last_vec=None, batch_size=1000, n_jobs=None):
        """Project trajectory onto generalized eigenvectors (optimized)
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            gevec: Generalized eigenvectors 
            mean_coords: Mean coordinates
            first_vec: First eigenvector to include (0-based)
            last_vec: Last eigenvector to include (0-based), None means all
            batch_size: Number of frames to process at once
            n_jobs: Number of parallel jobs to run (None uses all available cores)
            
        Returns:
            Projections array, shape (n_frames, n_vecs)
        """
        print("Projecting trajectory onto eigenvectors...")
        start_time = time.time()
        
        if n_jobs is None:
            n_jobs = multiprocessing.cpu_count()
            
        # Set environment variables for optimal performance
        os.environ['OMP_NUM_THREADS'] = str(n_jobs)
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        n_dims = n_atoms * 3
        
        # Set number of eigenvectors to use
        if last_vec is None:
            last_vec = gevec.shape[1] - 1
        
        n_vecs = last_vec - first_vec + 1
        gevec_subset = gevec[:, first_vec:last_vec+1]
        
        # Determine total number of frames
        n_frames = len(universe.trajectory)
        
        # Process trajectory in batches to minimize memory usage
        projections = np.zeros((n_frames, n_vecs), dtype=np.float64)
        
        # Process trajectory in batches
        for batch_start in range(0, n_frames, batch_size):
            batch_end = min(batch_start + batch_size, n_frames)
            batch_size_actual = batch_end - batch_start
            
            print(f"Processing frames {batch_start}-{batch_end-1} ({batch_size_actual} frames)")
            
            # Collect coordinates for this batch
            coords_array = np.zeros((batch_size_actual, n_atoms, 3), dtype=np.float64)
            
            for i, ts in enumerate(universe.trajectory[batch_start:batch_end]):
                coords_array[i] = atoms.positions
            
            # Project using matrix multiplication instead of loops for better BLAS acceleration
            # Reshape the coordinates array for efficient matrix operations
            flattened_coords = coords_array.reshape(batch_size_actual, -1)
            centered_coords = flattened_coords - mean_coords
            
            # Project onto eigenvectors using matrix multiplication
            batch_projections = np.dot(centered_coords, gevec_subset)
            
            # Store results
            projections[batch_start:batch_end] = batch_projections
        
        elapsed_time = time.time() - start_time
        print(f"Projection completed in {elapsed_time:.2f} seconds")
        
        return projections

    def compute_residue_indices(self, atoms, residues):
        """Compute residue atom indices efficiently"""
        residue_indices = []
        
        # Create a lookup dictionary for fast atom index retrieval
        atom_indices = {atom.index: i for i, atom in enumerate(atoms)}
        
        for residue in residues:
            # Get atoms in this residue that are also in the selection
            res_atoms = residue.atoms.intersection(atoms)
            # Get indices in the atoms selection
            idx = [atom_indices[atom.index] for atom in res_atoms]
            residue_indices.append(np.array(idx))
            
        return residue_indices
    
    @staticmethod
    @jit(nopython=True)
    def _compute_interaction_matrix(gevec_subset, kl_subset, residue_indices, n_res, n_atoms, n_vecs):
        """Numba-accelerated calculation of interaction matrix"""
        interaction_matrix = np.zeros((n_res, n_res), dtype=np.float64)
        atom_contrib = np.zeros(n_atoms, dtype=np.float64)
        
        # Pre-reshape gevec for better memory access patterns
        gevec_reshaped = gevec_subset.reshape(n_vecs, n_atoms, 3)
        
        # Loop over eigenvectors first for better cache locality
        for k in range(n_vecs):
            vec = gevec_reshaped[k]
            weight = kl_subset[k]
            
            # Update per-residue contributions
            for i in range(n_res):
                idx_i = residue_indices[i]
                
                for j in range(n_res):
                    idx_j = residue_indices[j]
                    
                    # Compute contribution between residue i and j
                    contrib = 0.0
                    for a_idx in range(len(idx_i)):
                        a = idx_i[a_idx]
                        for b_idx in range(len(idx_j)):
                            b = idx_j[b_idx]
                            contrib += np.dot(vec[a], vec[b]) * weight
                    
                    interaction_matrix[i, j] += contrib
                    
                    # Update per-atom contribution for diagonal elements
                    if i == j:
                        contrib_per_atom = contrib / len(idx_i)
                        for a_idx in range(len(idx_i)):
                            atom_contrib[idx_i[a_idx]] += contrib_per_atom
        
        return interaction_matrix, atom_contrib
    
    def compute_interaction_map(self, universe, selection, gevec, kl, first_vec=0, last_vec=None, n_jobs=None):
        """Compute per-residue contribution to conformational changes (optimized)
        
        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string
            gevec: Generalized eigenvectors
            kl: KL divergences
            first_vec: First eigenvector to include (0-based)
            last_vec: Last eigenvector to include (0-based), None means all
            n_jobs: Number of parallel jobs to run (None uses all available cores)
            
        Returns:
            Residue interaction matrix and per-atom contributions
        """
        print("Computing interaction map...")
        start_time = time.time()
        
        if n_jobs is None:
            n_jobs = multiprocessing.cpu_count()
            
        # Set environment variables for optimal performance
        os.environ['OMP_NUM_THREADS'] = str(n_jobs)
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        
        # Set number of eigenvectors to use
        if last_vec is None:
            last_vec = gevec.shape[1] - 1
            
        n_vecs = last_vec - first_vec + 1
        gevec_subset = gevec[:, first_vec:last_vec+1].copy()
        kl_subset = kl[first_vec:last_vec+1].copy()
        
        # Get unique residues
        residues = atoms.residues
        n_res = len(residues)
        
        # Pre-compute residue indices
        print("Pre-computing residue indices...")
        residue_indices = self.compute_residue_indices(atoms, residues)
        
        # Prepare numba-compatible format
        residue_indices_numba = []
        for indices in residue_indices:
            residue_indices_numba.append(np.array(indices, dtype=np.int64))
        
        # Compute interaction matrix
        print(f"Computing interactions for {n_vecs} eigenvectors and {n_res} residues...")
        interaction_matrix, atom_contrib = self._compute_interaction_matrix(
            gevec_subset.reshape(-1, n_atoms, 3), 
            kl_subset, 
            residue_indices_numba, 
            n_res, 
            n_atoms, 
            n_vecs
        )
        
        # Normalize atom contributions
        if np.max(atom_contrib) > 0:
            atom_contrib /= np.max(atom_contrib)
        
        elapsed_time = time.time() - start_time
        print(f"Interaction map computed in {elapsed_time:.2f} seconds")
            
        return interaction_matrix, atom_contrib
    
    def save_pdb_with_bfactors(self, universe, selection, bfactors, filename):
        """Save PDB with B-factors set to given values
        
        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string
            bfactors: B-factor values to assign
            filename: Output PDB filename
        """
        # Select atoms
        atoms = universe.select_atoms(selection)
        
        # Set B-factors
        atoms.tempfactors = bfactors
        
        # Write PDB
        atoms.write(filename)
        print(f"Wrote PDB with B-factors to {filename}")
    
    def plot_kl_divergence(self, kl, kl_m, acc_kl, filename=None):
        """Plot KL divergence components
        
        Args:
            kl: Total KL divergence for each component
            kl_m: Mean-shift contribution to KL
            acc_kl: Accumulated KL (percentage)
            filename: Output filename (optional)
        """
        plt.figure(figsize=(10, 6))
        
        x = np.arange(1, len(kl) + 1)
        
        plt.bar(x, kl, alpha=0.7, label='Total')
        plt.bar(x, kl_m, alpha=0.7, label='Due to mean change')
        
        plt.plot(x, acc_kl, 'r-', label='Accumulated (%)')
        
        plt.xlabel('Generalized Eigenvector')
        plt.ylabel('KL Divergence')
        plt.title('Generalized Eigenanalysis')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        if filename:
            plt.savefig(filename, dpi=300)
            print(f"KL divergence plot saved to {filename}")
        else:
            plt.show()
    
    def plot_projections(self, proj_a, proj_b, first_vec, filename=None):
        """Plot projections of trajectories onto eigenvectors
        
        Args:
            proj_a: Projections of state A
            proj_b: Projections of state B
            first_vec: First eigenvector index (for labeling)
            filename: Output filename (optional)
        """
        n_vecs = min(proj_a.shape[1], 4)  # Plot at most 4 eigenvectors
        
        plt.figure(figsize=(12, 8))
        
        for i in range(n_vecs):
            plt.subplot(2, 2, i+1)
            
            plt.hist(proj_a[:, i], bins=30, alpha=0.5, label='State A')
            plt.hist(proj_b[:, i], bins=30, alpha=0.5, label='State B')
            
            plt.xlabel(f'Projection on eigenvector {first_vec + i + 1}')
            plt.ylabel('Frequency')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300)
            print(f"Projection plot saved to {filename}")
        else:
            plt.show()


def main():
    """Main function to run RPCA analysis from command line"""
    parser = argparse.ArgumentParser(description='Relative Principal Components Analysis')
    
    # Input files
    parser.add_argument('-fa', '--traj_a', required=True, help='Trajectory file for state A')
    parser.add_argument('-fb', '--traj_b', required=True, help='Trajectory file for state B')
    parser.add_argument('-sa', '--top_a', required=True, help='Topology file for state A')
    parser.add_argument('-sb', '--top_b', required=True, help='Topology file for state B')
    
    # Selection options
    parser.add_argument('-sel', '--selection', default='protein and name CA', 
                        help='Atom selection string (default: protein and name CA)')
    
    # Output options
    parser.add_argument('-geig', '--geig_out', default='geigen.npy', 
                        help='Output file for generalized eigenpairs')
    parser.add_argument('-dkl', '--dkl_out', default='d_kl.png', 
                        help='Output file for KL divergence plot')
    parser.add_argument('-proj_a', '--proj_a_out', default='projection_a.npy',
                        help='Output file for state A projections')
    parser.add_argument('-proj_b', '--proj_b_out', default='projection_b.npy',
                        help='Output file for state B projections')
    parser.add_argument('-res', '--res_out', default=None,
                        help='Output file for residue interaction map')
    parser.add_argument('-respdb', '--respdb_out', default='contributions.pdb',
                        help='Output PDB with residue contributions as B-factors')
    
    # Analysis options
    parser.add_argument('-bt', '--begin_time', type=float, default=0,
                        help='Start time (ps) for analysis')
    parser.add_argument('-et', '--end_time', type=float, default=-1,
                        help='End time (ps) for analysis, -1 means until the end')
    parser.add_argument('-first', '--first_vec', type=int, default=0,
                        help='First eigenvector for analysis (0-based)')
    parser.add_argument('-last', '--last_vec', type=int, default=-1,
                        help='Last eigenvector for analysis, -1 means all')
    parser.add_argument('-algo', '--algorithm', type=int, default=0,
                        help='Algorithm for RPCA: 0=standard, 1=subspacing')
    parser.add_argument('-verbose', '--verbose', type=int, default=1,
                        help='Verbosity level (0-2)')
    
    # Performance options
    parser.add_argument('-n_jobs', '--n_jobs', type=int, default=None,
                        help='Number of parallel jobs (default: use all available cores)')
    parser.add_argument('-batch', '--batch_size', type=int, default=1000,
                        help='Batch size for trajectory processing')
    parser.add_argument('-fp32', '--use_float32', action='store_true',
                        help='Use single precision (float32) for faster computation but less accuracy')
    parser.add_argument('-mem', '--memory_efficient', action='store_true',
                        help='Use memory-efficient mode for large trajectories')
    
    args = parser.parse_args()
    
    # Set precision
    dtype = np.float32 if args.use_float32 else np.float64
    
    # Configure environment variables for optimal performance
    if args.n_jobs is not None:
        os.environ['OMP_NUM_THREADS'] = str(args.n_jobs)
        os.environ['MKL_NUM_THREADS'] = str(args.n_jobs)
    
    # Create RPCA object
    rpca = RPCAAnalysis()
    rpca.verbose = args.verbose
    
    # Print configuration
    print("=" * 50)
    print("RPCA Analysis Configuration:")
    print(f"  Trajectories: {args.traj_a}, {args.traj_b}")
    print(f"  Selection: {args.selection}")
    print(f"  Time range: {args.begin_time} to {args.end_time} ps")
    print(f"  Eigenvectors: {args.first_vec} to {args.last_vec}")
    print(f"  Algorithm: {args.algorithm}")
    print(f"  Precision: {dtype.__name__}")
    print(f"  Parallel jobs: {args.n_jobs}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Memory efficient: {args.memory_efficient}")
    print("=" * 50)
    
    # Measure total execution time
    start_time_total = time.time()
    
    # Read trajectories
    print(f"Reading state A trajectory: {args.traj_a}")
    universe_a = rpca.read_trajectory(args.traj_a, args.top_a, 
                                      args.begin_time, args.end_time)
    
    print(f"Reading state B trajectory: {args.traj_b}")
    universe_b = rpca.read_trajectory(args.traj_b, args.top_b,
                                      args.begin_time, args.end_time)
    
    # Perform GPA to find average structures
    print("\n" + "=" * 50)
    print("Performing GPA Analysis:")
    avg_coords_a = rpca.perform_gpa(universe_a, args.selection, n_jobs=args.n_jobs)
    avg_coords_b = rpca.perform_gpa(universe_b, args.selection, n_jobs=args.n_jobs)
    
    # Compute covariance matrices
    print("\n" + "=" * 50)
    print("Computing Covariance Matrices:")
    cov_a, mean_a = rpca.compute_covariance_matrix(universe_a, args.selection, avg_coords_a, n_jobs=args.n_jobs)
    cov_b, mean_b = rpca.compute_covariance_matrix(universe_b, args.selection, avg_coords_b, n_jobs=args.n_jobs)
    
    # Perform simultaneous diagonalization
    print("\n" + "=" * 50)
    print("Performing Simultaneous Diagonalization:")
    rpca_results = rpca.simultaneous_diagonalization(cov_a, cov_b, mean_a, mean_b)
    
    # Save results
    print("\n" + "=" * 50)
    print("Saving Results:")
    np.save(args.geig_out, rpca_results)
    print(f"Saved generalized eigenpairs to {args.geig_out}")
    
    # Plot KL divergence
    rpca.plot_kl_divergence(rpca_results['kl'], rpca_results['kl_m'], 
                          rpca_results['acc_kl'], args.dkl_out)
    
    # Set up eigenvector range for analysis
    if args.last_vec == -1:
        args.last_vec = rpca_results['rank'] - 1
    
    # Project trajectories
    print("\n" + "=" * 50)
    print("Projecting Trajectories:")
    proj_a = rpca.project_trajectory(universe_a, args.selection, 
                                   rpca_results['gevec'], mean_b,
                                   args.first_vec, args.last_vec,
                                   batch_size=args.batch_size,
                                   n_jobs=args.n_jobs)
    
    proj_b = rpca.project_trajectory(universe_b, args.selection,
                                   rpca_results['gevec'], mean_b,
                                   args.first_vec, args.last_vec,
                                   batch_size=args.batch_size,
                                   n_jobs=args.n_jobs)
    
    # Save projections
    np.save(args.proj_a_out, proj_a)
    np.save(args.proj_b_out, proj_b)
    print(f"Saved projections to {args.proj_a_out} and {args.proj_b_out}")
    
    # Plot projections
    rpca.plot_projections(proj_a, proj_b, args.first_vec, 'projections.png')
    
    # Compute interaction map
    if args.res_out:
        print("\n" + "=" * 50)
        print("Computing Interaction Map:")
        interaction_matrix, atom_contrib = rpca.compute_interaction_map(
            universe_b, args.selection, rpca_results['gevec'], rpca_results['kl'],
            args.first_vec, args.last_vec, n_jobs=args.n_jobs)
        
        np.save(args.res_out, interaction_matrix)
        print(f"Saved residue interaction map to {args.res_out}")
        
        # Save PDB with contributions as B-factors
        rpca.save_pdb_with_bfactors(universe_b, args.selection, 
                                  atom_contrib, args.respdb_out)
    
    # Report total execution time
    elapsed_time_total = time.time() - start_time_total
    print("\n" + "=" * 50)
    print(f"RPCA analysis completed in {elapsed_time_total:.2f} seconds!")
    
    # Summary of results
    print("\nSummary of RPCA Analysis:")
    print(f"  Total KL divergence: {rpca_results['sum_kl']:.4f}")
    print(f"  KL due to mean shift: {rpca_results['sum_kl_m']:.4f} ({100 * rpca_results['sum_kl_m'] / rpca_results['sum_kl']:.1f}%)")
    print(f"  KL due to covariance: {rpca_results['sum_kl'] - rpca_results['sum_kl_m']:.4f} ({100 * (rpca_results['sum_kl'] - rpca_results['sum_kl_m']) / rpca_results['sum_kl']:.1f}%)")
    print(f"  Rank of decomposition: {rpca_results['rank']}")
    print(f"  Number of eigenvectors analyzed: {args.last_vec - args.first_vec + 1}")
    print(f"  First 5 KL divergences: {', '.join([f'{kl:.4f}' for kl in rpca_results['kl'][:5]])}")
    
    # Suggest next steps
    print("\nSuggested next steps:")
    print("  - Examine the KL divergence plot to identify important components")
    print("  - Analyze projections to understand conformational differences")
    print("  - Inspect the PDB file with B-factors to identify key residues")
    print("  - Consider running with different atom selections for further analysis")
    print("=" * 50)


if __name__ == "__main__":
    main()
