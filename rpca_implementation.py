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
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import time


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
    
    def perform_gpa(self, universe, selection='protein', max_iterations=10, 
                    convergence=0.00001, ref_frame=0):
        """Perform Generalized Procrustes Analysis to find average structure
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            max_iterations: Maximum number of GPA iterations
            convergence: RMSD convergence criterion
            ref_frame: Initial reference frame number
            
        Returns:
            Average coordinates after GPA
        """
        print("Performing Generalized Procrustes Analysis (GPA)...")
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        
        # Initialize with the first frame as reference
        universe.trajectory[ref_frame]
        reference_coords = atoms.positions.copy()
        
        # Initialize average coordinates
        avg_coords = np.zeros_like(reference_coords)
        
        # Iterative GPA
        for iteration in range(max_iterations):
            print(f"GPA Iteration {iteration + 1}/{max_iterations}")
            
            # Reset average coordinates
            avg_coords.fill(0.0)
            total_frames = 0
            
            # Accumulate aligned coordinates
            for ts in universe.trajectory:
                # Align current frame to reference
                mobile_coords = atoms.positions
                
                # Center coordinates
                mobile_centered = mobile_coords - np.mean(mobile_coords, axis=0)
                ref_centered = reference_coords - np.mean(reference_coords, axis=0)
                
                # Calculate optimal rotation
                R, rmsd = self._calculate_rotation_matrix(mobile_centered, ref_centered)
                
                # Apply rotation
                rotated_coords = np.dot(mobile_centered, R)
                
                # Accumulate for average
                avg_coords += rotated_coords
                total_frames += 1
            
            # Calculate new average
            avg_coords /= total_frames
            
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
    
    def _calculate_rotation_matrix(self, mobile, reference):
        """Calculate optimal rotation matrix using SVD
        
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
        rotation = np.dot(Vt.T, np.diag([1, 1, d])).dot(U.T)
        
        # Calculate RMSD
        rotated = np.dot(mobile, rotation)
        rmsd = np.sqrt(np.mean(np.sum((rotated - reference)**2, axis=1)))
        
        return rotation, rmsd
    
    def compute_covariance_matrix(self, universe, selection='protein', average_coords=None):
        """Compute covariance matrix from trajectory
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            average_coords: Pre-computed average coordinates (optional)
            
        Returns:
            Covariance matrix and mean coordinates
        """
        print("Computing covariance matrix...")
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        n_dims = n_atoms * 3
        
        # Determine average coordinates if not provided
        if average_coords is None:
            print("Computing average coordinates...")
            average_coords = np.zeros((n_atoms, 3))
            n_frames = 0
            
            for ts in universe.trajectory:
                average_coords += atoms.positions
                n_frames += 1
                
            average_coords /= n_frames
        
        # Prepare for covariance calculation
        covariance_matrix = np.zeros((n_dims, n_dims))
        mean_coords = average_coords.flatten()
        n_frames = 0
        
        # Accumulate covariance matrix
        for ts in universe.trajectory:
            # Flatten coordinates
            coords = atoms.positions.flatten()
            
            # Compute deviation from mean
            deviation = coords - mean_coords
            
            # Accumulate covariance
            covariance_matrix += np.outer(deviation, deviation)
            n_frames += 1
            
        # Normalize by number of frames
        covariance_matrix /= n_frames
        
        return covariance_matrix, mean_coords
    
    def simultaneous_diagonalization(self, cov_a, cov_b, mean_a, mean_b):
        """Perform simultaneous diagonalization of two covariance matrices
        
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
        
        # Dimension of the problem
        n_dims = cov_a.shape[0]
        
        # Compute whitening transformation for matrix A
        evals_a, evecs_a = linalg.eigh(cov_a)
        
        # Filter out near-zero eigenvalues
        tol = np.sqrt(np.finfo(float).eps) * np.max(evals_a)
        mask = evals_a > tol
        rank = np.sum(mask)
        
        print(f"Matrix A has rank {rank} out of {n_dims} dimensions")
        
        # Compute whitening transformation
        W = evecs_a[:, mask] / np.sqrt(evals_a[mask])
        
        # Transform B with whitening transformation: W^T B W
        transformed_b = W.T @ cov_b @ W
        
        # Eigendecomposition of transformed B
        geigval, gevec_temp = linalg.eigh(transformed_b)
        
        # Full eigenvectors in original space
        gevec = W @ gevec_temp
        
        # Compute KL divergences
        # Mean difference vector
        mean_diff = mean_b - mean_a
        
        # KL divergence components
        kl = np.zeros(rank)
        kl_m = np.zeros(rank)
        
        for i in range(rank):
            # Compute mean contribution
            mean_proj = np.dot(gevec[:, i], mean_diff)
            kl_m[i] = 0.5 * mean_proj**2
            
            # Total KL divergence
            kl[i] = 0.5 * (-np.log(geigval[i]) + geigval[i] - 1) + kl_m[i]
        
        # Sort by KL divergence (decreasing order)
        sort_idx = np.argsort(-kl)
        kl = kl[sort_idx]
        kl_m = kl_m[sort_idx]
        geigval = geigval[sort_idx]
        gevec = gevec[:, sort_idx]
        
        # Accumulated KL divergence (percentage)
        kl_sum = np.sum(kl)
        acc_kl = np.cumsum(kl) / kl_sum * 100
        
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
    
    def project_trajectory(self, universe, selection, gevec, mean_coords, first_vec=0, last_vec=None):
        """Project trajectory onto generalized eigenvectors
        
        Args:
            universe: MDAnalysis Universe containing trajectory
            selection: Atom selection string for atoms to include in analysis
            gevec: Generalized eigenvectors 
            mean_coords: Mean coordinates
            first_vec: First eigenvector to include (0-based)
            last_vec: Last eigenvector to include (0-based), None means all
            
        Returns:
            Projections array, shape (n_frames, n_vecs)
        """
        print("Projecting trajectory onto eigenvectors...")
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        n_dims = n_atoms * 3
        
        # Set number of eigenvectors to use
        if last_vec is None:
            last_vec = gevec.shape[1] - 1
        
        n_vecs = last_vec - first_vec + 1
        
        # Allocate projections array
        n_frames = len(universe.trajectory)
        projections = np.zeros((n_frames, n_vecs))
        
        # Project each frame
        for i, ts in enumerate(universe.trajectory):
            # Flatten coordinates
            coords = atoms.positions.flatten()
            
            # Compute projection onto each eigenvector
            for j in range(n_vecs):
                vec_idx = first_vec + j
                projections[i, j] = np.dot(coords - mean_coords, gevec[:, vec_idx])
                
        return projections

    def compute_interaction_map(self, universe, selection, gevec, kl, first_vec=0, last_vec=None):
        """Compute per-residue contribution to conformational changes
        
        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string
            gevec: Generalized eigenvectors
            kl: KL divergences
            first_vec: First eigenvector to include (0-based)
            last_vec: Last eigenvector to include (0-based), None means all
            
        Returns:
            Residue interaction matrix and per-atom contributions
        """
        print("Computing interaction map...")
        
        # Select atoms for analysis
        atoms = universe.select_atoms(selection)
        n_atoms = len(atoms)
        
        # Set number of eigenvectors to use
        if last_vec is None:
            last_vec = gevec.shape[1] - 1
            
        n_vecs = last_vec - first_vec + 1
        
        # Get unique residues
        residues = atoms.residues
        n_res = len(residues)
        
        # Initialize interaction matrix
        interaction_matrix = np.zeros((n_res, n_res))
        
        # Compute per-atom contribution
        atom_contrib = np.zeros(n_atoms)
        
        # Iterate over included eigenvectors
        for k in range(n_vecs):
            vec_idx = first_vec + k
            vec = gevec[:, vec_idx].reshape(n_atoms, 3)
            weight = kl[vec_idx]
            
            # Update per-residue contributions
            for i, res_i in enumerate(residues):
                atoms_i = res_i.atoms.intersection(atoms)
                idx_i = np.array([atoms.indices[atoms.indices == a.index][0] for a in atoms_i])
                
                for j, res_j in enumerate(residues):
                    atoms_j = res_j.atoms.intersection(atoms)
                    idx_j = np.array([atoms.indices[atoms.indices == a.index][0] for a in atoms_j])
                    
                    # Compute contribution between residue i and j
                    contrib = 0
                    for a in idx_i:
                        for b in idx_j:
                            contrib += np.dot(vec[a], vec[b]) * weight
                    
                    interaction_matrix[i, j] += contrib
                    
                    # Update per-atom contribution for diagonal elements
                    if i == j:
                        for a in idx_i:
                            atom_contrib[a] += contrib / len(idx_i)
        
        # Normalize atom contributions
        if np.max(atom_contrib) > 0:
            atom_contrib /= np.max(atom_contrib)
            
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
    
    args = parser.parse_args()
    
    # Create RPCA object
    rpca = RPCAAnalysis()
    rpca.verbose = args.verbose
    
    # Read trajectories
    print(f"Reading state A trajectory: {args.traj_a}")
    universe_a = rpca.read_trajectory(args.traj_a, args.top_a, 
                                      args.begin_time, args.end_time)
    
    print(f"Reading state B trajectory: {args.traj_b}")
    universe_b = rpca.read_trajectory(args.traj_b, args.top_b,
                                      args.begin_time, args.end_time)
    
    # Perform GPA to find average structures
    avg_coords_a = rpca.perform_gpa(universe_a, args.selection)
    avg_coords_b = rpca.perform_gpa(universe_b, args.selection)
    
    # Compute covariance matrices
    cov_a, mean_a = rpca.compute_covariance_matrix(universe_a, args.selection, avg_coords_a)
    cov_b, mean_b = rpca.compute_covariance_matrix(universe_b, args.selection, avg_coords_b)
    
    # Perform simultaneous diagonalization
    rpca_results = rpca.simultaneous_diagonalization(cov_a, cov_b, mean_a, mean_b)
    
    # Save results
    np.save(args.geig_out, rpca_results)
    print(f"Saved generalized eigenpairs to {args.geig_out}")
    
    # Plot KL divergence
    rpca.plot_kl_divergence(rpca_results['kl'], rpca_results['kl_m'], 
                          rpca_results['acc_kl'], args.dkl_out)
    
    # Set up eigenvector range for analysis
    if args.last_vec == -1:
        args.last_vec = rpca_results['rank'] - 1
    
    # Project trajectories
    proj_a = rpca.project_trajectory(universe_a, args.selection, 
                                   rpca_results['gevec'], mean_b,
                                   args.first_vec, args.last_vec)
    
    proj_b = rpca.project_trajectory(universe_b, args.selection,
                                   rpca_results['gevec'], mean_b,
                                   args.first_vec, args.last_vec)
    
    # Save projections
    np.save(args.proj_a_out, proj_a)
    np.save(args.proj_b_out, proj_b)
    print(f"Saved projections to {args.proj_a_out} and {args.proj_b_out}")
    
    # Plot projections
    rpca.plot_projections(proj_a, proj_b, args.first_vec, 'projections.png')
    
    # Compute interaction map
    if args.res_out:
        interaction_matrix, atom_contrib = rpca.compute_interaction_map(
            universe_b, args.selection, rpca_results['gevec'], rpca_results['kl'],
            args.first_vec, args.last_vec)
        
        np.save(args.res_out, interaction_matrix)
        print(f"Saved residue interaction map to {args.res_out}")
        
        # Save PDB with contributions as B-factors
        rpca.save_pdb_with_bfactors(universe_b, args.selection, 
                                  atom_contrib, args.respdb_out)
    
    print("RPCA analysis complete!")


if __name__ == "__main__":
    main()
