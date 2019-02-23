#ifndef _covariance_h
#define _covariance_h

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "linear_algebra.h"
#include "rmpbc.h"
#include "typedefs.h"



/////////////////////////////////////////////////////
/**************************************************************
structure for the information of generalized eigenanalysis
**************************************************************/
			 
typedef struct 
{
	int ndim;			// nr. of dims
	int 	rank;		// rank of mat_a = nr. of elements of kl, geival ,and nr. of vectors in gevec
	double  *mat_a;		// matrix A    gevec^t*A*gevec=I  
	double  *mat_b;		// matrix B    gevec^t*B*gevec=diag(geigval)
	double  *m_a;		// average of state A
	double  *m_b;		// average of state B
	double  *geigval;	// generalized eigenvalues
	double  *gevec;		// generalized eigenvectors in col. wise 1D array
	double  *kl;		// KL divergences of the eigenvectors
	double  *acc_kl;	// accumulated KL divergences of the eigenvectors
	double  *kl_m;		// KL divergences of the eigenvectors due to the mean change
	int 	ordered;    // are geigval, gevec and kl reordered according to kl values
	double  sum_kl;		// sum of kl
	double  sum_kl_m;	// sum of kl_m
	int     alg;		// 1 mean-fluctuation subspacing, 0 otherwise
} Geigdecomposition_t ;			 
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void blas_compute_covariance
			(const char 	*trxfile,   // traj file to read
			rvec  		*xref,      // refrence coor structure to be used for fitting 
			double 		*mat,  // pointer where the cov mat is allocated using 1D array
			int 		nfit,   // nr atoms used for fitting
			atom_id 	*ifit,  // index of ftting atoms
			int 		natom,   // nr atoms used for used for covariance analysis
			atom_id 	*index,  // index of atoms used for covariance analysis
			gmx_bool 	bPBC,
			gmx_rmpbc_t     gpbc,
			double     *m,    /// array of the average atoms
			float first_time,
			float last_time
		     )	;	 

//////////////////////////////////////////////////////////////////	
void  cov2datfile(double *mat, int ndim, const char *covmatfile);
////////////////////////////////////////////////////////////////////////
//////////////// wrtting and reading covariance matrix to a binary file
///int write_read_sym_matrix(const char *outfile,double **mat , int *ndim, int W);
int write_cov_matrix(const char *outfile,double **mat , double **average,int *ndim,  int verbose);
int read_cov_matrix(const char *outfile,double **mat , double **average,int *ndim);		
////////////////////////////////////////////////////////////////////////	  
int write_read_eigen(const char *outfile,eigdecomposition_t *eig , double **average, int W);
int write_eigen(const char *outfile,eigdecomposition_t *eig , double **average, int verbose);
int read_eigen(const char *outfile,eigdecomposition_t *eig , double **average);
////////////////////////////////////////////////////////////////////////
int write_read_Geigen(const char *outfile,Geigdecomposition_t *eig ,int W);
int read_Geigen(const char *outfile,Geigdecomposition_t *eig );
int write_Geigen(const char *outfile,Geigdecomposition_t *eig , int verbose);
/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//void write_eigenvalues(const char *eigvalfile, int ndim, double *eigval)		;
void write_array2xvg(const char *eigvalfile, int ndim, double *eigval, char *title ,char *xaxis, char *yaxis  )	;
////////////////////////////////////////////////////////////////////////


 void sort_decreasing( double *x , int n, int *x_indx);
 double max_array(double *a, int n);

////////////////////////////////////////////


void relative_pca(Geigdecomposition_t *G,int ordered ,int algo,int verbose );
//////////////////////////////////////////////////////////////////////////

void relative_pca_subspacing(Geigdecomposition_t *G,int ordered ,int verbose );
//////////////////////////////////////////////////////////////////////////

void write_Geigen2xvg(const char *Geigfile, Geigdecomposition_t *G );	
//////////////////////////////////////////////////////////////////////////

void interaction_map(const char *resfile, Geigdecomposition_t *eig, t_atoms *atoms,atom_id *ifit,int first ,int last,  const char *pdbBfactorsfile);
//////////////////////////////////////////////////////////////////////////

void RPCA_project(const char 	*trxfile,   // traj file to read
					int          state_A,   /// project state A or B
					const char *projfile ,    // write the projections to a dat file
					Geigdecomposition_t *eig, 
					int Beg_evec,
					int End_evec,
					int 		natom,
					atom_id 	*index,
					gmx_bool 	bPBC,
					gmx_rmpbc_t     gpbc,
					float  timeBeg,
					float timeEnd,
					int skip           /// nr of frames to be skipped between projected frames
					);


////////////////////////////////////////////////////////////////////////
void PCA_project(const char 	*trxfile,   // traj file to read
					const char *projfile ,    // write the projections to a dat file
					eigdecomposition_t *eig,
					double             *m,   //  the average
					int Beg_evec,
					int End_evec,
					int 		natom,
					atom_id 	*index,
					gmx_bool 	bPBC,
					gmx_rmpbc_t     gpbc,
					float first_time,
					float last_time
					);
///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void interpolate_coor(char *pdbfname,int r, int c , double *G ,double *m_x, double *ev ,char *title, int natom , t_atoms *atoms,atom_id *index, double sd );
///////////////////////////////////////////////////////////////////////////

void show_motions_along_vectors(	const char *pdbfile ,    // write the structures to the file
							Geigdecomposition_t *geig,
							int 		natom,
							atom_id 	*index,
							int first_vec,  /// index of first vector
							int last_vec,   /// index of the last vector
							t_atoms *atoms,
							int write_individual);  /// write pdb for the changes along each vector besides the total

/////////////////////////////////////////////////////////////////////////////////
void show_increased_decresed_motions(	const char *pdbfile ,    // write the projections to a dat file
							Geigdecomposition_t *geig, 
							int 		natom,
							atom_id 	*index,
							double tol,
							int nr_vec_up,  /// number of vectors (high eigenvalues) to be used if tol =0
							int nr_vec_low,  /// number of vectors (high eigenvalues) to be used if tol =0
							t_atoms *atoms,
							int write_individual) ; /// write pdb for the changes along each vector besides the total
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
#endif	/* _covariance_h*/
