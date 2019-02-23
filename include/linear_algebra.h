/*
 * linear_algebra.h
 *
 *  Created on: Mar 29, 2018
 *      Author: mahmad
 */

#ifndef INCLUDE_LINEAR_ALGEBRA_H_
#define INCLUDE_LINEAR_ALGEBRA_H_
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
/////////////////////////////////////////////////////////////////////////

#define BLAS_CALL(x) cblas_##x

/////////////////////////////////////////////////////////////////////////
typedef struct
{
	int     ndim;		// nr. of dims
	int 	rank;		// rank of mat
	double  *mat;		// matrix
	double  *eval;		//  eigenvalues
	double  *vec;		//  eigenvectors in col. wise 1D array
	int     precision;	// 1 double ,0 single
} eigdecomposition_t ;
///////////////////////////////////////////

//////////////////////////// cast array from float to double
void cast_F2D(int n, float *Xf, double *Xd) ;

///////////////////////////cast array from  double float
void cast_D2F(int n, float *Xf, double *Xd);
////////////////////////////////////////////
int IsIdentity(double *X, int n, double tol, int verbose);
//// check if the symm mat X is diagonal
int IsDiagonal(double *X, int n, double tol, int verbose);
////////////////////////////////////////////
void mat2fvec(int m, int n, double ** mat, double * fvec) ;
void fvec2mat(int m, int n, double * fvec, double ** mat) ;
void One_DIM_matrix_transpose(int m, int n, double *X, double *Y);
void write_Fvec_matrix2file( double *fmat, int nr, int nc, const char *outfile );
void lapack_ggemm_mat_prod( double* A, int m, int k , double* B,int n, double* C);

////////////////////////////////////////////////////////////////////////

void eigensolver_1D(	float *   af,     		// 1D matrix  as a col.wise vector of size n*n
					int      n,        		// Side of the matrix to calculate eigenvalues for.
					int      index_lower, 	// Index of first eigenvector to determine.
					int      index_upper,   // Last eigenvector determined is index_upper-1.
					float *   eigenvalues,   // 1 D array of the eigenvalues
					float *   eigenvectors     // 1 D array of the eigenvectors filled col.wise
				) ;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
eigensolver(	float **   a2d,     		// matrix as 2D array
				int      n,        		// Side of the matrix to calculate eigenvalues for.
				int      index_lower, 	// Index of first eigenvector to determine.
				int      index_upper,   // Last eigenvector determined is index_upper-1.
				float *   eigenvalues,   // 1 D array of the eigenvalues
				float **   eigenvectors     // 2 D array of the eigenvectors in the rows such that eigenvector i is allocated in eigenvec[i-1][]

			);


////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// double  precision eignen analysis //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void double_eigensolver_1D(	double *   af,     		// 1D matrix  as a col.wise vector of size n*n
					int      n,        		// Side of the matrix to calculate eigenvalues for.
					int      index_lower, 	// Index of first eigenvector to determine.
					int      index_upper,   // Last eigenvector determined is index_upper-1.
					double *   eigenvalues,   // 1 D array of the eigenvalues
					double *   eigenvectors     // 1 D array of the eigenvectors filled col.wise
				);

///////////////////////////////////////////////////////////////////////////////////////
///////////////// warpper for double pre. cacl of eigenpairs of 2D matrix ////////////
void double_eigensolver(	double **   a2d,     		// matrix as 2D array
					int      n,        		// Side of the matrix to calculate eigenvalues for.
					int      index_lower, 	// Index of first eigenvector to determine.
					int      index_upper,   // Last eigenvector determined is index_upper-1.
					double *   eigenvalues,   // 1 D array of the eigenvalues
					double **   eigenvectors     // 2 D array of the eigenvectors in the rows such that eigenvector i is allocated in eigenvec[i-1][]
				);
/////////////////////////////////////////////////////////////////////////////////////////

void eigen(double *mat,int n, eigdecomposition_t *eig );

void eigen2matrix(eigdecomposition_t *eig, double **mat) ;

//////////////////////////////////////////////////////////////////////////////////////////
//////// find the generalized inverse of a symmetric matrix X of dim n*n           ///////
///////   the inverse pinv_X is of dim n*n
void symmpinv(int n, double *X, double *pinv_X) ;
///////////// generalized inverse for float (internally done using double). output is float
void symmpinvF2F(int n, float *X, float *pinv_X);
///////////// generalized inverse for float (internally done using double). output is double
void symmpinvF2D(int n, float *X, double *pinv_X);

//////////////////////
void Mahalanobis_transformation(int n, double *X, double *mah_X, int *rank);
void withening_transformation(int n, double *X, double **W, int *rank, int verbose) ;
///////////////////////////////////////////////////////////////////////////////////////////////
///// compute x^t*MAT*y   where MAT is symm.
double   vecX_symat_vecY( double *mat, int n  , double *x , double *y  );
void quadratic_Mat_multiply(double *B, int n  , double *W , int r , double *WtBW ) ;
double mah_norm(double *ginv, int n  , double *x  ) ;
float Fmah_norm(float *ginv, int n  , float *x  ) ;


void Simultaneous_Diagonalization(int n, double *A, double *B, double **geigval, double **gevec, int *rank , int verbose ) ;

////////////////////////////////////////////////////////////////////////////////////////////
///// complement projection on vector v
///// if ref is not null: projection w.r.t. ref
void Proj_on_vec(double *v, int n , double *mat ,double *proj, double *ref , int complement);


////////////////////////////////////////////////////////////////////////

void Simultaneous_Diagonalization_subspacing(int n, double *A, double *B, double **geigval, double **gevec, int *rank ,double *dav  ,int verbose);
void Simultaneous_Diagonalization_last3(int n, double *A, double *B, double **geigval, double **gevec, int *rank ,double *dav ,int verbose) ;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////  projection matrix on the subspace spanned by the vectors of a matrix A (r x c)
///// P = A(A^tA)^-1 A^t
void Projection_matrix(int r, int c , double *A ,double *P);

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#endif /* INCLUDE_LINEAR_ALGEBRA_H_ */
