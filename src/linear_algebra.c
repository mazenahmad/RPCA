

#include  <math.h>
#include <float.h>
#include <time.h>
#include <cblas.h>
#include "linear_algebra.h"
#include "arrays.h"



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void cast_F2D(int n, float *Xf, double *Xd)
{
	int i;
	for (i=0; i< n  ; i++) Xd[i]= (double) Xf[i];
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void cast_D2F(int n, float *Xf, double *Xd)
{
	int i;
	for (i=0; i< n  ; i++) Xf[i]= (float) Xd[i];
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//// check if the symm mat X is identity
int IsIdentity(double *X, int n, double tol ,int verbose)
{
	int i,j, out=1;
	for (i=0; i< n; i++) for (j=0; j< n; j++) 
	{
		
		if( i!=j && (X[i*n +j] > tol ) ) 
		{	
			if(verbose) printf("Non diagonal (%d  %d) is not zero= %f\n",i, j, X[i*n +j]);
			out =0;
		}
		if(i==j &&  (fabs(X[i*n +j] -1) > tol )) 
		{
			if(verbose) printf("diagonal (%d  %d) = %f\n",i, j, X[i*n +j]);
			out =0;
		}
	}
	return out ;
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//// check if the symm mat X is diagonal
int IsDiagonal(double *X, int n, double tol, int verbose)
{
	int i,j, out=1;
	for (i=0; i< n; i++) for (j=0; j< n; j++) 
	{
		
		if( i!=j && (X[i*n +j] > tol ) ) 
		{	
			if(verbose) printf("Non diagonal (%d  %d) is not zero= %f\n",i, j, X[i*n +j]);
			out =0;
		}
	}
	return out ;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
 void fvec2mat(int m, int n, double * fvec, double ** mat)
 {
        /*
                Turns a matrix in (fortran) column vector format into a matrix given by a double  pointer
                The matrices have dimension m*n
        */
        int i,j;
        for (i=0; i<m; i++) for (j=0; j<n; j++)  mat[i][j] = fvec[i+m*j];
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
 void mat2fvec(int m, int n, double ** mat, double * fvec) 
 {
        /*
                Turns a matrix (given by a double  pointer) into its Fortran vector format 
                (single vector, colunmwise). The matrix "mat" needs to be an m*n matrix
                The vector "vec" is a vector of lenght mn
        */
        int i,j;
        for (i=0; i<m; i++) for (j=0; j<n; j++) fvec[i+m*j] = mat[i][j];
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void One_DIM_matrix_transpose(int m, int n, double *X, double *Y) 
{
        ////
        ////   Transposes an m*n matrix: Y = X^T
        ////   Matrix can be in either C or Fortran vector format
        int i,j;
        for (i=0; i<m; i++) for (j=0; j<n; j++) Y[n*i+j] = X[i+m*j];
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void write_Fvec_matrix2file( double *fmat, int nr, int nc, const char *outfile )
{
	
	int i,j;
	FILE           *out;
	out = fopen(outfile, "w");
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++)	fprintf(out, "%g \t", fmat[i+nr*j]);
		fprintf(out, " \n");
	}
	fclose(out);

}
////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
void lapack_ggemm_mat_prod( double* A, int nrows_A, int ncol_A , double* B,int ncol_B, double* C)
{
    ////// dimensions: A is nrows_A*ncol_A, B is ncol_A*ncol_B, so C has to be nrows_A*ncol_B
    ///// Calculates C = a*A*B + b*C: a =1,  b=0
    ///// A, B and C need to be matrices in Fortran vector format (single vector, columnwise)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows_A, ncol_B, ncol_A, 1.0, A, nrows_A, B, ncol_A, 0.0, C, nrows_A);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
///////////////////         single precision eignen analysis      ////////////////
///////////////////////////////////////////////////////////////////////////////////

void eigensolver_1D(	float *   X,     	// 1D matrix  as a col.wise vector of size n*n
					int      n,        		// Side of the matrix to calculate eigenvalues for.
					int      index_lower, 	// Index of first eigenvector to determine.
					int      index_upper,   // Last eigenvector determined is index_upper-1.
					float *   eigenvalues,   // 1 D array of the eigenvalues
					float *   eigenvectors     // 1 D array of the eigenvectors filled col.wise
				)
{

	extern  void  ssyevr_ (const char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                 float *A, int *LDAp, float *VLp, float *VUp,
                 int *ILp, int *IUp, float *ABSTOLp, int *Mp,
                 float *W, float *Z, int *LDZp, int *ISUPPZ,
                 float *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                 int *INFOp);
	ONE_DIM_ARRAY(af, float, (n*n));
	cblas_scopy ((n*n),  X, 1, af, 1);
	int i, j;
	int           lwork, liwork;
    int           il, iu, m, iw0, info;
    int       *   iwork;
    const char *  jobz;
    if (index_lower < 0)
    {
        index_lower = 0;
    }
    if (index_upper >= n)
    {
        index_upper = n-1;
    }

    // Make jobz point to the character "V" if eigenvectors should be calculated, otherwise "N" (only eigenvalues).
	jobz = "V";
  	ONE_DIM_ARRAY(isuppz, int, 2*n);
    ////  ask the routine how much workspace it needs 
    lwork  = -1;
    liwork = -1;
    // Convert indices to fortran standard 
    index_lower++;
    index_upper++;
	float          w0, abstol;
	float       *  work;
    float          vl, vu;
	vl = vu = 0;
	ONE_DIM_ARRAY(eigenvalues_tmp, float, n );
	ONE_DIM_ARRAY(eigenvectors_tmp, float, (n*n) );
	ssyevr_(jobz,   "I",    "L",  &n, af,  &n,  &vl, &vu, &index_lower, &index_upper, &abstol, &m, eigenvalues_tmp, eigenvectors_tmp, &n, isuppz, &w0,  &lwork, &iw0,  &liwork, &info);
    lwork  = w0;
	ALLO_ONE_DIM_ARRAY(work,float ,lwork);
	if (info != 0)
    {
		free(isuppz);
        PRINT_FATAL_ERR_MES("Internal errror in LAPACK diagonalization.");
    }
    liwork = iw0;
	ALLO_ONE_DIM_ARRAY(iwork, int, liwork);
    abstol = 0;
	
	ssyevr_(jobz, "I", "L", &n, af, &n, &vl, &vu, &index_lower, &index_upper, &abstol, &m, eigenvalues_tmp, eigenvectors_tmp, &n, isuppz, work, &lwork, iwork, &liwork, &info);
    ////////////// reorder the eigenvectors and the eigen vaules
	free(af);
	for (i=0; i< n; i++) for (j=0; j< n; j++)  eigenvectors[(n-i-1)*n +j] = eigenvectors_tmp[i*n +j];
	free(eigenvectors_tmp);	
	for (i=0; i< n; i++) eigenvalues[i]= eigenvalues_tmp[n-i-1];
	free(eigenvalues_tmp);
	free(work);
  	free(isuppz);
    free(iwork);
    if (info != 0) PRINT_FATAL_ERR_MES("Internal errror in LAPACK diagonalization.");
		
 }  
  
////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// double  precision eigen  analysis ///////////////////
////////////////////////////////////////////////////////////////////////////////////
 
void double_eigensolver_1D(
					double *   X,     		// 1D matrix  as a col.wise vector of size n*n 
					int      n,        		// Side of the matrix to calculate eigenvalues for.
					int      index_lower, 	// Index of first eigenvector to determine. (c type <=> first =0)
					int      index_upper,   // Last eigenvector determined is index_upper-1.
					double *   eigenvalues,   // 1 D array of the eigenvalues
					double *   eigenvectors     // 1 D array of the eigenvectors filled col.wise   
				)
{
	extern  void  dsyevr_ (const char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                 double *A, int *LDAp, double *VLp, double *VUp,
                 int *ILp, int *IUp, double *ABSTOLp, int *Mp,
                 double *W, double *Z, int *LDZp, int *ISUPPZ,
                 double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                 int *INFOp);
	ONE_DIM_ARRAY(af, double, (n*n));
	cblas_dcopy ((n*n),  X, 1, af, 1);
	int i, j;
	int           lwork, liwork;
    int           il, iu, m, iw0, info;
    int       *   iwork;
    const char *  jobz;
    if (index_lower < 0)
    {
        index_lower = 0;
    }
    if (index_upper >= n)
    {
        index_upper = n-1;
    }

    // Make jobz point to the character "V" if eigenvectors should be calculated, otherwise "N" (only eigenvalues).
	jobz = "V";
  	ONE_DIM_ARRAY(isuppz, int, 2*n);
    ////  ask the routine how much workspace it needs 
    lwork  = -1;
    liwork = -1;
    // Convert indices to fortran standard 
    index_lower++;
    index_upper++;
	double          w0, abstol;
	double       *  work;
    double          vl, vu;
	vl = vu = 0;
	ONE_DIM_ARRAY(eigenvalues_tmp, double, n );
	ONE_DIM_ARRAY(eigenvectors_tmp, double, (n*n) );
	dsyevr_(jobz,   "I",    "L",  &n, af,  &n,  &vl, &vu, &index_lower, &index_upper, &abstol, &m, eigenvalues_tmp, eigenvectors_tmp, &n, isuppz, &w0,  &lwork, &iw0,  &liwork, &info);
    lwork  = w0;
	ALLO_ONE_DIM_ARRAY(work,double ,lwork);
	if (info != 0)
    {
		free(isuppz);
        PRINT_FATAL_ERR_MES("Error in LAPACK dsyevr_ function.");
    }
    liwork = iw0;
	ALLO_ONE_DIM_ARRAY(iwork, int, liwork);
    abstol = 0;
	
	dsyevr_(jobz, "I", "L", &n, af, &n, &vl, &vu, &index_lower, &index_upper, &abstol, &m, eigenvalues_tmp, eigenvectors_tmp, &n, isuppz, work, &lwork, iwork, &liwork, &info);
    ////////////// reorder the eigenvectors and the eigen vaules
	free(af);
	for (i=0; i< n; i++) for (j=0; j< n; j++)  eigenvectors[(n-i-1)*n +j] = eigenvectors_tmp[i*n +j];
	free(eigenvectors_tmp);	
	for (i=0; i< n; i++) eigenvalues[i]= eigenvalues_tmp[n-i-1];
	free(eigenvalues_tmp);
	free(work);
  	free(isuppz);
    free(iwork);
    if (info != 0) PRINT_FATAL_ERR_MES("Error in LAPACK dsyevr_ function.");
 }  
 
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void eigen(double *mat,int n, eigdecomposition_t *eig )
{
	eig->mat = mat;
	int i ;
	eig->ndim =n;
	ALLO_ONE_DIM_ARRAY( (eig->eval), double,n);
	ALLO_ONE_DIM_ARRAY(eig->vec, double, (n* n));
	double tol;
	double_eigensolver_1D(	eig->mat, n,  0, n, eig->eval,   eig->vec);
	eig->rank= n ;
	tol =  sqrt(DBL_EPSILON) * eig->eval[0] ;
	if (tol < 0.0) tol= 0.0;
	for (i=(n -1) ; i >= 0 ; i--) 
	{
		if (eig->eval[i]  < tol) eig->rank -=1;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// convert an eigenpair back to the symm matrix  //////////////////////

void eigen2matrix(eigdecomposition_t *eig, double **mat) 
{
    int n = eig->ndim;
    int i, j;
	double tol, ev, *px;
	ALLO_ONE_DIM_ARRAY(*mat, double, (n*n));
	//// this is how tolerence is defined in MASS library in R
	tol =  sqrt(DBL_EPSILON) * eig->eval[0] ;
	if (tol < 0.0) tol= 0.0;
	for (i=0; i< (n*n); i++) 	(*mat)[i]=0;
	for (i=0; i< n; i++)
	{
		if (eig->eval[i] > tol)
		{
			ev= eig->eval[i] ;
			px = (eig->vec) + i*n ;			
			// DSYR   performs the symmetric rank 1 operation    A := alpha*x*x**T + A,
			//  where alpha is a double scalar, x is an n element vector and A is an  n by n symmetric matrix.
			cblas_dsyr( CblasColMajor ,CblasLower, n ,ev, px ,1,(*mat) ,n );
		}
	}
	//// fill the upper part of the sym mat
	for (i=0; i< n; i++) for (j=i; j< n; j++) (*mat)[j*n +i]=(*mat)[i*n +j] ;
/*****************************
	/////////////////////////// for testing //////////////////////////////////
	printf("testing\n");
	eigdecomposition_t eig2;
    eigen((*mat), n,&eig2);
	for (i = 0; i < n; i++) if ((eig->eval[i] - eig2.eval[i]) > tol) printf("%e  %e \n", eig->eval[i] , eig2.eval[i]);
	int v;
	for (v = 0; v < n; v++) for (i = 0; i < n; i++) if((eig->vec[v*n +i] - eig2.vec[v*n+i]) > tol) 
														printf("V= %d %e  %e \n",v, eig->vec[v*n +i] , eig2.vec[v*n+i]);
		
***********************************/	
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//////// find the generalized inverse of a symmetric matrix X of dim n*n     ///////
///////   the inverse pinv_X is of dim n*n
void symmpinv(int n, double *X, double *pinv_X) 
{
    int i, j;
	double tol, ev, *px;
	ONE_DIM_ARRAY(eigval, double,n );
	ONE_DIM_ARRAY(eigvec, double, (n*n));
	double_eigensolver_1D(	X, n,  0, 	n,   eigval,   eigvec);
	tol =  sqrt(DBL_EPSILON) * eigval[0] ;
	if (tol < 0.0) tol= 0.0;
	for (i=0; i< (n*n); i++) 	pinv_X[i]=0;
	for (i=0; i< n; i++)
	{
		if (eigval[i] > tol)
		{
			ev= 1/eigval[i] ;
			px = eigvec + i*n ;
			// DSYR   performs the symmetric rank 1 operation    A := alpha*x*x**T + A,
			//  where alpha is a double scalar, x is an n element vector and A is an  n by n symmetric matrix.
			cblas_dsyr( CblasColMajor ,CblasLower, n ,ev, px ,1,pinv_X ,n );
		}
	}
	//// fill the upper part of the sym mat
	for (i=0; i< n; i++) for (j=i; j< n; j++) pinv_X[j*n +i]=pinv_X[i*n +j] ;
	free(eigval);
	free(eigvec);
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///////////// generalized inverse for float (internally done using double). output is float
void symmpinvF2F(int n, float *X, float *pinv_X)
{
	int i;
	ONE_DIM_ARRAY(at, double, (n*n));
	for (i=0; i<(n*n) ; i++) at[i]= (double) X[i];
	/// comput generalized inverse
	ONE_DIM_ARRAY(ginv, double, (n*n));
	symmpinv(n, at, ginv);
	free(at);
	for (i=0; i<(n*n) ; i++) pinv_X[i]= (float) ginv[i];
	free(ginv);
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///////////// generalized inverse for float (internally done using double). output is double
void symmpinvF2D(int n, float *X, double *pinv_X)
{
	int i;
	ONE_DIM_ARRAY(at, double, (n*n));
	for (i=0; i<(n*n) ; i++) at[i]= (double) X[i];
	/// comput generalized inverse
	symmpinv(n, at, pinv_X);
	free(at);
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// find the Mahalanobis transformation = the inverse of sqr root of a matrix of dim n*n ///////
// /// equations(A.6.15) & (A.6.10) of Maridia 1979									    ///////
////////  rank: in the nr. of  nonzero eigenvalues                                       //////
void Mahalanobis_transformation(int n, double *X, double *mah_X, int *rank) 
{
    //clock_t begin = clock();
    int i, j;
	double tol, ev, *px;
	ONE_DIM_ARRAY(eigval, double,n );
	ONE_DIM_ARRAY(eigvec, double, (n*n));
	double_eigensolver_1D(	X, n,  0, 	n,   eigval,   eigvec);
	//free(Y);
	*rank = n;
	//// this is how tolerence is defined in MASS library in R
	tol =  sqrt(DBL_EPSILON) * eigval[0] ;
	if (tol < 0.0) tol= 0.0;
	//printf("eps= %g , tol= %g \n",DBL_EPSILON, tol );
	for (i=0; i< (n*n); i++) 	mah_X[i]=0;
	for (i=0; i< n; i++)
	{
		if (eigval[i] > tol)
		{
			ev= 1/sqrt(eigval[i]) ;
			px = eigvec + i*n ;
			// DSYR   performs the symmetric rank 1 operation    A := alpha*x*x**T + A,
			//  where alpha is a double scalar, x is an n element vector and A is an  n by n symmetric matrix.
			cblas_dsyr( CblasColMajor ,CblasLower, n ,ev, px ,1,mah_X ,n );
		}
		else *rank -=1 ;
	}
	//// fill the upper part of the sym mat
	for (i=0; i< n; i++) for (j=i; j< n; j++) mah_X[j*n +i]=mah_X[i*n +j] ;
	free(eigval);
	free(eigvec);
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// find the Whitening Transformation from the covariance matrix X of dim n*n ///////
//////// equations (2.88) & (2.89) of Fukunaga(1990)   			   ////////////
////////  rank: in the nr. of  nonzero eigenvalues                     //////
//void withening_transformation(int n, double *X, double *W_X, int *rank) 
void withening_transformation(int n, double *X, double **W, int *rank , int verbose) 
{
	int i, j, ndx, wndx;
	double tol, ev, *px;
	ONE_DIM_ARRAY(eigval, double,n );
	ONE_DIM_ARRAY(evec, double, (n*n));
	double_eigensolver_1D(	X, n,  0, 	n,   eigval,   evec);
	*rank = n;
	//// this is how tolerence is defined in MASS library in R
	tol =  sqrt(DBL_EPSILON) * eigval[0] ;
	if (tol < 0.0) tol= 0.0;
	if (verbose > 1) printf("tolerence for zero value is : %g \n ", tol);
	for (i=0; i< n; i++) if (eigval[i] <= tol) *rank -=1 ;
	ALLO_ONE_DIM_ARRAY((*W), double, (n*(*rank)));
	//printf("eps= %g , tol= %g \n",DBL_EPSILON, tol );
	for (i=0; i< (*rank); i++)
	{
		ndx= i*n ;
		wndx= ((*rank)-i-1)*n;
		ev= 1/sqrt(eigval[i]) ;
		for (j=0; j< n; j++) (*W)[wndx +j] =  ev* evec[ndx +j] ;
	}
	free(eigval);
	free(evec);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///// compute x^t*MAT*y   where MAT is symm (n,n). x and y are n elements vectors 
////////////////////////////////////////////////////////////////////////////////////
double   vecX_symat_vecY( double *mat, int n  , double *x , double *y  )
{
	ONE_DIM_ARRAY(tt, double, n);
	cblas_dsymv ( CblasColMajor,  CblasLower, n, 1, mat, n, y, 1, 0, tt, 1);
	double d = cblas_ddot (n, x, 1, tt, 1);
	free(tt);
	return d;
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void quadratic_Mat_multiply(double *B, int n  , double *W , int r , double *WtBW ) 
{
	/***************************************************************************
	compute the the multiply of 3 mats    WtBW = t(w)%*%B %*%w   
	B is symmetric n*n 
	W is n*r
	WtBW is r*r
	*****************************************************************************/
	///////// first compute  BW= B%*% W
	ONE_DIM_ARRAY(BW, double, (n* r));
	///  DSYMM  performs one of the matrix-matrix operations
	///        C := alpha*A*B + beta*C,
	/// or     C := alpha*B*A + beta*C,
	/// cblas_dsymm (const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const int M, const int N, const double alpha,
	///          const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
	///  where alpha and beta are scalars,  A is a symmetric matrix and  B and C are  m by n matrices.
	/// M  specifies the number of rows of the matrix  C. N specifies the number of columns of the matrix C
	cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, n, r, 1, B, n, W, n, 0, BW, n);
	//////// compute W^tBW  
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, r, r, n, 1.0, W, n, BW, n, 0.0, WtBW, r);
	free(BW);
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
double mah_norm(double *ginv, int n  , double *x  ) 
{
	ONE_DIM_ARRAY(y, double, n);
	cblas_dsymv ( CblasColMajor,  CblasLower, n, 1, ginv, n, x, 1, 0, y, 1);
	double d = cblas_ddot (n, x, 1, y, 1);
	free(y);
	return d;
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
float Fmah_norm(float *ginv, int n  , float *x  ) 
{
	ONE_DIM_ARRAY(y, float, n);
	cblas_ssymv ( CblasColMajor,  CblasLower, n, 1, ginv, n, x, 1, 0, y, 1);
	float d = cblas_sdot (n, x, 1, y, 1);
	free(y);
	return d;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
///////////   Simultaneous Diagonalization: Introduction to Statistical  ///////////
//////////    Pattern Recognition Fukunaga(1990) p.31                   ////////////
////////      B with w.r. to A   
////////      A and B have to be in col.wise form
/*********************************************************************************
solve generalized eigen value problem in general case 
Generalized eigenvectors are allocated (inside the function) to gevec (n*r).
Generalized eigenvalues are allocated (inside the function) to geigval (n*r).
Note: privide the pointer to the pointer for both of gevec and geigval to enable 
      the allocation inside the function
G= gevec :
G^t %*% A %*% G =I
G^t %*% B %*% G = diag(geigval)
**********************************************************************************/

void Simultaneous_Diagonalization(int n, double *A, double *B, double **geigval, double **gevec, int *rank , int verbose) 
{
	//////////////// Whitening Transformation w  ////////////////////
	//ONE_DIM_ARRAY(W, double, (n*n));
	double *W;
	int r;
	withening_transformation(n, A, &W, &r , verbose) ;
	*rank=r;
	ONE_DIM_ARRAY(WtBW, double, (r*r));
	quadratic_Mat_multiply(B,  n  , W , r , WtBW ) ;
	//write_Fvec_matrix2file( WtBW, r, r, "WtBW.dat");
	ALLO_ONE_DIM_ARRAY( (*geigval),double,r) ;
	ONE_DIM_ARRAY(eigvec_b, double, (r*r));
	double_eigensolver_1D(	WtBW, r,  0, 	r,   *geigval,   eigvec_b);
	free(WtBW);
	ALLO_ONE_DIM_ARRAY( (*gevec),double,(n*r)) ;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, r, r, 1.0, W, n, eigvec_b, r, 0.0, *gevec, n);
	free(eigvec_b);
	free(W);
	/*********************************** for testing ***********************************/
	/************************************************************************************
	test if gevec can diagonalize A -> identity and B-> diagonal
	test if G^TBG == diag(gevals)
	************************************************************************************/
	if (verbose)
	{
		int i;
		printf("Rank of GE = %d \n", r);
		printf(" Nr. of  Dims = %d \n", n);
		printf("==========================\n");
		printf("Testing if the Simultaneous Diagonalization is fine:...\n");
		ONE_DIM_ARRAY(WtXW, double, (r*r));
		quadratic_Mat_multiply(A, n  ,  *gevec , r, WtXW );
		//for (i = 0; i <10 ; i++) printf("  %f ", WtXW[i*r +i]);
		int id = IsIdentity(WtXW,  r, 0.0000001, 0);
		if(! id)  fprintf(stderr,"Error: G^tAG != I\n");
		for (i = 0; i <(r*r); i++) WtXW[i]=0;
		quadratic_Mat_multiply(B, n  ,  *gevec , r , WtXW );
		int diag = IsDiagonal( WtXW,  r, 0.0000001 , 0) ;
		if(!diag) fprintf(stderr,"Error: G^tBG !=diag(geval) \n");
		int tt=1;
		for (i = 0; i <r; i++)  
		{	
			if(fabs( WtXW[i*r +i]-(*geigval)[i]) >  0.0000001)
			{
				printf("geval %d differs:  %f <-->  %f \n", i+1, WtXW[i*r +i],(*geigval)[i]  );
				tt=0;
			}
		}
		if (tt && diag && id) printf("Simultaneous Diagonalization is OK:\n G^t*A*G=I \n G^t*B*G=diag(geval) \n");
		 else fprintf(stderr,"Error:Simultaneous Diagonalization is NOT fine!\n");
		printf("==========================\n");
		free(WtXW);
	}
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///// complement projection on vector v
///// if ref is not null: projection w.r.t. ref (symmetric)
void Proj_on_vec(double *v, int n , double *mat ,double *proj, double *ref , int complement)
{
	int i,j;
	double scale ,*Pc,*Pct;
	ALLO_ONE_DIM_ARRAY(Pc, double, (n*n));
	for (i = 0; i <(n*n); i++) Pc[i]=0;
	if(ref)
	{
		scale = vecX_symat_vecY( ref, n  , v ,v);
	}
	else
	{
		scale= cblas_ddot (n, v, 1, v, 1);
	}
	if (complement) scale *= -1;
	cblas_dsyr( CblasColMajor ,CblasLower, n ,(1/scale), v ,1,Pc ,n );
	for (i=0; i< n; i++) for (j=i; j< n; j++) Pc[j*n +i]=Pc[i*n +j] ;
	if(ref)
	{
		ALLO_ONE_DIM_ARRAY(Pct, double, (n*n));
		cblas_dsymm(CblasColMajor, CblasRight, CblasLower, n, n, 1, ref, n, Pc, n, 0, Pct, n);
		for (i = 0; i < (n*n); i++) Pc[i] =Pct[i];
		
	} 
	if (complement) for (i = 0; i < n; i++) Pc[i*n +i] +=1;
	quadratic_Mat_multiply(mat,  n  , Pc , n , proj ) ;
	free(Pc);
	if(ref) free(Pct);
}
////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void Simultaneous_Diagonalization_subspacing(int n, double *A, double *B, double **geigval, double **gevec, int *rank ,double *dav  ,int verbose)
{
    //////////////// Whitening Transformation w  ////////////////////
	ONE_DIM_ARRAY(g_m, double, n);
	ONE_DIM_ARRAY(ginv,double,(n*n)) ;
	symmpinv(n, A, ginv);
	cblas_dgemv(CblasColMajor,CblasTrans,n, n, 1,ginv, n, dav, 1, 0, g_m, 1);
	free(ginv);
	double normalize=1/sqrt(vecX_symat_vecY( A, n , g_m  ,g_m   ));
	cblas_dscal(n, normalize, g_m, 1);
    int r ,i , j;
	/////////////////////////////////////////////////////////
	ONE_DIM_ARRAY(Agm, double, n);
    cblas_dsymv( CblasColMajor,  CblasLower, n, 1, A, n, g_m, 1, 0, Agm, 1);
    ONE_DIM_ARRAY(Ac, double, (n*n));
    Proj_on_vec(Agm, n , A ,Ac, NULL, 1);
    eigdecomposition_t eig_a;
    eigen(Ac, n,&eig_a);
    ///////////////////
    r= eig_a.rank;
    ONE_DIM_ARRAY(W, double, (n*r));
    double scale;
    for (i = 0; i <(eig_a.rank); i++)
    {
    	scale = 1/sqrt(eig_a.eval[i]);
    	for (j=0; j < n; j++) W[i*n +j]= scale* eig_a.vec[i*n +j];
    }
    free(eig_a.vec);
    free(eig_a.eval);
    //printf("W^t*Dav\n");
    // for (i = 0; i < 50; i++)  	printf( "%d  %f \n",i,cblas_ddot (n, W + i*n, 1,dav, 1) );
	/////////////////////////////////////////////////////////
	ONE_DIM_ARRAY(Bgm, double, n);
	cblas_dsymv( CblasColMajor,  CblasLower, n, 1, B, n, g_m, 1, 0,Bgm , 1);
	ONE_DIM_ARRAY(WBgm, double, r);
	cblas_dgemv(CblasColMajor,CblasTrans,n, r, 1, W, n, Bgm, 1, 0, WBgm, 1);
	free(Bgm);
    ONE_DIM_ARRAY(WtBW, double, (r*r));
    quadratic_Mat_multiply(B,  n  , W , r , WtBW ) ;
    ONE_DIM_ARRAY(WtBWc, double, (r*r));
    //residal_proj_mat(WtBWam, r , WtBW ,WtBWc, NULL, 1);
	Proj_on_vec(WBgm, r , WtBW ,WtBWc, NULL, 1);
	free(WtBW);
	free(WBgm);
    eigdecomposition_t eig_c;
    eigen(WtBWc, r,&eig_c);

	//printf(" Rank2 = %d, rank1= %d \n", eig_c.rank , r );
	free(WtBWc);
	*rank = eig_c.rank+1 ;
	r= *rank;

    ALLO_ONE_DIM_ARRAY( (*gevec),double,(n*r)) ;
    ALLO_ONE_DIM_ARRAY( (*geigval),double,r) ;
    lapack_ggemm_mat_prod( W, n, r,eig_c.vec ,eig_c.rank, *gevec+n);
	for (i = 0; i < n; i++) (*gevec)[i]= g_m[i];
    free(W);
 
    for (i = 0; i < r; i++)  (*geigval)[i]=vecX_symat_vecY( B, n  , *gevec +n*i ,*gevec +n*i );
    /****************
    printf("GtAG\n");
    for (i = 0; i < 30; i++)  printf( "%d  %f \n",i,vecX_symat_vecY( A, n  , *gevec +n*i ,g_m ));
    printf("GtBG\n");
    for (i = 0; i < 30; i++)  printf( "%d  %f \n",i,vecX_symat_vecY( B, n  , *gevec +n*i ,g_m ));
    /****************
	free(g_m);
	free(eig_c.vec);
	free(eig_c.eval);
    /*********************************** for testing ***********************************/
    /************************************************************************************
    test if gevec can diagonalize A -> identity and B-> diagonal
    test if G^TBG == diag(gevals)
    ************************************************************************************/
    if (verbose)
    {
        int i;
        printf("Rank of GE = %d \n", r);
        printf(" Nr. of  Dims = %d \n", n);
        printf("==========================\n");
        printf("Testing if the Simultaneous Diagonalization is fine: \n");
        //r= r-1;
        ONE_DIM_ARRAY(WtXW, double, (r* r));
        quadratic_Mat_multiply(A, n  ,  *gevec  , r, WtXW );
        //for (i = 0; i <10 ; i++) printf("  %f ", WtXW[i*r +i]);
        int id = IsIdentity(WtXW,  r, 0.0000001, 0);
        if(! id)  fprintf(stderr,"Error: G^tAG != I\n");
        for (i = 0; i <(r*r); i++) WtXW[i]=0;
        quadratic_Mat_multiply(B, n  ,  *gevec , r , WtXW );
        int diag = IsDiagonal( WtXW,  r, 0.0000001 , 0) ;
        if(!diag) fprintf(stderr,"Error: G^tBG !=diag(geval) \n");
        if ( diag && id) printf("Simultaneous Diagonalization is OK:\n G^t*A*G=I \n G^t*B*G=diag(geval) \n");
        else fprintf(stderr,"Simultaneous Diagonalization is NOT fine!\n");
        printf("==========================\n");
        free(WtXW);
    }
}
/////////////////////////////////////////////////////////////////////////
///////similar to the above ...not used
//////////////////////////////////////////////////////////////////////////
void Simultaneous_Diagonalization_last3(int n, double *A, double *B, double **geigval, double **gevec, int *rank ,double *dav  ,int verbose)
{
    //////////////// Whitening Transformation w  ////////////////////
    double *W1;
    int r ,i , j;
    withening_transformation(n, A, &W1, &r , verbose) ;
    *rank=r;
    ONE_DIM_ARRAY(a_m, double, r);
    cblas_dgemv(CblasColMajor,CblasTrans,n, r, 1, W1, n, dav, 1, 0, a_m, 1);
    cblas_dscal(r, (1/cblas_dnrm2(r,a_m,1)), a_m, 1);
    ONE_DIM_ARRAY(g_m, double, n);
    cblas_dgemv(CblasColMajor,CblasNoTrans,n, r, 1, W1, n, a_m, 1, 0, g_m, 1);
    double m= cblas_ddot (n, g_m, 1, dav, 1);
    double aa = vecX_symat_vecY( B, n , g_m  ,g_m   );
    printf( "Dav_kl_var= %f\n", 0.5* (-log(aa) + aa -1 ));
    printf( "Dav_kl= %f\n", 0.5*m*m);
	/////////////////////////////////////////////////////////
	ONE_DIM_ARRAY(Agm, double, n);
    cblas_dsymv( CblasColMajor,  CblasLower, n, 1, A, n, g_m, 1, 0, Agm, 1);    
    ONE_DIM_ARRAY(Ac, double, (n*n));
    Proj_on_vec(Agm, n , A ,Ac, NULL, 1);
    eigdecomposition_t eig_a;
    eigen(Ac, n,&eig_a);
	ONE_DIM_ARRAY(G0tAG0, double, (eig_a.rank*eig_a.rank)); 
	quadratic_Mat_multiply(A,  n  , eig_a.vec , eig_a.rank , G0tAG0 ) ;
	double *W2;
	withening_transformation(eig_a.rank, G0tAG0, &W2, &r , verbose) ;
	ONE_DIM_ARRAY(W, double, (n*r));
	lapack_ggemm_mat_prod( eig_a.vec, n, eig_a.rank, W2,r, W);
	/////////////////////////////////////////////////////////   
	ONE_DIM_ARRAY(Bgm, double, n);
	cblas_dsymv( CblasColMajor,  CblasLower, n, 1, B, n, g_m, 1, 0,Bgm , 1);
	ONE_DIM_ARRAY(WBgm, double, r);
	cblas_dgemv(CblasColMajor,CblasTrans,n, r, 1, W, n, Bgm, 1, 0, WBgm, 1);
    ONE_DIM_ARRAY(WtBW, double, (r*r));
    quadratic_Mat_multiply(B,  n  , W , r , WtBW ) ;
    ONE_DIM_ARRAY(WtBWc, double, (r*r));
    //residal_proj_mat(WtBWam, r , WtBW ,WtBWc, NULL, 1);
	Proj_on_vec(WBgm, r , WtBW ,WtBWc, NULL, 1);
    eigdecomposition_t eig_c;
    //eig_c.mat=WtBWc;
    eigen(WtBWc, r,&eig_c);
	printf(" Rank2 = %d, rank1= %d \n", eig_c.rank , r );
	free(WtBWc);
	*rank = eig_c.rank+1 ;
	r= *rank;
    free(WtBW);
    ALLO_ONE_DIM_ARRAY( (*gevec),double,(n*r)) ;
    ALLO_ONE_DIM_ARRAY( (*geigval),double,r) ;
    lapack_ggemm_mat_prod( W, n, r,eig_c.vec ,eig_c.rank, *gevec+n);
	for (i = 0; i < n; i++) (*gevec)[i]= g_m[i]; 
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, r, r, 1.0, W, n, G1, r, 0.0, *gevec, n);
     free(W);
//free(G1);
/*****************/
 
    for (i = 0; i < r; i++)  (*geigval)[i]=vecX_symat_vecY( B, n  , *gevec +n*i ,*gevec +n*i ); 
    for (i = 0; i < 30; i++)  printf( "%d  %f \n",i,vecX_symat_vecY( A, n  , *gevec +n*i ,g_m )); 
    printf("BBB\n");
    for (i = 0; i < 30; i++)  printf( "%d  %f \n",i,vecX_symat_vecY( B, n  , *gevec +n*i ,g_m )); 
	free(a_m);
	free(g_m);
    /*********************************** for testing ***********************************/
    /************************************************************************************
    test if gevec can diagonalize A -> identity and B-> diagonal
    test if G^TBG == diag(gevals)
    ************************************************************************************/
    if (verbose)
    {
        int i;
        printf("Rank of GE = %d \n", r);
        printf(" Nr. of  Dims = %d \n", n);
        printf("Testing if the Simultaneous Diagonalization is fine: \n");
        //r= r-1;
        ONE_DIM_ARRAY(WtXW, double, (r* r));
        quadratic_Mat_multiply(A, n  ,  *gevec  , r, WtXW );
        //for (i = 0; i <10 ; i++) printf("  %f ", WtXW[i*r +i]);
        printf("\n");
        int id = IsIdentity(WtXW,  r, 0.0000001, 0);
        printf("WtAW= %d\n", id);
        for (i = 0; i <(r*r); i++) WtXW[i]=0;
        quadratic_Mat_multiply(B, n  ,  *gevec , r , WtXW );
        int diag = IsDiagonal( WtXW,  r, 0.0000001 , 0) ;
        printf("WtBW= %d\n", diag);
        if ( diag && id) printf("Simultaneous Diagonalization is fine!\n G^t*A*G=I \n G^t*B*G=diag(geval) \n");
        else printf("Simultaneous Diagonalization is NOT fine!\n");
        free(WtXW);
    }
}
 
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////  projection matrix on the subspace spanned by the vectors of a matrix G (r x c)
///// P = G(G^tG)^-1 G^t
void Projection_matrix(int r, int c , double *G ,double *P)
{
		ONE_DIM_ARRAY(GtG, double, (c*c));
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, c, c, r, 1.0, G, r, G, r, 0.0, GtG, c);
		ONE_DIM_ARRAY(GtG_inv, double, (c*c));
		symmpinv(c, GtG, GtG_inv) ;
		free(GtG);
		ONE_DIM_ARRAY(G_GtG_inv, double, (r*c));
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, r, c, c, 1.0, G, r, GtG_inv, c, 0.0, G_GtG_inv, r);
		free(GtG_inv);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, r, c, r, 1.0, G_GtG_inv, r, G, r, 0.0, P, r);
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows_A, ncol_B, ncol_A, 1.0, A, nrows_A, B, ncol_A, 0.0, C, nrows_A);
		free(G_GtG_inv);
}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
