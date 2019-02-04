#include <cblas.h>
#include <math.h>
#include "covariance.h" 
#include "arrays.h"
#include "fitting.h"
#include "xdrfile_xtc.h"
#include "linear_algebra.h"
#include "strTools.h"


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
void blas_compute_covariance
			(const char 	*trxfile,   // traj file to read
			rvec  		*xref,      // reference coor structure to be used for fitting 
			double 		*mat,  // pointer where the cov mat is allocated using 1D array
			int 		nfit,   // nr atoms used for fitting
			atom_id 	*ifit,  // index of fitting atoms
			int 		natom,   // nr atoms used for used for covariance analysis
			atom_id 	*index,  // index of atoms used for covariance analysis
			double     *m ,   /// array of the average atoms
			float first_time,
			float last_time
		     )
{

    double            inv_nframes;
    real            t;
    int             nframes;
    int ndim, i, j, d;
    nframes = 0;
    ndim = natom*3;
    fprintf(stderr, "\nConstructing covariance matrix (%dx%d) ...\n", (int)ndim, (int)ndim);
    xtc_t f;
	Read_first_xtc( trxfile,&f, first_time );
	//////// make a reference of selected fitting atoms
	rvec *ref_sel;

	if (xref && (nfit < f.natoms))
	{

		ALLO_ONE_DIM_ARRAY( ref_sel, rvec, nfit);
		for (i = 0; i < nfit; i++) for (d=0; d < 3; d++) ref_sel[i][d]= xref[ifit[i]][d];
	}
	else
	{
	  ref_sel= xref;
	}
	rvec cog;
	if (xref) get_COG(ref_sel ,NULL,nfit,  cog );

	for (i=0; i< (ndim*ndim) ; i++)  mat[i]=0;
	for(i=0; i < ndim; i++) m[i]=0;
	ONE_DIM_ARRAY( xt, double, ndim);
	int ndx;	
    do
    {
        nframes++;
        if(xref)
		{
			//bb=calc_rmsd( xread, ifit, xref, ifit , nfit);
			//center_x(nfit, ifit, nat,  xread);
			//rotate_x(nat,  xref, xread, nfit, ifit);
			//center_x2ref_diffrenct_atoms(f.x, ifit, ref_sel, NULL , nfit, f.natoms );
        	center_x(f.x, ifit,  nfit,f.natoms, cog );
			rotate_B2R(                  f.x, ifit, ref_sel, NULL , nfit, f.natoms );
			//ee=calc_rmsd( xread, ifit, xref, ifit , nfit);
			//printf("RMSD change: %g --> %g \n",bb,ee);
		}
		if(index)
		{
			for (i = 0; i < natom; i++)
			{
				ndx=index[i];
				xt[i*3  ]=(double) f.x[ndx][0] ;
				xt[i*3+1]=(double) f.x[ndx][1] ;
				xt[i*3+2]=(double) f.x[ndx][2] ;
			}
		}
		else
		{
			for (i = 0; i < natom; i++) 
			{
				xt[i*3  ]=(double) f.x[i][0] ;
				xt[i*3+1]=(double) f.x[i][1] ;
				xt[i*3+2]=(double) f.x[i][2] ;
			}
		}
		cblas_daxpy (ndim, 1, xt, 1,m, 1);
		// SSYR   performs the symmetric rank 1 operation    A := alpha*x*x**T + A,
		//  where alpha is a double scalar, x is an n element vector and A is an  n by n symmetric matrix.
		cblas_dsyr( CblasColMajor ,CblasLower, ndim ,1, xt ,1,mat ,ndim );
	}
    while( Read_next_xtc(&f, last_time));
    reset_frame( &f);
    inv_nframes = 1.0/ ((double) nframes);
	cblas_dscal(ndim, inv_nframes, m, 1);
	/// ma = inv_nframes* sum(xx^T) - xavfit*xavfit^T
	cblas_dscal(ndim*ndim , inv_nframes, mat, 1);
	cblas_dsyr( CblasColMajor ,CblasLower, ndim ,-1, m ,1,mat ,ndim );
	//// fill the upper part of the sym mat
    for (i=0; i< ndim; i++) for (j=i; j< ndim; j++) mat[j*ndim +i]=mat[i*ndim +j] ;
	free(xt);
	printf("*********************************************\n");
	printf("Computation of the covariance matrix is done!\n");
	printf("Number of dimensions: %d\n", ndim);
	printf("Number of used data points (conformations): %d\n", nframes);
	if(nframes < ndim)
	{
		fprintf(stderr,"Note!!!:Number of used data points (%d) is less than\nthe number of dimensions(%d)\n", nframes, ndim);
		fprintf(stderr,"The rank of the matrix will be reduced by  %d \n", nframes- ndim);
	}
	printf("*********************************************\n");
}

//////////////////////////////////////////////////////////////////
void  cov2datfile(double *mat, int ndim, const char *covmatfile)
{
	printf("Saving the covariance matrix to the file %s\n", covmatfile);
	int i,j;
	FILE           *out;
	out = fopen(covmatfile, "w");

	fprintf(out, "##### covariance of the atoms  (x1y1z1x2y2z2...)\n");
	for (j = 0; j < ndim; j++)
	{
		for (i = 0; i < ndim; i++)	fprintf(out, "%g \t", mat[ndim*j+i]);
		fprintf(out, " \n");
	}
	fclose(out);
}

//////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/********************************************************************
read and write average and the symm. (covariance) matrix to a binary file
note: the structure of xdr defines reading and writing which
have the same function inside xdr routines. 
*********************************************************************/
 int write_read_cov_matrix(const char *outfile,double **mat , double **average,int *ndim, int W)
{
	XDRFILE *xd;
	if (W) xd = xdrfile_open(outfile,"w");
	else   xd = xdrfile_open(outfile,"r");
	if(NULL == xd) PRINT_ERR_OPENFILE(outfile);
    if (xdrfile_read_int(ndim,1,xd) != 1)
		return exdrINT;
	int n= *ndim;
	if(!W)
	{
		ALLO_ONE_DIM_ARRAY((*average), double, n);
		ALLO_ONE_DIM_ARRAY((*mat), double, (n*n));
	}
	if (xdrfile_read_double((*average),n,xd) != n)	return exdrDOUBLE;
	if (xdrfile_read_double((*mat),n*n,xd) != n*n)	return exdrDOUBLE;
	xdrfile_close(xd);
	return exdrOK;
}

////////////////////////////////////////////////////////////////////////////

int write_cov_matrix(const char *outfile,double **mat , double **average,int *ndim,int verbose)
{
	int out;
	if ( exdrOK != write_read_cov_matrix(outfile ,mat ,average,ndim, 1 ))
		out=0;
	else out=1;
	int i, n2 , n=*ndim;
	double *mat2, *av2;
	if(verbose > 1)
	{
		printf("**** Testing writing the cov. matrix via rereading...\n");
		if (  exdrOK != write_read_cov_matrix(outfile ,&mat2 ,&av2 ,&n2, 0 ) ) fprintf(stderr,"Error in re-reading the covariance matrix\n");
		if( n!= n2) printf("different ndim when rereading the covariance %d --- %d\n", n, n2);
		double tol=0.0000000001;
		for(i=0;i< n; i++) if( ((*average)[i]-av2[i] ) > tol) printf("different read- write average: %g--- %g \n", (*average)[i],av2[i]);
		for(i=0;i< n*n; i++) if( ((*mat)[i]-mat2[i] ) > tol) printf("different read- write covariance: %g--- %g \n", (*mat)[i],mat2[i]);
		free(mat2);
		free(av2);
	}
	return out;
}

////////////////////////////////////////////////////////////////////////////

int read_cov_matrix(const char *outfile,double **mat , double **average,int *ndim)
{
	if ( exdrOK != write_read_cov_matrix(outfile ,mat ,average ,ndim, 0 ))
		return 0;
	else return 1;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**************************************************************************
read and write eigen value and vectors and mean from/to a binary file *.mvv
note: the structure of xdr defines reading and writing which
have the same function inside xdr routines. 
****************************************************************************/
 int write_read_eigen(const char *outfile,eigdecomposition_t *eig , double **average, int W)
{
	XDRFILE *xd;
	if (W) xd = xdrfile_open(outfile,"w");
	else   xd = xdrfile_open(outfile,"r");
	if(NULL == xd) PRINT_ERR_OPENFILE(outfile);
    if (xdrfile_read_int(&(eig->ndim),1,xd) != 1)
		return exdrINT;
	int n= eig->ndim;
	if(!W)
	{
		ALLO_ONE_DIM_ARRAY( (*average), double,n);
		ALLO_ONE_DIM_ARRAY( (eig->eval), double,n);
		ALLO_ONE_DIM_ARRAY(eig->vec, double, (n* n));
	}
	if (xdrfile_read_double( (*average) ,n,xd) != n)
            return exdrDOUBLE;
    if (xdrfile_read_double(eig->eval,n,xd) != n)
            return exdrDOUBLE;
    if (xdrfile_read_double(eig->vec,n*n,xd) != n*n )
			return exdrDOUBLE;
	xdrfile_close(xd);
	return exdrOK;
}
/**************************************/
int write_eigen(const char *outfile,eigdecomposition_t *eig , double **average, int verbose)
{
	int res= 1;
	if( exdrOK != write_read_eigen(outfile,eig , average, 1)) res=0 ;
	int i;
	double tol= 0.000000001;
	if (verbose > 1)//// test
	{
		printf("Testing writting eigenpairs via rereading\n");
		ONE_DIM_ARRAY( eig2, eigdecomposition_t,1);
		double *m2;
		write_read_eigen("eigen.mvv",eig2 , &m2 ,0);
		if(eig->ndim != eig2->ndim) printf("Eroor: diff. ndim\n");
		int n= eig->ndim;
		for( i=0 ; i < n; i++) if ((m2[i] - (*average)[i]) > tol ) printf("Eroor: diff. averages\n");
		//for( i=0 ; i < 10; i++)  printf("averages diff. %g \n", m2[i] - (*average)[i] );
		for( i=0 ; i < n; i++) if ((eig2->eval[i] - eig->eval[i]) > tol ) printf("Eroor: diff. evalues\n");
		//for( i=10 ; i < 30; i++)  printf("evals diff. %g - %g = %g \n", eig2->eval[i] , eig->eval[i], eig2->eval[i] - eig->eval[i] );
		for( i=0 ; i < n*n; i++) if ((eig2->vec[i] - eig->vec[i]) > tol ) printf("Eroor: diff. %d evector is diff. %g --- %g\n",i, eig2->vec[i] , eig->vec[i] );
		free(eig2);
		free(m2);
	}
	return res;
}
///////////////////////////////////////////////
int read_eigen(const char *outfile,eigdecomposition_t *eig , double **average)
{
	if( exdrOK != write_read_eigen(outfile,eig , average, 0)) return 0;
	else return 1;
}
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************
read and write mean of state A , B ,generalized eigenvalues and vectors from/to a
binary file *.mmvv
note: the structure of xdr defines reading and writing which
have the same function inside xdr routines. 
************************************************************************************/
 int write_read_Geigen(const char *outfile,Geigdecomposition_t *eig ,int W)
{
	XDRFILE *xd;
	if (W) xd = xdrfile_open(outfile,"w");
	else   xd = xdrfile_open(outfile,"r");
	if(NULL == xd) PRINT_ERR_OPENFILE(outfile);
    if (xdrfile_read_int(&(eig->ndim),1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&(eig->rank),1,xd) != 1)
		return exdrINT;
	int n = eig->ndim;
	int r = eig->rank;
	if(!W)
	{
		ALLO_ONE_DIM_ARRAY( eig->m_a, double,n);
		ALLO_ONE_DIM_ARRAY( eig->m_b, double,n);
		ALLO_ONE_DIM_ARRAY( eig->kl, double,r);
		ALLO_ONE_DIM_ARRAY( eig->kl_m, double,r);
		ALLO_ONE_DIM_ARRAY( eig->acc_kl, double,r);
		ALLO_ONE_DIM_ARRAY( eig->geigval, double,r);
		ALLO_ONE_DIM_ARRAY(eig->gevec, double, (n* r));
	}
	if (xdrfile_read_double( eig->m_a ,n,xd) != n)
            return exdrDOUBLE;
	if (xdrfile_read_double( eig->m_b ,n,xd) != n)
            return exdrDOUBLE;
	if (xdrfile_read_double(eig->kl,r,xd) != r)
            return exdrDOUBLE;
	if (xdrfile_read_double(eig->kl_m,r,xd) != r)
	            return exdrDOUBLE;
	if (xdrfile_read_double(eig->acc_kl,r,xd) != r)
	            return exdrDOUBLE;
    if (xdrfile_read_double(eig->geigval,r,xd) != r)
            return exdrDOUBLE;
    if (xdrfile_read_double(eig->gevec,n*r,xd) != n*r )
			return exdrDOUBLE;
	xdrfile_close(xd);
	return exdrOK;
}
/////////////////////////////////////////////////////////////////////////////////////

int read_Geigen(const char *outfile,Geigdecomposition_t *eig )
{
	int out, result= 1;
	out = write_read_Geigen(outfile, eig ,0);
	if( exdrOK != out)
	{
		printf("xdr error = %d\n",out );
		result = 0;
	}
	return result;
}

/////////////////////////////////////////////////////////////////////
int write_Geigen(const char *outfile,Geigdecomposition_t *eig , int verbose)
{
	int aa, result= 1;
	aa= write_read_Geigen(outfile, eig ,1);
	if( exdrOK != aa)
	{
		printf("exdr error = %d\n",aa );
		result = 0;
	}
	int i, j;
	double tol= 0.000000001;
	//// for testing 
	int out=1;
	if (verbose > 1)
	{
		printf("Testing writing of Geigen via rereading...\n");
		ONE_DIM_ARRAY( eig2, Geigdecomposition_t ,1); 
		if (!read_Geigen(outfile,eig2 ) )
		{
			fprintf(stderr," Can not read G eigenpairs from the file %s \n",outfile );
			out= 0;
		}
		if(eig->ndim != eig2->ndim)
		{
			fprintf(stderr,"Error: diff. ndim\n");
			out= 0;
		}
		if(eig->rank != eig2->rank)
		{
			fprintf(stderr,"Error: diff. rank\n");
			out= 0;
		}
		for( i=0 ; i < eig2->ndim; i++) if ((eig->m_a[i] - eig2->m_a[i]) > tol )
		{
			fprintf(stderr,"reading error: diff. mean A\n");
			out= 0;
		}
		for( i=0 ; i < eig2->ndim; i++) if ((eig->m_b[i] - eig2->m_b[i]) > tol )
		{
			fprintf(stderr,"reading error: diff. mean B\n");
			out= 0;
		}
		for( i=0 ; i < eig2->rank; i++) if ((eig->kl[i] - eig2->kl[i]) > tol )
		{
			fprintf(stderr,"reading error: diff. KL values B\n");
			out= 0;
		}
		for( i=0 ; i < eig2->rank; i++) if ((eig->geigval[i] - eig2->geigval[i]) > tol )
		{
			fprintf(stderr,"reading error: diff. geigval B\n");
			out= 0;
		}
		for( i=0 ; i < (eig2->rank * eig2->ndim); i++) if ((eig->gevec[i] - eig2->gevec[i]) > tol )
		{
			fprintf(stderr,"reading error: diff. gevec B\n");
			out= 0;
		}
		if(!out) fprintf(stderr,"rereading problem!\n");
		else fprintf(stderr,"Rereading is fine.\n");
	}	
	return result;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void write_array2xvg(const char *eigvalfile, int ndim, double *eigval, char *title ,char *xaxis, char *yaxis , char *header )
{		
	//fprintf(stderr,"\nWriting eigenvalues to %s\n",eigvalfile);

	FILE   *out;
	out=fopen(eigvalfile,"w");
    //fprintf(out,"# by the following command:\n# %s\n#\n",command_line());
	fprintf(out,"# by the following command:\n# %s\n#\n",header);
	fprintf(out,"@    title \"%s\"\n", title);
	fprintf(out,"@    xaxis  label \"%s\"\n",xaxis );
	fprintf(out,"@    yaxis  label \"%s\"\n",yaxis);
    fprintf(out,"@TYPE xy\n");
   	int i;
	for (i=0; (i<ndim); i++)  fprintf (out,"%10d %g\n",i+1,eigval[i]);
	fclose(out); 
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
 void sort_decreasing( double *x , int n,  int *x_indx)
{
	int i ,j , ff;
	double tmp;
	for (i = 0; i < n; ++i) x_indx[i]=i; 
	for (i = 0; i < n; ++i)
    {
	    tmp= x[x_indx[i]];
		for (j = i + 1; j < n; ++j)
        {
            if ( tmp < x[x_indx[j]])
            {
                ff= x_indx[i];
				x_indx[i] = x_indx[j];
				x_indx[j] =  ff;
				tmp= x[x_indx[i]];
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
double max_array(double *a, int n)
{
	int i;
	double max=a[0];
	for (i = 1; i < n; ++i) if(a[i]> max) max=a[i];
	return max;
}
///////////////////////////////////////////////////////////////////////////////
void relative_pca(Geigdecomposition_t *G,int ordered ,int algo,int verbose  )
{
	int ndim = G->ndim;
	int i,v, ndx;
	if (verbose)
	{
		printf("*********************************************\n");
		printf("Running the Simultaneous diagonalization...\n");
		printf("*********************************************\n");
	}
	ONE_DIM_ARRAY(dav, double , ndim );
	for (i = 0; i < ndim; i++) dav[i]=   G->m_b[i] - G->m_a[i] ;
	if(algo==0)
	{
		printf("RPCA will be performed using algorithm 0 (no sub-spacing)\n");
		Simultaneous_Diagonalization(ndim, G->mat_a, G->mat_b, &(G->geigval), &(G->gevec), &(G->rank),verbose);
	}
	else if (algo==1)
	{
		printf("RPCA will be performed using algorithm 1 (average variance sub-spacing)\n");
		Simultaneous_Diagonalization_subspacing(ndim, G->mat_a, G->mat_b, &(G->geigval), &(G->gevec), &(G->rank),dav,verbose);
	}
	int rank =G->rank;
	ALLO_ONE_DIM_ARRAY((G->kl), double , rank );
	ALLO_ONE_DIM_ARRAY((G->acc_kl), double , rank );
	ALLO_ONE_DIM_ARRAY((G->kl_m), double , rank );
	double *pv , m ;
	for (i = 0; i < rank; i++)
	{
		m=0;
		pv= G->gevec +i* ndim;
		m=cblas_ddot(ndim, pv, 1, dav, 1);
		G->kl_m[i]= 0.5*m*m;
		G->kl[i] = 0.5*( -log(G->geigval[i] ) + G->geigval[i]-1 ) + G->kl_m[i];
	}
	G->sum_kl= 0 , G->sum_kl_m=0;
	for (i = 0; i < rank; i++)
	{
		G->sum_kl   += G->kl[i];
		G->sum_kl_m += G->kl_m[i];
	}
	if (verbose)
	{
		printf("total DKL= %f\n",G->sum_kl );
		printf("total DKL due to variance change = %f\n",G->sum_kl - G->sum_kl_m);
		printf("total Mean change contr. to DKL = %f\n",G->sum_kl_m );
	}
	///////////////
	int *klInd;
	double *tmp_sorted;
	if (ordered) 
	{
		G->ordered=1;
		printf("*********************************************\n");
		printf("Reordering according to the values of DKL\n" );
		printf("*********************************************\n");
		ALLO_ONE_DIM_ARRAY(klInd, int ,rank );
		sort_decreasing( G->kl ,rank ,klInd);
		ALLO_ONE_DIM_ARRAY(tmp_sorted, double , rank );
		for (i = 0; i <rank; i++) tmp_sorted[i]= G->kl[klInd[i]];
		for (i = 0; i <rank; i++) G->kl[i]=tmp_sorted[i];
		for (i = 0; i <rank; i++) tmp_sorted[i]= G->kl_m[klInd[i]];
		for (i = 0; i <rank; i++) G->kl_m[i]=tmp_sorted[i];
		for (i = 0; i <rank; i++) tmp_sorted[i]= G->geigval[klInd[i]];
		for (i = 0; i <rank; i++) G->geigval[i]=tmp_sorted[i];
		free(tmp_sorted);
		ALLO_ONE_DIM_ARRAY(tmp_sorted, double , (ndim*rank) );
		for (i = 0; i < (ndim*rank); i++) tmp_sorted[i]=G->gevec[i] ;
		for (v = 0; v <rank; v++) 
		{
			ndx= klInd[v];
			for (i = 0; i < ndim; i++) G->gevec[v*ndim +i]= tmp_sorted[ndx*ndim +i];
		}
		free(klInd);
		free(tmp_sorted);
		free(dav);
	}
	G->acc_kl[0] = G->kl[0];
	for (i = 1; i < rank; i++) G->acc_kl[i]= G->acc_kl[i-1] + G->kl[i] ;
	for (i = 0; i < rank; i++) G->acc_kl[i]= G->acc_kl[i] *(100/ G->sum_kl) ;
	printf("*********************************************\n");
	printf("Simultaneous diagonalization is done!\n" );
	printf("*********************************************\n");
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void write_Geigen2xvg(const char *Geigfile, Geigdecomposition_t *eig ,  char *header)
{		
	FILE   *out;
	out=fopen(Geigfile,"w");
	fprintf(out,"# by the following command:\n# %s\n#\n", header);
    //fprintf(out,"# by the following command:\n# %s\n#\n",command_line());
	fprintf(out,"@    title \"Generalized Eigen-Analysis\"\n");
	fprintf(out,"@    xaxis  label \"Generalized Eigen Component\"\n" );
	fprintf(out,"@    yaxis  label \"Divergence(kt)\"\n");
    fprintf(out,"@TYPE xy\n");
	fprintf(out,"@ view 0.15, 0.15, 0.75, 0.85\n");
	fprintf(out,"@ legend on\n");
	fprintf(out,"@ legend box on\n");
	fprintf(out,"@ legend loctype view\n");
	fprintf(out,"@ legend 0.78, 0.8\n");
	fprintf(out,"@ legend length 2\n");
	fprintf(out,"@ s0 legend \"Total\"\n");
	fprintf(out,"@ s1 legend \"Due to mean change\"\n");
	fprintf(out,"@ s2 legend \"Accumulated (percent)\"\n");
	fprintf(out,"@ s3 legend \"Generalized eigenvalues (nm)\"\n");
   	int i;
	for (i=0; i< eig->rank; i++)  fprintf (out,"%10d	%12.3f	%12.3f	%12.3f	%e \n",
									            i+1,eig->kl[i],eig->kl_m[i], eig->acc_kl[i], eig->geigval[i]);
	fclose(out); 
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void RPCA_project(const char 	*trxfile,   // traj file to read
					int          state_A,   /// project state A or B
					const char *projfile ,    // write the projections to a dat file
					Geigdecomposition_t *eig, 
					int Beg_evec,
					int End_evec,
					int 		natom,
					atom_id 	*index,
					float first_time,
					float last_time,
					int skip           /// nr of frames to be skipped between projected frames
					)
{
	
	int ndim= eig->ndim;
	if(natom != ((int)(ndim/3)))
	{
		PRINT_FATAL_ERR_MES("Nr of atoms in the index is different from ndim/3 in Geigen structure");
	}
	printf("Projection of the conformations in file %s on the Generalized eigenvectors %d to %d...\n", trxfile,Beg_evec, End_evec );
	printf("saving the projections to the file %s\n", projfile);
	FILE    *out;
	out = fopen(projfile, "w");
	fprintf(out, "##### Projection on PCAs %d to %d of the conformations in %s (x1y1z1x2y2z2...)\n",Beg_evec+1,End_evec, trxfile);
	int   v, a, d ,i, j;
	rvec *ref ;
	ONE_DIM_ARRAY(ref_a, rvec ,natom);
	ONE_DIM_ARRAY(ref_b, rvec ,natom);
	for (i = 0; i < natom; i++) 
	{
		for (j=0; j< 3 ;j++)
		{ 	
			ref_a[i][j]= 	(real) eig->m_a[i*3 +j] ;
			ref_b[i][j]= 	(real) eig->m_b[i*3 +j] ;
		}
	}
	if (state_A) ref= ref_a;
	else        ref= ref_b;
    int             nframes0;
    xtc_t f;
    Read_first_xtc( trxfile,&f, first_time );
    int nr_vec= End_evec-Beg_evec;
	ONE_DIM_ARRAY(progection, real, (End_evec-Beg_evec));
	ONE_DIM_ARRAY(x, double ,ndim);
	int frame_count=0;
	int co;
    do
    {
    	frame_count++;
    	if(skip && frame_count < skip)
    	{
    		continue;
    	}
    	else if (skip && frame_count== skip)
    	{
    		frame_count=0;
    	}

		center_x2ref_diffrenct_atoms( f.x, index, ref, NULL , natom, f.natoms );
		rotate_B2R(                   f.x, index, ref, NULL , natom, f.natoms );
		
		for (i = 0; i < natom; i++)
		{
			for(d=0; d< 3; d++)  x[i*3 +d]= (double)(f.x[index[i]][d] - eig->m_b[i*3 +d]);
		}
		//ZERO_ONE_DIM_ARRAY(progection,(End_evec-Beg_evec));
		co=0;
		for(v= Beg_evec; v < End_evec; v++)
		{
			progection[co]= cblas_ddot (ndim, x, 1, eig->gevec + v*ndim, 1);
			co +=1;
		}
		fprintf(out, "%12.3f   ", f.time);
		for(i=0; i<(End_evec-Beg_evec) ;i++ ) fprintf(out, "%.3f ", progection[i] );
		fprintf(out, " \n");
    }
    while( Read_next_xtc(&f, last_time));
    reset_frame( &f);
    fclose(out);
	free(progection);
	free(x);
	free(ref_a);
	free(ref_b);
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void PCA_project(const char 	*trxfile,   // traj file to read
					const char *projfile ,    // write the projections to a dat file
					eigdecomposition_t *eig,
					double             *m,   //  the average
					int Beg_evec,
					int End_evec,
					int 		natom,
					atom_id 	*index,
					float first_time,
					float last_time
					)
{

	int ndim= eig->ndim;
	if(natom != ((int)(ndim/3)))
	{
		PRINT_FATAL_ERR_MES("Nr of atoms in the index is different from ndim/3 in eigen structure");
	}
	printf("Projection of the conformations in file %s on the eigenvectors %d to %d...\n", trxfile,Beg_evec, End_evec );
	printf("saving the projections to the file %s\n", projfile);
	FILE    *out;
	out = fopen(projfile, "w");
	fprintf(out, "##### Projection on PCAs %d to %d of the conformations in %s (x1y1z1x2y2z2...)\n",Beg_evec+1,End_evec, trxfile);
	int   v, a, d ,i, j;
	ONE_DIM_ARRAY(ref, rvec ,natom);
	for (i = 0; i < natom; i++) for (j=0; j< 3 ;j++) ref[i][j]= 	(real) m[i*3 +j] ;

    xtc_t f;
    Read_first_xtc( trxfile,&f, first_time );
    int nr_vec= End_evec-Beg_evec;
	ONE_DIM_ARRAY(progection, real, (End_evec-Beg_evec));
	ONE_DIM_ARRAY(x, double ,ndim);
	int co;
    do
    {
		center_x2ref_diffrenct_atoms( f.x, index, ref, NULL , natom, f.natoms );
		rotate_B2R(                   f.x, index, ref, NULL , natom, f.natoms );

		for (i = 0; i < natom; i++)
		{
			for(d=0; d< 3; d++)  x[i*3 +d]= (double)(f.x[index[i]][d] - m[i*3 +d]);
		}
		//ZERO_ONE_DIM_ARRAY(progection,(End_evec-Beg_evec));
		co=0;
		for(v= Beg_evec; v < End_evec; v++)
		{
			progection[co]= cblas_ddot (ndim, x, 1, eig->vec + v*ndim, 1);
			co +=1;
		}
		fprintf(out, "%12.3f   ", f.time);
		for(i=0; i<(End_evec-Beg_evec) ;i++ ) fprintf(out, "%.3f ", progection[i] );
		fprintf(out, " \n");
    }
    while( Read_next_xtc(&f, last_time));
    reset_frame( &f);
    fclose(out);
	free(progection);
	free(x);
	free(ref);
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void interaction_map(const char *resfile, Geigdecomposition_t *eig, t_atoms *atoms,atom_id *index,int Beg_evec,int	End_evec, const char *pdbBfactorsfile)
{
	printf("The interaction map will be constructed using the generalized eigenvectors %d to %d...\n", Beg_evec, End_evec );

	int ndim= eig->ndim ;
	int natom = (int) ( ndim /3);
	int i ,j ;
	////// count nr of residues in the atoms
	int nres=1 ;
	for (i=1; i< natom; i++ )	if (atoms->atom[ index[i] ].resind != atoms->atom[ index[i-1] ].resind ) nres +=1;
	//////////////dissection the atoms accor. to their residues
	typedef struct
	{
		int 	 index;    /// index of the fist atom in the residue
		int 	natoms;			// nr. of atoms
		int    	resi;        // residue nr in the pdb file
		char resname[5];
		char   chain;
	} resindx_t ;

	ONE_DIM_ARRAY(residues, resindx_t, nres );
	int  r=-1;
	int curr_res,old_res=-1;
	for (i=0; i< natom; i++ )
	{
		curr_res= atoms->atom[ index[i] ].resind;
		if (curr_res != old_res )
		{
			r++;
			residues[r].natoms=0;
			residues[r].index = i;
			residues[r].resi=atoms->resinfo[curr_res].nr  ;
			residues[r].chain=atoms->resinfo[curr_res].chainid;
			strcpy(residues[r].resname ,  *(atoms->resinfo[curr_res].name));
			old_res=curr_res;
		}
		residues[r].natoms +=1;
	}
	///////////////////// matrix for the average change
	ONE_DIM_ARRAY(dm, double,ndim );
	for (i = 0; i < ndim; i++) dm[i] = eig->m_b[i]-  eig->m_a[i];
	ONE_DIM_ARRAY(MM, double,ndim*ndim );
	for (i = 0; i < (ndim*ndim); i++) MM[i]=0;
	cblas_dsyr( CblasColMajor ,CblasLower, ndim ,1, dm ,1,MM ,ndim );
	//// fill the upper part of the sym mat
	for (i=0; i< ndim; i++) for (j=i; j< ndim; j++) MM[j*ndim +i]=MM[i*ndim +j] ;
	//////////////////////////
	ONE_DIM_ARRAY(contributions, double,(nres * nres) );
	for (i = 0; i < (nres *nres); i++) contributions[i]=0;
	int r1,r2 ,p1,p2, n1,n2;
	int v=0;
	for (r1=0; r1< nres; r1++ )
	{
		n1=residues[r1].natoms *3;
		p1= residues[r1].index *3 ;
		for (r2=r1; r2< nres; r2++ )
		{
			n2=residues[r2].natoms *3;
			p2= residues[r2].index *3 ;
			for (v = Beg_evec; v < End_evec; v++)
			{
				for (i = 0; i < n2; i++)
				{
					contributions[r1*nres + r2] += cblas_ddot (n1, eig->gevec+ v*ndim +p1, 1,  eig->mat_b +p1 +ndim*(p2+i), 1) * eig->gevec[v*ndim +p2+i];
					contributions[r1*nres + r2] += cblas_ddot (n1, eig->gevec+ v*ndim +p1, 1,  MM         +p1 +ndim*(p2+i), 1) * eig->gevec[v*ndim +p2+i];
				}
			}
		}
	}
	/////////// fill the upper part of the matrix
	for (r1=0; r1< nres; r1++ )for (r2=r1; r2< nres; r2++ ) contributions[r2*nres + r1]= contributions[r1*nres + r2];
	////// write the output
	FILE  *out = fopen(resfile, "w");
	for (r=0; r< nres; r++ ) fprintf(out,"%c-%s%d ", residues[r].chain,residues[r].resname,residues[r].resi );
	fprintf(out,"\n");
	for (r1=0; r1< nres; r1++ )
	{
		for (r2=0; r2< nres; r2++ )
		{
			fprintf(out, "%.6f  ",contributions[r1*nres + r2]);
			//if(r1 == r2 ) fprintf(out, "%.6f  ",contributions[r1*nres + r2]);
			//else fprintf(out, "%.6f  ",2*contributions[r1*nres + r2]);
		}
		fprintf(out,"\n");
	}
	fclose(out);

	/*******************
	double tot= 0;
	for (r1=0; r1< nres; r1++ )
		for (r2=r1; r2< nres; r2++ )
		{
			if(r1 == r2 ) tot += contributions[r1*nres + r2];
			else tot += 2*contributions[r1*nres + r2];
		}
	for (r1=0; r1< eig->rank ; r1++ ) tot -= eig->geigval[r1];
	printf("Difference =%g \n",tot);
	//printf("sub=%g <---> blas=%g \n",tot, vecX_symat_vecY( eig->mat_b,ndim  , eig->gevec , eig->gevec ) + vecX_symat_vecY( MM,ndim  , eig->gevec , eig->gevec ));
	**************************/
	rvec *ref;
	double *local_contr, max;
	int firstAtom;
	if(pdbBfactorsfile)
	{
		printf("Pdb file %s is written with beta factors of the\n residue score of the contributions to the conformational changes \n", pdbBfactorsfile);
		ALLO_ONE_DIM_ARRAY(ref, rvec,natom );
		for (i = 0; i < natom; i++) for (j=0; j< 3 ;j++) ref[i][j]= 	(real) eig->m_b[i*3 +j] ;
		ALLO_ONE_DIM_ARRAY(local_contr, double,natom );
		for (r=0; r< nres; r++ )
		{
			firstAtom=residues[r].index;
			for (i=0; i< residues[r].natoms ; i++ )local_contr[firstAtom +i]= contributions[r*nres + r];
		}
		max= max_array(local_contr, natom);
		for (i = 0; i < natom; i++) local_contr[i] /=max;
		writepdb(pdbBfactorsfile,"Contributions in B-factors", 0,natom,	atoms,index,ref,NULL,local_contr,0);
		free(ref);
		free(local_contr);
	}

	free(MM);
	free(dm);
	free(residues);
	free(contributions);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////  Interpolation of the average structure  m_x along the directions of the vectors of a matrix G (r x c)
////   done by reconstructing x from y around the average of m_x in the y space
///// x_hat= m_x + G(G^tG)^-1 G^t x= G(G^tG)^-1 (y-m_y)= P(y-m_y) + m_x 
///// P = G(G^tG)^-1 
///// (y-m_y) is taken from -3* SD(y) to 3* SD(y)
////// needed in the next function show_dynamic_changes()
void interpolate_coor(char *pdbfname,int r, int c , double *G ,double *m_x, double *ev ,char *title, int natom , t_atoms *atoms,atom_id *index ,double sd)
{
	ONE_DIM_ARRAY(GtG, double, (c*c));
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, c, c, r, 1.0, G, r, G, r, 0.0, GtG, c);
	ONE_DIM_ARRAY(GtG_inv, double, (c*c));
	symmpinv(c, GtG, GtG_inv) ;
	free(GtG);
	ONE_DIM_ARRAY(P, double, (r*c));
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, r, c, c, 1.0, G, r, GtG_inv, c, 0.0, P, r);
	free(GtG_inv);
	//// project the average m_x on the y space m_y= G^t m_x
	//ONE_DIM_ARRAY(m_y, double, c);
	//cblas_dgemv(CblasColMajor , CblasTrans, c, r,1,G,c, m_x, 1, 0,m_y,1);
	/// generate y around m_y -+ 3 SD(y)
	ONE_DIM_ARRAY(x_hat, double, r);
	ONE_DIM_ARRAY(x_fr, rvec, natom);
	ONE_DIM_ARRAY(y, double, c);
	FILE *fout = fopen(pdbfname, "w");
	if (fout  == NULL)
	{
		fprintf(stderr, "Cannot open file %s for writing, exit\n",	pdbfname);
		exit(1);
	}
	int i, model=0;
	/**********************************************************/
	for (i = 0; i < natom; i++)
	{
			x_fr[i][0]=  m_x[i*3] ;
			x_fr[i][1]=  m_x[i*3 +1] ;
			x_fr[i][2]=  m_x[i*3 +2] ;
	}
	writepdblines2file(fout, title, model,natom,atoms,index,x_fr,NULL,NULL);
	fprintf(fout,"ENDMDL\n");
	/*****************************************************/
	double fac= -sd;
	while( fac <= sd)
	{
		model++;
		for (i = 0; i < c; i++) y[i]= fac*ev[i];
		cblas_dgemv(CblasColMajor , CblasNoTrans, r, c,1,P,r, y, 1, 0,x_hat,1);
		for (i = 0; i < natom; i++)
		{
			x_fr[i][0]= x_hat[i*3]+ m_x[i*3] ;
			x_fr[i][1]= x_hat[i*3 +1]+ m_x[i*3 +1] ;
			x_fr[i][2]= x_hat[i*3 +2]+ m_x[i*3 +2] ;
		}
		writepdblines2file(fout, title, model,natom,atoms,index,x_fr,NULL,NULL);
		fac = fac + 0.2;
	}
	fprintf(fout,"END   \n");
	fclose(fout);
	free(P);
	free(x_hat);
	free(x_fr);
	free(y);
}
////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void show_motions_along_vectors(	const char *pdbfile ,    // write the structures to the file
							Geigdecomposition_t *geig,
							int 		natom,
							atom_id 	*index,
							int first_vec,  /// index of first vector
							int last_vec,   /// index of the last vector
							t_atoms *atoms,
							int write_individual)  /// write pdb for the changes along each vector besides the total
{

	int ndim= geig->ndim;
	if(natom != ((int)(ndim/3)))
	{
		PRINT_FATAL_ERR_MES("Nr of atoms in the index is different from ndim/3 in Geigen structure");
	}
	int i,j;
	int nr_gvectors = last_vec - first_vec ;
	int ndx1, ndx2;
	char buffer[200] ;
	char titlebuffer[300] ;
	sprintf(buffer,"%s_b_%dto%d.pdb",pdbfile, first_vec+1,last_vec);
	printf(" motions along vectors %d to %d are saved to the file %s\n",first_vec+1,last_vec+1,buffer);
	sprintf(titlebuffer,"dynamics along vector %d to %d",first_vec+1,last_vec+1);
	interpolate_coor(buffer, ndim, nr_gvectors , geig->gevec + (ndim*first_vec),geig->m_b,geig->geigval + first_vec,titlebuffer, natom , atoms,index,2);
	////// vector-wise
	if(write_individual)
	{
		for (i = first_vec; i < last_vec; i++)
		{
			sprintf(buffer,"%s_motion_%d.pdb",pdbfile, i+1 );
			sprintf(titlebuffer,"increased dynamics along vector%d of lambda %e",i+1,geig->geigval[i] );
			interpolate_coor(buffer, ndim, 1 ,  geig->gevec + (i*ndim),geig->m_b,geig->geigval +i,titlebuffer, natom , atoms,index,2);
		}
	}

}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void show_increased_decresed_motions(	const char *pdbfile ,    // write the projections to a dat file
							Geigdecomposition_t *geig, 
							int 		natom,
							atom_id 	*index,
							double tol,
							int nr_vec_up,  /// number of vectors (high eigenvalues) to be used if tol =0
							int nr_vec_low,  /// number of vectors (high eigenvalues) to be used if tol =0
							t_atoms *atoms,
							int write_individual)  /// write pdb for the changes along each vector besides the total
{
	
	int ndim= geig->ndim;
	if(natom != ((int)(ndim/3)))
	{
		PRINT_FATAL_ERR_MES("Nr of atoms in the index is different from ndim/3 in Geigen structure");
	}
	
	/////////////// sort g-eigenvalues 
	ONE_DIM_ARRAY(gevalInd, int ,(geig->rank) );
	sort_decreasing( geig->geigval ,geig->rank ,gevalInd);
	////////// 
	int i,j, ind_upper, ind_lower ;
	int indiv= write_individual;
	if(tol > 0)
	{
		printf("The sets of eigenvectors are defined via the tolerance value (tol > 0) \n");
		i=0;
		while(geig->geigval[gevalInd[i]] > (1+tol)) i++;
		ind_upper= i;
		while(geig->geigval[gevalInd[i]] > (1-tol)) i++;
		ind_lower= i;
		indiv=0;
	}
	else
	{	
		ind_upper= nr_vec_up;
		ind_lower= geig->rank - nr_vec_low;
	}
	////////////
	int nr_gvectors , v;
	int ndx1, ndx2;
	char buffer[200] ;
	char titlebuffer[700] ;
	////////////////////// increased motions around average state b////////////////////////
	printf(" increased motions saved to the file %s\n", pdbfile);
	nr_gvectors=ind_upper;
	ONE_DIM_ARRAY(G, double ,(nr_gvectors*ndim) );
	ONE_DIM_ARRAY(ev, double ,(nr_gvectors) );	
	for (i = 0; i < nr_gvectors; i++) 
	{
		ev[i]=geig->geigval[gevalInd[i]];
		ndx1=ndim *gevalInd[i];
		ndx2=ndim *i;
		for (j = 0; j< ndim; j++) G[ndx2+j]=geig->gevec[ndx1 +j];
	}
	sprintf(buffer,"%s_increased_b_%dto%d.pdb",pdbfile, 1,nr_gvectors );
	sprintf(titlebuffer,"increased dynamics along vector %d to %d gevalues from %e to %e",1,i+1,ev[0] ,ev[nr_gvectors-1] );
	interpolate_coor(buffer, ndim, nr_gvectors , G,geig->m_b, ev,titlebuffer, natom , atoms,index,2);
	////// vector-wise
	if(indiv)
	{
		for (i = 0; i < nr_gvectors; i++)
		{
			sprintf(buffer,"%s_increased_b_%d.pdb",pdbfile, i+1 );
			sprintf(titlebuffer,"increased dynamics along vector%d of lambda %e",i+1, ev[i] );
			interpolate_coor(buffer, ndim, 1 , G + (i*ndim),geig->m_b, ev+i,titlebuffer, natom , atoms,index,2);
		}
	}
	free(G);
	free(ev);
	//////////////////////////  decreased motions  /////////////////////////////////
	printf(" decreased motions saved to the file %s\n", pdbfile);
	nr_gvectors= geig->rank  - ind_lower ;
	ALLO_ONE_DIM_ARRAY(G, double ,(nr_gvectors*ndim) );
	ALLO_ONE_DIM_ARRAY(ev, double ,(nr_gvectors) );
	for (i = 0; i < nr_gvectors; i++)
	{
		v=gevalInd[i+ind_lower];
		ev[i]=1;
		ndx1=ndim *v;
		ndx2=ndim *i;
		for (j = 0; j< ndim; j++) G[ndx2+j]=geig->gevec[ndx1 +j];
	}
	sprintf(buffer,"%s_decreased_a_-%dto-%d.pdb",pdbfile,nr_gvectors,1);
	sprintf(titlebuffer,"decreased dynamics along vector %d to %d gevalues from %e to %e",1,i+1,geig->geigval[gevalInd[geig->rank-1]] ,geig->geigval[gevalInd[ind_lower]]);
	interpolate_coor(buffer, ndim, nr_gvectors , G,geig->m_a, ev,titlebuffer, natom , atoms,index,8);
	if(indiv)
	{
		for (i = 0; i < nr_gvectors; i++)
		{
			sprintf(buffer,"%s_decreased_a_-%d.pdb",pdbfile, nr_gvectors -i);
			sprintf(titlebuffer,"decreased dynamics along vector%-d of lambda %e", nr_gvectors -i, geig->geigval[gevalInd[i+ind_lower]] );
			interpolate_coor(buffer, ndim, 1 , G + (i*ndim),geig->m_b, ev+i,titlebuffer, natom , atoms,index,8);
		}
	}
	free(G);
	free(ev);
	/**************
	//////////////////////////  unchanged motions  /////////////////////////////////
	if(tol > 0)
	{
		nr_gvectors=ind_lower - ind_upper;
		printf(" unchanged motions along %d vectors are saved to the file %s\n", nr_gvectors,pdbfile);
		ALLO_ONE_DIM_ARRAY(G, double ,(nr_gvectors*ndim) );
		ALLO_ONE_DIM_ARRAY(ev, double ,(nr_gvectors) );
		for (i = 0; i < nr_gvectors; i++)
		{
			v=gevalInd[i+ind_upper];
			ev[i]=geig->geigval[v];
			ndx1=ndim *v;
			ndx2=ndim *i;
			for (j = 0; j< ndim; j++) G[ndx2+j]=geig->gevec[ndx1 +j];
		}
		iterpolate_coor("dynamic_unchanged_a.pdb", ndim, nr_gvectors , G,geig->m_a, ev,"unchanged_dynamics_a", natom , atoms,index,10);
		free(G);
		free(ev);
	}
	***************/
		
	free(gevalInd);	
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


