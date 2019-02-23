#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "covariance.h" 
#include "arrays.h" 
#include "linear_algebra.h"
#include "xdrfile_xtc.h"
#include "statutil.h"

////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    	const char     *desc[] = 
	{
        "Perform simultaneous diagonalization of the covariance matrices. Optionally, the subspacing optimal algorithm can be used."
    };
	static int   verbose=1,  algo=0;
   	t_pargs    pa[] =
	{
			{ "-algo",  TRUE, etINT, {&algo}, "Algorithm used for relative pca: 1, sub spacing, 0 non sub spacing"},
			{ "-verbose",  TRUE,  etINT, {&verbose}, "print additional details" }
    };

    const char     *covfile_a, *covfile_b, *eigvalfile , *eigvecfile  ,*dklfile;
    t_filenm  fnm[] = 
	{
		///// input files
        { efMTX, "-covar_a",  "covar_a.mtx", ffREAD },
		{ efMTX, "-covar_b",  "covar_b.mtx", ffREAD },
	    //// output files
 		{ efXVG, "-eval",  "Geigenval.xvg", ffOPTWR },
		{ efXVG, "-DKL",  "d_kl.xvg", ffOPTWR },
		{ efTRN, "-geig",  "geigen.geig", ffOPTWR }
    };
	#define asize(a) (sizeof(a)/sizeof((a)[0]))
    #define NFILE asize(fnm)

	output_env_t    oenv;
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE, NFILE,
					   fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	covfile_a = opt2fn("-covar_a", NFILE, fnm);
	covfile_b = opt2fn("-covar_b", NFILE, fnm);
	///////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// control the output files //////////////////////////////
	eigvalfile = opt2fn("-eval", NFILE, fnm);
	eigvecfile = opt2fn("-geig", NFILE, fnm);
	dklfile    = opt2fn("-DKL", NFILE, fnm);
	///////////////////////////////////////////////////////////////////////////////////////
	///////////////// read the covariance matrix and mean of state A   ////////////////////
	int ndim;
	double *mat_a, *m_a ;
	printf("#############################################\n");
	printf("Reading the mean & the covariance matrix\n of state A from the file:\n %s\n",covfile_a );
	if (  ! read_cov_matrix(covfile_a ,&mat_a, &m_a, &ndim)  )
	{
		fprintf(stderr,"Error in reading the mean and the covariance matrix of state A \n");
		exit(EXIT_FAILURE);
	}
	printf("#############################################\n");
	int ndim_b;
	double *mat_b, *m_b ;
	printf("Reading the mean & the covariance matrix\n of state B from the file:\n %s\n",covfile_b );
	if (  ! read_cov_matrix(covfile_b ,&mat_b, &m_b, &ndim_b)  )
	{
		fprintf(stderr,"Error in reading the mean and the covariance matrix of state B \n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
	if(ndim_b != ndim)
	{
		fprintf(stderr,"Error: Nr of dims in the covariance matrix file of state A \nis different from the number in the file of state B\n");
		exit(EXIT_FAILURE);
	}
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////      RPCA  analysis     ///////////////////////////////////
	printf("#############################################\n");
	printf("******* Simultaneously diagonalizing ********\n");
	printf("*********************************************\n");
	ONE_DIM_ARRAY(eig, Geigdecomposition_t,1 );
	eig->mat_a = mat_a;
	eig->mat_b = mat_b;
	eig->ndim = ndim;
	eig->m_a =m_a;
	eig->m_b =m_b;
	relative_pca(eig,  1,algo,verbose );

	write_array2xvg(eigvalfile ,eig->rank,eig->geigval,"Eigenvalues of the covariance matrix" ,"Eigenvector index","(nm\\S2\\N)" );
	write_Geigen2xvg(dklfile, eig )	;
	printf("*********************************************\n");
	printf("Writting the Generalized eigenpairs and the means\nto the file %s\n", eigvecfile);
	if( !write_Geigen(eigvecfile,eig , 2)) fprintf(stderr," Can not Write G eigenpairs to the file %s \n",eigvecfile );
	printf("#############################################\n");
	//////////////////////////////////////////////////////////////////////
return 0;
}


