 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "statutil.h"
#include "tpxio.h"
#include "pbc.h"
#include "index.h"
#include "rmpbc.h"
#include "arrays.h"
#include "fitting.h"
#include "covariance.h" 
#include "linear_algebra.h"
#include "xdrfile_xtc.h"
#include "strTools.h"

#include "CWfitting.h"

//
//#include "gromacs_tools.h"
//#include "arg_parser.h"
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    	const char     *desc[] = 
	{
        	"\n\nThis tool can do the following tasks:\n \
			 1. perform GPA fitting for the conformations in the trajactory in the file (-f). The output is the average structure (-av). \n \
			 2. fit the conformations to the GPA average. The output trajectory  (-o) \n    \
			 3. compute the covariance matrix and write it to the binary file ( covar covar.mtx) \n    \
			 4. perform eigen-decomposition of the covariance matrix and write the mean and eigenpairs to the file (-eigen eigen.mvv) \n\n "
        	
    };

	static gmx_bool outFit = FALSE, bPBC = TRUE , pca=FALSE ,covar=FALSE ;  
    static int   nogpa=0,  TC=10, verbose=0;
	static real  timeBeg=0, timeEnd=-1, convergence = 0.00001 ;
	t_pargs    pa[] = 
	{
		{ "-nogpa",  FALSE,  etBOOL, {&nogpa}, "fit structures to the structure file instesd of the default GPA fitting (to the gpa average)" },
		{ "-bt",  TRUE, etREAL, {&timeBeg},  "start reading the frame at this time" },
		{ "-et",  TRUE, etREAL, {&timeEnd},  "end reading the frame at this time; -1 mean no upper limit" },
		{ "-Max_iterations",  TRUE, etINT, {&TC},  "The max nr. GPA of iteration" }, 
		{ "-convergence",  TRUE, etREAL, {&convergence},  "The convergence threshold of the total rmsd of all the conformations" },
	    { "-pbc",  FALSE,  etBOOL, {&bPBC},    	"Apply corrections for periodic boundary conditions" },
		{ "-verbose",  TRUE, etINT, {&verbose},  "Perform addational check....." }, 
		
    };

    const char     *fitfile, *trxfile, *ndxfile, *covarfile , *eigenfile ;
    const char     *out_file,  *averfile, *eigvalfile;
    int             i, d;

    t_filenm  fnm[] = 
	{
        { efTRX, "-f",  NULL, ffREAD },
		{ efTRO, "-o",  "GPA_fitted.xtc", ffOPTWR },
		{ efTRO, "-cw",  "CW_fitted.xtc", ffOPTWR },
        { efTPS, NULL,  NULL, ffREAD },
        { efNDX, NULL,  NULL, ffOPTRD },
        { efSTO, "-av", "average.pdb", ffOPTWR },
		{ efMTX, "-covar",  "covar.mtx", ffOPTWR },
		{ efTRN, "-eigen",  "eigen.mvv", ffOPTWR },
		{ efXVG, "-eval",  "eigenval", ffOPTWR } 
	};
	#define NFILE asize(fnm)
    /////////////////////////////////////////////////////////////////////////
    output_env_t    oenv;
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
						NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	fitfile    = ftp2fn(efTPS, NFILE, fnm);
    trxfile    = ftp2fn(efTRX, NFILE, fnm);
	averfile   = ftp2fn(efSTO, NFILE, fnm);
	ndxfile    = ftp2fn_null(efNDX, NFILE, fnm);
	out_file   = ftp2fn_null(efTRO, NFILE, fnm);
	covarfile = opt2fn_null("-covar", NFILE, fnm);
	eigenfile = opt2fn_null("-eigen", NFILE, fnm);
	if(eigenfile) eigvalfile = opt2fn("-eval", NFILE, fnm);
	/////////////////////////////////////////////////////////////////////////
	////////////////////          selections                  ///////////////
	char str[STRLEN];
	t_topology      top;
	int             ePBC;
	t_atoms        *atoms;
	rvec            *xrefT;
	matrix          box;
    read_tps_conf(fitfile, str, &top, &ePBC, &xrefT, NULL, box, TRUE);
    atoms = &top.atoms;
    atom_id        *cov_index, *ifit, *outindex;
    int             cov_natoms, nfit,outnatoms;
    char            *cov_name, *fitname, *outname;
	printf("\n Choose a group for the least squares fit\n\n");
	get_index(atoms, ndxfile, 1, &nfit, &ifit, &fitname);
    if (nfit < 3)   fprintf(stderr, "Need >= 3 points to fit!\n");
    printf("\n Choose a group for the outputs: %s, and fitted traj %s(optional):\n\n", averfile, out_file);
    get_index(atoms, ndxfile, 1, &outnatoms, &outindex, &outname);
    if(covarfile || eigenfile)
    {
    	printf("\nChoose a group for the covariance analysis\n");
    	printf("IT IS HIGHLY  RECOMMENDED TO BE THE SAME FITTING GROUP\n\n");
    	get_index(atoms, ndxfile, 1, &cov_natoms, &cov_index, &cov_name);

    }

	//////////// Prepare reference frame ////////////
    gmx_rmpbc_t     gpbc = NULL;
    if (bPBC)
    {
       	gpbc = gmx_rmpbc_init(&top.idef, ePBC, atoms->nr, box);
        gmx_rmpbc(gpbc, atoms->nr, box, xrefT);
    }
    center_x(xrefT, ifit,  nfit,atoms->nr, NULL );
    printf( "*****************************\n Calculating the average structure\n");
	int nat;
	read_xtc_natoms(trxfile,&nat);
	if (nat != atoms->nr) fprintf(stderr,"\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n", atoms->nr, nat);
	rvec *xav, *xref;
	if (!nogpa)
	{
		//////////////////////  MAIN LOOP Of  GPA   /////////////////////////////
		/////////////////////////////////////////////////////////////////////////
		printf("\n\n*******************************************************************\n");
		printf(    "***********          Performing GPA fitting        ****************\n");
		printf(  "\n*******************************************************************\n");
		printf(    "* The conformations are read from the file %s *\n", trxfile);
		gpa(trxfile, &xav,  nfit,ifit, bPBC,gpbc, TC, convergence, verbose, timeBeg, timeEnd);
		/////  write the average structure ////////
		printf("\n The average structure of the selected group (%s) \n will be written to file: %s\n",outname, averfile );
		if (outnatoms > nat)
		{
			fprintf(stderr,"Error: Number of selected atoms for writing (%d) is larger than number in the trajectory (%d)\n",outnatoms,nat);
			fprintf(stderr,"Error: average %s will not be written\n", averfile);

		}
		else   writepdb(averfile,"Average structure", 0,outnatoms,	atoms,outindex,xav,outindex,NULL,0);
		printf("GPA is done\n");
		////////////////////////////////////////////////////////////////////////////////////////						
		/////////////////////////////////// output the fitted traj /////////////////////////////
		if (out_file)
		{
			printf("\n\n*******************************************************************\n");
			printf("\n Fitting the conformations to the (GPA) average \n and Writing trajectory to the file %s\n", out_file );
			if (outnatoms > nat)
			{
				fprintf(stderr,"Error: Number of selected atoms for writing (%d) is larger than number in the trajectory (%d)\n",outnatoms,nat);
				fprintf(stderr,"Error: fitted trajectory  %s will not be written\n", out_file);
			}
			else write_fitted_traj(trxfile, out_file, xav, nfit,ifit, outnatoms, outindex, bPBC, gpbc, timeBeg, timeEnd);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////     Constructing covariance matrix   ///////////////////////
	int ndim = cov_natoms*DIM;
	double *mat, *m_a;  
	if(covarfile || eigenfile)
	{
		ALLO_ONE_DIM_ARRAY(mat, double,  (ndim*ndim));	
		ALLO_ONE_DIM_ARRAY(m_a, double,ndim );
		printf("*********************************************************************\n");
		printf("***************    Computing covariance matrix     ******************\n");
		printf("Group %s will be used for the covariance analysis\n",cov_name);
		if (nogpa)
		{
			printf("The structure in the file %s will be used as a reference for fitting the structures(nogpa=on) \n", fitfile);
			ALLO_ONE_DIM_ARRAY(xref, rvec ,nat);
			for (i = 0; i < nat; i++) for(d = 0; d < DIM; d++) xref[i][d] =xrefT[i][d];
		    center_x(xref, ifit,  nfit,nat, NULL );
		}
		else /// the default
		{
			center_x(xav, ifit,  nfit,nat, NULL );
			xref= xav;
			printf("The GPA average will be used as a reference for fitting the structures \n");
		}
				
		blas_compute_covariance(trxfile, xref, mat, nfit, ifit,cov_natoms, cov_index,bPBC, gpbc,  m_a,  timeBeg, timeEnd) ;
	}
	if(covarfile )
	{
		printf("Writing the covariance matrix to the file %s \n", covarfile);
		if ( !write_cov_matrix(covarfile ,&mat ,&m_a ,&ndim,  verbose)  ) 
			        fprintf(stderr,"Error in writing the covariance matrix\n");
	}
	printf("*********************************************************************\n");
	//////////////////////////////////////////////////////////////////////////////////
	/////////////////            eigen decomposition            //////////////////////
	eigdecomposition_t *eig;
	if (eigenfile)
	{
		printf("************       Performing eigen decomposition       ************* \n");
		ALLO_ONE_DIM_ARRAY(eig, eigdecomposition_t,1 );
		eigen(mat, ndim,eig);
		printf("Writing the eigenvalues to the file: %s\n",eigvalfile);
		write_array2xvg(eigvalfile ,ndim,eig->eval,"Eigenvalues of the covariance matrix" ,"Eigenvector index","(nm\\S2\\N)" );
    	///// write the eigenvectors
		printf("Writing the average ,eigenvalues and eigenvectors to the binary file: %s\n",eigenfile);
		write_eigen(eigenfile,eig , &m_a , 0);
		printf("*********************************************************************\n");
		free(m_a);
		free(eig);
	}
	if(covarfile || eigenfile) free(mat);
	free(xav);
	if (nogpa) free(xref);
	return 0;
}




