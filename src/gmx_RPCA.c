#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "typedefs.h"
#include "vec.h"
#include "index.h"
#include "rmpbc.h"
#include "fitting.h"
#include "covariance.h" 
#include "arrays.h" 
#include "CWfitting.h"
#include "linear_algebra.h"
#include "strTools.h"
#include "xdrfile_xtc.h"


////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    	const char     *desc[] = 
	{
        "A tool for Relative Principal Components Analysis"
    };
    static gmx_bool   bPBC = TRUE  , cw_fit=TRUE ;
	static int   verbose=1, fit_b=1,  algo=0;
    static int     first=1,last=-1,  TC=10, skip=0 ;
	static real timeBeg=0, timeEnd=-1, convergence = 0.00001 ;
	t_pargs    pa[] = 
	{
		{ "-fit_b",  TRUE, etINT, {&fit_b}, "0: use the structure file of B for the average. 1:Use GPA to get the average of state B 2: useCW fitting to the average of the state A."},
		{ "-bt",  TRUE, etREAL, {&timeBeg},  "start reading the frame at this time" },
		{ "-et",  TRUE, etREAL, {&timeEnd},  "end reading the frame at this time; -1 mean no upper limit" },
		{ "-algo",  TRUE, etINT, {&algo}, "Algorithm used for relative pca: 1, sub spacing, 0 non sub spacing"},
		{ "-convergence",  TRUE, etREAL, {&convergence},  "The convergence threshold of the total rmsd of all the conformations" },
		{ "-cwfit",  TRUE, etINT, {&cw_fit}, "use the cw fitting between the averages"},
        { "-first", FALSE, etINT, {&first}, "First eigenvector for analysis (-1 is select)" },
		{ "-last",  FALSE, etINT, {&last}, "Last eigenvector for analysis (-1 is till the last)" },
		{ "-skip",  FALSE, etINT, {&skip}, "skip this number of frames between every selected frame for projection" },
		{ "-pbc",  FALSE,  etBOOL, {&bPBC}, "Apply corrections for periodic boundary conditions" },
		{ "-verbose",  TRUE,  etINT, {&verbose}, "print additional details" }
    };
    int ndim, i, j, d;
    const char     *fitfile_a,*fitfile_b,*trxfile_a ,*trxfile_b,*ndxfile_a ,*ndxfile_b ;
    const char     *averfile_a, *averfile_b;
    const char     *covfile_a, *covfile_b, *eigvalfile , *eigvecfile  ,*dklfile,  *projfile_a , *projfile_b ,*resfile, *respdb;
    t_filenm  fnm[] = 
	{
		///// input files
        { efTRX, "-fa",  NULL, ffREAD },
		{ efTRX, "-fb",  NULL, ffREAD }, 
        { efTPS, "-sa",  "structure_a.tpr", ffREAD },
		{ efTPS, "-sb",  "structure_b.tpr", ffREAD },
		{ efMTX, "-covar_a",  "covar_a.mtx", ffREAD },
	    { efNDX, "-na",  "index_a.ndx", ffOPTRD },
		{ efNDX, "-nb",  "index_b.ndx", ffOPTRD },
		///// output files
		{ efMTX, "-covar_b",  "covar_b.mtx", ffOPTWR },
        { efSTO, "-av_a", "average_a.pdb", ffOPTWR },
		{ efSTO, "-av_b", "average_b.pdb", ffOPTWR },
 		{ efXVG, "-eval",  "Geigenval.xvg", ffOPTWR },
		{ efXVG, "-DKL",  "d_kl.xvg", ffOPTWR },
		{ efTRN, "-geig",  "geigen.geig", ffOPTWR },
		{ efDAT, "-proj_a", "projection_a.dat", ffOPTWR},
		{ efDAT, "-proj_b", "projection_b.dat", ffOPTWR},
		{ efDAT, "-res", "residue.dat", ffOPTWR},
		{ efSTO, "-respdb", "contributions.pdb", ffOPTWR },
    };
    #define NFILE asize(fnm)
	output_env_t    oenv;
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE, NFILE,
					   fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	fitfile_a    = opt2fn("-sa", NFILE, fnm);
	fitfile_b    = opt2fn("-sb", NFILE, fnm);
	trxfile_a    = opt2fn("-fa", NFILE, fnm);
	trxfile_b    = opt2fn("-fb", NFILE, fnm);
	covfile_a = opt2fn("-covar_a", NFILE, fnm);
	ndxfile_a    = opt2fn_null("-na", NFILE, fnm);
	ndxfile_b    = opt2fn_null("-nb", NFILE, fnm);
	///////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// control the output files //////////////////////////////
	covfile_b = opt2fn("-covar_b", NFILE, fnm);
	averfile_a = opt2fn("-av_a", NFILE, fnm);
	averfile_b = opt2fn("-av_b", NFILE, fnm);
	eigvalfile = opt2fn("-eval", NFILE, fnm);
	eigvecfile = opt2fn("-geig", NFILE, fnm);
	dklfile    = opt2fn("-DKL", NFILE, fnm);
	resfile    = opt2fn_null("-res", NFILE, fnm);
	respdb    = opt2fn_null("-respdb", NFILE, fnm);
	projfile_a   = opt2fn_null("-proj_a", NFILE, fnm);
	projfile_b   = opt2fn_null("-proj_b", NFILE, fnm);

	///////////////////////////////////////////////////////////////////////////////////////
	////////////////////                read structure files                ///////////////
	///////////////////////////////////////////////////////////////////////////////////////
    t_topology      top_a, top_b;
    int             ePBC_a, ePBC_b;
    t_atoms        *atoms_a, *atoms_b;
    rvec           *xread,   *xrefT_a, *xrefT_b;
	matrix          box_a, box_b;
	gmx_rmpbc_t     gpbc_a = NULL,  gpbc_b =NULL;
	char            str[STRLEN] ;
	read_tps_conf(fitfile_a, str, &top_a, &ePBC_a, &xrefT_a, NULL, box_a, TRUE);
    atoms_a = &top_a.atoms;
	read_tps_conf(fitfile_b, str, &top_b, &ePBC_b, &xrefT_b, NULL, box_b, TRUE);
    atoms_b = &top_b.atoms;
    if (bPBC)
    {
        gpbc_a = gmx_rmpbc_init(&top_a.idef, ePBC_a, atoms_a->nr, box_a);
        gmx_rmpbc(gpbc_a, atoms_a->nr, box_a, xrefT_a);
		gpbc_b = gmx_rmpbc_init(&top_b.idef, ePBC_b, atoms_b->nr, box_b);
        gmx_rmpbc(gpbc_b, atoms_b->nr, box_b, xrefT_b);
    }
    //ONE_DIM_ARRAY(xref_b, rvec ,natoms_b);
    //for (i = 0; i < natoms_b; i++) for (d = 0; d < DIM; d++) xref_b[i][d] =xrefT_b[i][d];
	///////////////////////////////////////////////////////////////////////////////////////
	////////////////////       read index files and do the selections       ///////////////
	///////////////////////////////////////////////////////////////////////////////////////
	int     nfit_a, nfit_b ;
    atom_id   *ifit_a, *ifit_b;
	char  *fitname_a,*fitname_b ;
	printf("\n Choose a group for the covariance analysis( and least squares fit) of state A: \n");
   	get_index(atoms_a, ndxfile_a, 1, &nfit_a, &ifit_a, &fitname_a);
	printf("\n Choose a group for the covariance analysis( and least squares fit) of state B \n");
	printf("\n This goup should be the same which was selected for state A: \n");
	get_index(atoms_b, ndxfile_b, 1, &nfit_b, &ifit_b, &fitname_b);
	///////////////////////////////////////////////////////////////////////////////////
	////////// check the correspondence between the selected groups in A & B //////////
	printf("#############################################\n");
    printf("Checking the correspondence between the selected \nresidues and atoms of the two states... \n");
	if (nfit_a != nfit_b)
	{
		fprintf(stderr,"Error:Numbers of selected atoms in A and B are different!\n");
		exit(EXIT_FAILURE);	
	} 
    if (nfit_a < 3)
	{
		fprintf(stderr,"Error:Need >= 3 points to fit!\n");
		exit(EXIT_FAILURE);
	}
	int a,b ,r1, r2;
	for(i=0; i < nfit_a; i++) 
	{
		a= ifit_a[i];
		r1= atoms_a->atom[ a ].resind;
		b= ifit_b[i];
		r2= atoms_b->atom[ b ].resind;
		if ( (strcmp(*(atoms_a->atomname[a]) ,*(atoms_b->atomname[b]) ) != 0) ||
			 (strcmp( *(atoms_a->resinfo[r1].name),*(atoms_b->resinfo[r2].name)) != 0)	)
		{
			fprintf(stderr,"Error: Residue:atom %s%d:%s%d form state A is not corresponding to %s%d:%s%d from state B!\n",
					*(atoms_a->resinfo[r1].name),r1+1,*(atoms_a->atomname[a]),a+1, *(atoms_b->resinfo[r2].name),r2+1, *(atoms_b->atomname[b]),b+1);
			exit(EXIT_FAILURE);
		}
	}

	printf("Correspondence is OK!\n" );
	printf("#############################################\n");
	////////////////////////////////////////////////////////////////////////
	/////////// read the covariance matrix and mean of state A /////////////
	ndim = nfit_a*DIM;
	int n_tmp;
	double *mat_a, *m_a ;
	printf("#############################################\n");
	printf("Reading the mean & the covariance matrix\n of state A from the file:\n %s\n",covfile_a );
	if (  ! read_cov_matrix(covfile_a ,&mat_a, &m_a, &n_tmp)  ) 
	{
		fprintf(stderr,"Error in reading the mean and the covariance matrix of state A \n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
	if(n_tmp != ndim)
	{
		fprintf(stderr,"Error: Nr of dims in the covariance matrix file is different from the Nr. of the selected group of the structure file\n");
		exit(EXIT_FAILURE);
	}
	ONE_DIM_ARRAY(ref_a, rvec ,nfit_a);
	for (i = 0; i < nfit_a; i++) for (j=0; j< 3 ;j++) ref_a[i][j]= 	(real) m_a[i*3 +j] ;
	////////////// generalized inverse is needed for CW fitting
	printf("#############################################\n");
	printf("Computing the generalized inverse of the  covariance \n matrix A which is needed for CW fitting \n");
	ONE_DIM_ARRAY(ginv,double,(ndim*ndim)) ;
	symmpinv(ndim, mat_a, ginv);
	////////////////////////////////////////////////////////////////////////
	/////////// read the covariance matrix and mean of state B /////////////
	ONE_DIM_ARRAY(mat_b,double,(ndim*ndim)) ;
	ONE_DIM_ARRAY(m_b,double,ndim) ;
	int natoms_b;
	read_xtc_natoms(trxfile_b,&natoms_b);
	rvec *xav_b;



	if (fit_b==0)
	{
		///////////////////       average is known    //////////////////////
		printf("#############################################\n");
		printf("The average of state B  is obtained from the structure file %s\n", fitfile_b);
		printf("MAKE SURE that the atoms in the index corresponding to the traj)\n");
		ALLO_ONE_DIM_ARRAY(xav_b, rvec ,natoms_b);
		for (i = 0; i < natoms_b; i++) for (d = 0; d < 3; d++) xav_b[i][d] =xrefT_b[i][d];
	}
	if(fit_b==1)
	{
		///////////////////       GPA  of B        //////////////////////
		printf("#############################################\n");
		printf("********** GPA fitting of state B ***********\n");
		printf("#############################################\n");
		printf("obtain the average of state B via GPA fitting\nfor the conformations in the file:\n%s\n", trxfile_b);
		gpa(trxfile_b,   &xav_b,  nfit_b,ifit_b, bPBC,gpbc_b, TC, convergence, verbose, timeBeg, timeEnd);
	}



	if(fit_b ==0 || fit_b ==1)
	{
		/////////////////// between states fitting //////////////////////
		printf("#############################################\n");
		printf("*********************************************\n");
		if(verbose) printf("RMSD (B2A) before= %g\n",calc_rmsd(  ref_a, NULL,xav_b,ifit_b, nfit_b));
		center_x2ref_diffrenct_atoms(xav_b, ifit_b,ref_a ,NULL, nfit_b, natoms_b );
		rotate_B2R(                  xav_b, ifit_b,ref_a, NULL ,nfit_b, natoms_b );
		if(verbose) printf("RMSD after= %g\n",calc_rmsd(ref_a, NULL,xav_b,ifit_b, nfit_b));
		if( cw_fit)
		{
			printf("Covariance weighted fitting of the average structure\n of B to the average structure of A\n" );
			CW_fitting(xav_b,ifit_b, ref_a, NULL ,nfit_b, natoms_b,ginv,verbose);
		}
		printf("#############################################\n");
		printf("Computing the covariance matrix of state B ...\n" );
		blas_compute_covariance(trxfile_b, xav_b, mat_b, nfit_b, ifit_b,nfit_b, ifit_b,bPBC, gpbc_b,  m_b,  timeBeg, timeEnd) ;
	}
	if( fit_b==2)
	{
		Covariance_CW_fitting_traj(trxfile_b ,NULL, ref_a, ginv,nfit_b,ifit_b, nfit_b,ifit_b, bPBC,  gpbc_b, atoms_b,mat_b,m_b, NULL, verbose,  timeBeg, timeEnd);
	}
	///// make a reference for next fitting of the selected atoms only
	ONE_DIM_ARRAY(ref_b, rvec ,nfit_b);
	for (i = 0; i < nfit_b; i++) for (j=0; j< 3 ;j++) ref_b[i][j]= 	(real) m_b[i*3 +j] ;
	//Covariance_CW_fitting_traj(trxfile_a ,NULL, ref_a, ginv,nfit_b,ifit_b, nfit_b,ifit_b, bPBC,  gpbc_b, oenv,atoms_b,mat_b,m_b, NULL, verbose);
	///////////////////////////////////////////////////////////////////
	if(covfile_b )
	{
		printf("Writing the average and the covariance matrix \nof state B to the file %s \n", covfile_b);
		if ( !write_cov_matrix(covfile_b ,&mat_b ,&m_b ,&ndim,  verbose)  )
			fprintf(stderr,"Error in writing the covariance matrix\n");
	}
	///////////////////////////////////////////////////////////////////
	//////////////////// save averages structures ///////////////////// 
	printf("#############################################\n");
	printf("saving the average structures of states A& B\nto the files %s & %s\n",averfile_a, averfile_b);
	writepdb(averfile_a,"Average structure of state A", 0,nfit_b,	atoms_b,ifit_b,ref_a,NULL,NULL,0);
	writepdb(averfile_b,"Average structure of state B", 0, nfit_b,	atoms_b,ifit_b,ref_b,NULL,NULL,0);
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
	//else  relative_pca_subspacing(eig,  1,verbose );
	write_array2xvg(eigvalfile ,eig->rank,eig->geigval,"Eigenvalues of the covariance matrix" ,"Eigenvector index","(nm\\S2\\N)" );
	write_Geigen2xvg(dklfile, eig )	;
	printf("*********************************************\n");
	printf("Writting the Generalized eigenpairs and the means\nto the file %s\n", eigvecfile);
	if( !write_Geigen(eigvecfile,eig , verbose)) fprintf(stderr," Can not Write G eigenpairs to the file %s \n",eigvecfile );
	printf("#############################################\n");
	//ONE_DIM_ARRAY(eig2, Geigdecomposition_t,1 );
	//if (!read_Geigen(eigvecfile,eig2 ) )  fprintf(stderr," Can not read G eigenpairs from the file %s \n",eigvecfile );
	//////////////////////////////////////////////////////////////////////
	/////////////////////// interactions map   ///////////////////////////
	int Beg_evec, End_evec;
	if( first < 0) Beg_evec= 0;
	else Beg_evec= first -1;
	if(last > eig->rank || last == -1 ) End_evec= eig->rank;
	else End_evec= last;
	printf("#############################################\n");
	if(projfile_a && projfile_b)
	{
		RPCA_project(trxfile_a,TRUE, projfile_a ,eig,  Beg_evec,End_evec,nfit_b,ifit_a,bPBC,gpbc_b,  timeBeg, timeEnd, skip);
		RPCA_project(trxfile_b,FALSE, projfile_b ,eig, Beg_evec,End_evec,nfit_b,ifit_b,bPBC,gpbc_b, timeBeg, timeEnd, skip);
	}
	printf("#############################################\n");
	if(resfile)
	{
		printf("Computing the interaction map\n");
		interaction_map(resfile, eig,atoms_b,ifit_b,Beg_evec, End_evec,respdb);
	}
return 0;

}


