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
        "Relative Principal Components Analysis:"
    };
    static gmx_bool   bPBC = TRUE , sep=FALSE;
    static int     first=1,last=-1,nv=1  ,verbose=1,  skip=0 ;
	static real timeBeg=0, timeEnd=-1 , tol=0;
	t_pargs    pa[] = 
	{

		{ "-bt",  TRUE, etREAL, {&timeBeg},  "start reading the frame at this time" },
		{ "-tol",  TRUE, etREAL, {&tol},  "Tolerance to consider generalized eigenvalue =1" },
		{ "-et",  TRUE, etREAL, {&timeEnd},  "end reading the frame at this time; -1 mean no upper limit" },
		{ "-nvec", FALSE, etINT, {&nv}, "number of vectors to interpolate the latent variable when showing the dynamics" },
		{ "-sep",  FALSE,  etBOOL, {&sep}, "Write pdb file for the motion along every vector from nvec" },
        { "-first", FALSE, etINT, {&first}, "First eigenvector for analysis (-1 is select)" },
		{ "-last",  FALSE, etINT, {&last}, "Last eigenvector for analysis (-1 is till the last)" },
		{ "-skip",  FALSE, etINT, {&skip}, "skip this number of frames between every selected frame for projection" },
		{ "-pbc",  FALSE,  etBOOL, {&bPBC}, "Apply corrections for periodic boundary conditions" },
		{ "-verbose",  TRUE,  etINT, {&verbose}, "print additional details" }
    };
    int ndim, i, j, d;
    const char     *fitfile_a,*fitfile_b,*trxfile_a ,*trxfile_b,*ndxfile_a ,*ndxfile_b ;
    const char     *averfile_a, *averfile_b;
    const char      *covfile_b, *eigvalfile , *eigvecfile  ,*dklfile,  *projfile_a , *projfile_b ,*resfile, *respdb, *dynamicspdb;
    t_filenm  fnm[] = 
	{
		///// input files
		{ efTRN, "-geig",  "geigen.geig", ffREAD },
		{ efMTX, "-covar_b",  "covar_b.mtx", ffOPTRD },
        { efTRX, "-fa",  NULL, ffREAD },
		{ efTRX, "-fb",  NULL, ffREAD }, 
        { efTPS, "-sa",  "structure_a.tpr", ffREAD },
		{ efTPS, "-sb",  "structure_b.tpr", ffREAD },
	    { efNDX, "-na",  "index_a.ndx", ffOPTRD },
		{ efNDX, "-nb",  "index_b.ndx", ffOPTRD },
		///// output files
		{ efDAT, "-proj_a", "projection_a.dat", ffOPTWR},
		{ efDAT, "-proj_b", "projection_b.dat", ffOPTWR},
		{ efDAT, "-res", "residue.dat", ffOPTWR},
		{ efSTO, "-respdb", "contributions.pdb", ffOPTWR },
		{ efSTO, "-dyn", "dynamics.pdb", ffOPTWR },
		
    };
	#define asize(a) (sizeof(a)/sizeof((a)[0]))
    #define NFILE asize(fnm)
	output_env_t    oenv;
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE, NFILE,
					   fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	eigvecfile = opt2fn("-geig", NFILE, fnm);
	covfile_b = opt2fn_null("-covar_b", NFILE, fnm);
    fitfile_a    = opt2fn("-sa", NFILE, fnm);
	fitfile_b    = opt2fn("-sb", NFILE, fnm);
	trxfile_a    = opt2fn("-fa", NFILE, fnm);
	trxfile_b    = opt2fn("-fb", NFILE, fnm);
	ndxfile_a    = opt2fn_null("-na", NFILE, fnm);
	ndxfile_b    = opt2fn_null("-nb", NFILE, fnm);
	///////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// control the output files //////////////////////////////
	resfile    = opt2fn_null("-res", NFILE, fnm);
	respdb    = opt2fn_null("-respdb", NFILE, fnm);
	projfile_a   = opt2fn_null("-proj_a", NFILE, fnm);
	projfile_b   = opt2fn_null("-proj_b", NFILE, fnm);
	dynamicspdb = opt2fn_null("-dyn", NFILE, fnm);
	if (!resfile && !respdb && !projfile_a && ! projfile_b && ! dynamicspdb)
	{
		fprintf(stderr,"Nothing to do!\n");
		return 0;
	}
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
  	///////////////////////////////////////////////////////////////////////////////////////
	////////////////////       read index files and do the selections       ///////////////
	///////////////////////////////////////////////////////////////////////////////////////
	int     nfit, nfit_b ;
    atom_id   *ifit_a, *ifit_b;
	char  *fitname_a,*fitname_b ;
	printf("\n Choose a group for the covariance analysis( and least squares fit) of state A: \n");
   	get_index(atoms_a, ndxfile_a, 1, &nfit, &ifit_a, &fitname_a);
	printf("\n Choose a group for the covariance analysis( and least squares fit) of state B \n");
	printf("\n This goup should be the same which was selected for state A: \n");
	get_index(atoms_b, ndxfile_b, 1, &nfit_b, &ifit_b, &fitname_b);
	///////////////////////////////////////////////////////////////////////////////////
	////////// check the correspondence between the selected groups in A & B //////////
	printf("#############################################\n");
    printf("Checking the correspondence between the selected \nresidues and atoms of the two states... \n");
	if (nfit != nfit_b)
	{
		fprintf(stderr,"Error:Numbers of selected atoms in A and B are different!\n");
		exit(EXIT_FAILURE);	
	} 
    if (nfit < 3)
	{
		fprintf(stderr,"Error:Need >= 3 points to fit!\n");
		exit(EXIT_FAILURE);
	}
	int a,b ,r1, r2;
	for(i=0; i < nfit; i++)
	{
		a= ifit_a[i];
		r1= atoms_a->atom[ a ].resind;
		b= ifit_b[i];
		r2= atoms_a->atom[ b ].resind;
		if ( (strcmp(*(atoms_a->atomname[a]) ,*(atoms_b->atomname[b]) ) != 0) ||
			 (strcmp( *(atoms_a->resinfo[r1].name),*(atoms_b->resinfo[r2].name)) != 0)	)
		{
			fprintf(stderr,"Error: Residue:atom %s:%s form state A is not corresponding to %s:%s from state B!\n",
							*(atoms_a->atomname[a]),*(atoms_a->resinfo[r1].name),*(atoms_b->atomname[b]), *(atoms_b->resinfo[r2].name));
			exit(EXIT_FAILURE);
		}
	}


	printf("Correspondence is OK!\n" );
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////      RPCA  analysis     ///////////////////////////////////
	printf("#############################################\n");
	ONE_DIM_ARRAY(eig, Geigdecomposition_t,1 );
	printf("*********************************************\n");
	printf("Reading the Generalized eigenpairs and the means\nfrom the file %s\n", eigvecfile);
	if (!read_Geigen(eigvecfile,eig ) )
	{
		fprintf(stderr," Can not read G eigenpairs from the file %s \n",eigvecfile );
		exit(EXIT_FAILURE);
	}
	if(nfit !=((eig->ndim)/3))
	{
		fprintf(stderr,"Error: Nr of dims in the eigenpairs file %s is different from the Nr. in the selected group of the structure files\n",eigvecfile);
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
	printf("#############################################\n");

/**********************
int r = eig->rank;
int n = eig->ndim;
ONE_DIM_ARRAY(WtXW, double, (r*r));
quadratic_Mat_multiply( eig->mat_a, n  ,   eig->gevec , r, WtXW );
int id = IsIdentity(WtXW,  r, 0.0000001, 0);
for (i = 0; i <(r*r); i++) WtXW[i]=0;
quadratic_Mat_multiply( eig->mat_b, n  ,  eig->gevec , r , WtXW );
int diag = IsDiagonal( WtXW,  r, 0.0000001 , 0) ;
int tt=1;
for (i = 0; i <r; i++)
{
	if(fabs( WtXW[i*r +i]-( eig->geigval)[i]) >  0.0000001)
	{
	printf("geval %d differs:  %f <-->  %f \n", i+1, WtXW[i*r +i],( eig->geigval)[i]  );
	tt=0;
	}
	}
if (tt && diag && id) printf("Simultaneous Diagonalization is fine!\n G^t*A*G=I \n G^t*B*G=diag(geval) \n");
else printf("Simultaneous Diagonalization is NOT fine!\n");
free(WtXW);
*******************************/

//write_Geigen2xvg("DK2.xvg", eig )	;
	printf("#############################################\n");
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
		RPCA_project(trxfile_a,TRUE, projfile_a ,eig,  Beg_evec,End_evec,nfit,ifit_a,bPBC,gpbc_b,  timeBeg, timeEnd, skip);
		RPCA_project(trxfile_b,FALSE, projfile_b ,eig,  Beg_evec,End_evec,nfit,ifit_b,bPBC,gpbc_b, timeBeg, timeEnd, skip);
	}
	printf("#############################################\n");
	if(resfile && covfile_b)
	{
		int n_tmp;
		printf("#############################################\n");
		printf("Reading the mean & the covariance matrix\n of state B from the file:\n %s\n",covfile_b );
		if (  ! read_cov_matrix(covfile_b ,&(eig->mat_b),  &(eig->m_b), &n_tmp)  )
		{
			fprintf(stderr,"Error in reading the mean and the covariance matrix of state B \n");
			exit(EXIT_FAILURE);
		}
		printf("Done\n");
		if(n_tmp != eig->ndim)
		printf("Computing the interaction map\n");
		interaction_map(resfile, eig,atoms_b,ifit_b,Beg_evec, End_evec,respdb);
	}
	else if(resfile && !covfile_b)
		fprintf(stderr,"Error: the requested interaction maps can not be performed without providing the needed covariance matrix B \n");
	//////// this is for the example in the paper
		/***************************
			printf("Project state b on its pca\n");
			int n_tmp;
			read_cov_matrix(covfile_b ,&(eig->mat_b),  &(eig->m_b), &n_tmp);
			eigdecomposition_t eigb;
			eigen(eig->mat_b, eig->ndim, &eigb);
			write_array2xvg("pca_eigenval.xvg" ,eigb.ndim,eigb.eval,"Eigenvalues of the covariance matrix" ,"Eigenvector index","(nm\\S2\\N)" );
			for (i = 0; i <20; i++) printf("%d %f\n", d+1, eigb.eval[i]);
			eig->gevec= eigb.vec ;
			RPCA_project(trxfile_a,TRUE, "PCA_PROJ_A.dat" ,eig,  Beg_evec,End_evec,nfit_b,ifit_b,bPBC,gpbc_b,  timeBeg, timeEnd, skip);
			RPCA_project(trxfile_b,FALSE, "PCA_PROJ_B.dat" ,eig,  Beg_evec,End_evec,nfit_b,ifit_b,bPBC,gpbc_b, timeBeg, timeEnd, skip);
		****************************************/
		printf("show the dynamical changes\n");
//show_dynamic_changes("file.pdb",eig,nfit,ifit_b,0.1,atoms_b );
show_increased_decresed_motions("dynamic",eig,nfit,ifit_b,tol,nv,nv,atoms_b, sep );
//show_motions_along_vectors("interpolate",eig,nfit,ifit_b,0,nv,atoms_b, TRUE );   // write the structures to the file

return 0;
}


