#ifndef _CWFITTING
#define _CWFITTING

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "statutil.h"
#include "typedefs.h"
#include "vec.h"
#include "rmpbc.h"


//////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
/// form total rotation matrix given the rotation angles around x, y and z axis
 double mah_dist(  rvec *xav, rvec *x ,int nfit ,atom_id  ifit[] ,double * ginv);

  double diff_frames_mah_dist(  rvec *xav, rvec *x ,int nfit ,atom_id ifit_av[] , atom_id ifit_x[] ,double * ginv);
 
void rotmatxyz(double  Qx, double Qy, double Qz, double **R);
/////////////////////////////////////////////////////////////////////
typedef struct {
  int nfit;
  double  **x;   
  double  **m ;
  double  *ginv;
} farguments_t;
/////////////////////
void rot_trans_frame( rvec *x, int np,double *vtr);

////////////////////////////////////////////////////////////////////
void for_minimizer_mah_rot_tran( int nvar , double *v ,  double *mah_dist,farguments_t *state);

///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////////////  fit x to xref using covariance weighted fitting  //////////////
//////////////////////////////////////////////////////////////////////////////
void CW_fitting(			rvec *x,
					atom_id 	*ifit_x,      // atoms in x for fitting
					rvec *xref,                   // reference structure used  for fitting (usually the average) 
					atom_id 	*ifit_r,      // atoms in xref for fitting
		 			int 		nfit,
		 			int natoms,      // nr atoms in x
		                        double      * ginv,   // generalized inverse of the covariance matrix
					int 	verbose     // print the change of MAH dist
				);
 
 
 //////////////////////////////////////////////////////////////////////////////
void Covariance_CW_fitting_traj(	
						const char 	*trxfile,   // traj file to read
						const char 	*out_file, // if not NULL , the fitted traj will be written to the file out_file
						rvec  		*xav,      // average structure used as reference for fitting 
						double      * ginv,   // generalized inverse of the covariance matrix
						int 		nfit,
						atom_id 	*ifit,      // atoms for fitting
						int 		natoms,  // nr atoms to write
						atom_id 	*index, /// atoms to write
						gmx_bool 	bPBC,
						gmx_rmpbc_t     gpbc,
						t_atoms        *atoms,
						double          *mat,  // output the new covariance matrix (nfit*3 *nfit*3 ) after CW fitting
						double           *m,   // the computed average after cw fitting 
						atom_id 	*ifit_xav ,     // atoms for fitting, NULL means all the atoms in xav are taken as a reference
						int verbose,
						float first_time,
						float last_time);
			 
			 
/////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 //////   apply covariance weigthed fitting and write the trajactory to the file out_file
 ////////////////////////////////////////////////////////////////
void CW_fitting_traj(	
						const char 	*trxfile,   // traj file to read
						const char 	*out_file, //  the fitted traj will be written to the file out_file
						rvec  		*ref,      // average structure used as refrence for fitting 
						double      * ginv,   // generalized inverse of the covariance matrix
						int 		nfit,
						atom_id 	*ifit,      // atoms for fitting
						int 		natoms,  // nr atoms to write
						atom_id 	*index, /// atoms to write
						gmx_bool 	bPBC,
						gmx_rmpbc_t     gpbc,
						output_env_t    oenv,
						t_atoms        *atoms,
						atom_id 	*ifit_ref,      // atoms for fitting if refrence is different
						float first_time,
						float last_time
					);





float log_Eval(double *eigval , int n );


//void nlmin_slim( fcn_p fcn , int n , double *x, void *state, double *fpls ,int itnlim ,  int verbose );

/////////////////////////////////////////////////////////////////////

#endif	/* _CWFITTING */
