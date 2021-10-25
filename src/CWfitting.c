////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
#include <float.h>
#include <math.h>
#include <cblas.h>
#include "linear_algebra.h"
#include "fitting.h"
#include "arrays.h"
#include "uncmin.h"
#include "CWfitting.h"
#include "xdrfile_xtc.h"


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// compute mah. distance b. two frames x and xav
// using the indexed atoms in ifit and the generalized inverse of the 
//covariance mat. of the indexed atoms
//////////////////////////////////////////////////////////////////////////////
 double mah_dist(  rvec *xav, rvec *x ,int nfit ,atom_id  ifit[] ,double * ginv)
 {
	ONE_DIM_ARRAY(dx, double, (nfit*3) );
	int i,ndx;
	int ndim = nfit*3;
	for (i = 0; i < nfit; i++)
	{
		ndx= ifit[i];
		dx[3*i]    = x[ndx][0] - xav[ndx][0] ;
		dx[3*i +1] = x[ndx][1] - xav[ndx][1];
		dx[3*i +2] = x[ndx][2] - xav[ndx][2];
	}
	double d = mah_norm( ginv, ndim  , dx  ) ;
	//printf("inside:   %f    %f    %f \n", dx[0], dx[1] , dx[2]);
	free(dx);
	return d;
	
 }
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// compute mah dist b. two frames x and xav
// using the indexed atoms in ifit and the ginverse of the covariance mat
// for the indexed atoms
//////////////////////////////////////////////////////////////////////////////
 double diff_frames_mah_dist(  rvec *xav, rvec *x ,int nfit ,atom_id ifit_av[] , atom_id ifit_x[] ,double * ginv)
 {
	ONE_DIM_ARRAY(dx, double, (nfit*3) );
	int i,ndx,ndv;
	int ndim = nfit*3;
	for (i = 0; i < nfit; i++)
	{
		ndx= ifit_x[i];
		ndv=ifit_av[i];
		dx[3*i]    = x[ndx][0] - xav[ndv][0] ;
		dx[3*i +1] = x[ndx][1] - xav[ndv][1];
		dx[3*i +2] = x[ndx][2] - xav[ndv][2];
	}
	double d = mah_norm( ginv, ndim  , dx  ) ;
	free(dx);
	return d;
	
 }
////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///// rotate the frame x using angles in vtr[0-2] and translate using the
///// displacement in vtr[3-5]
///// np is the number of atoms in the frame
//////////////////////////////////////////////////////////////////////////////
void rot_trans_frame( rvec *x, int np,double *vtr)
{
 
	int i;
	ONE_DIM_ARRAY( v, float, 6);
	for (i = 0; i < 6; i++) v[i] = (float) vtr[i];
	TWO_DIM_ARRAY(R, float, 3, 3);
	float sx = sinf(v[0]); float sy = sinf(v[1]);	float sz = sinf(v[2]);
	float cx = cosf(v[0]); float cy = cosf(v[1]);	float cz = cosf(v[2]);
	//// rotation mat
	R[0][0] =  cy*cz ;    R[0][1] =  sx*sy*cz -cx*sz;           R[0][2] =  sx*sz + cx*sy*cz;
	R[1][0] =  cy*sz ;    R[1][1] =  cx*cz  + sx*sy*sz;         R[1][2] =  cx*sy*sz - sx*cz;
	R[2][0] =  -sy   ;    R[2][1] =  sx*cy;                     R[2][2] =  cx*cy;
	float Rx ,Ry ,Rz ;
	for (i = 0; i < np; i++)
	{
		Rx = (float)v[3] + R[0][0]*x[i][0] + R[0][1]*x[i][1] + R[0][2]*x[i][2] ;
		Ry = (float)v[4] + R[1][0]*x[i][0] + R[1][1]*x[i][1] + R[1][2]*x[i][2] ;
		Rz = (float)v[5] + R[2][0]*x[i][0] + R[2][1]*x[i][1] + R[2][2]*x[i][2] ;
		x[i][0] = Rx;
		x[i][1] = Ry;
		x[i][2] = Rz;
	}
	free(R);
	free(v);
}
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/// form the total rotation matrix given the rotation angles around x, y and z axis
/***************************************************************************************
total rotation matrix
http://danceswithcode.net/engineeringnotes/rotations_in_3d/rotations_in_3d_part1.html
	    Rx 	1,       0 ,      0 ,
			0, cos(Qx) , -sin(Qx),
			0, sin(Qx), cos(Qx)
					   
	    Ry  cos(Qy),  0,  sin(Qy),
			0,        1,         0,
			-sin(Qy), 0,   cos(Qy)
				   
	   Rz   cos(Qz),  -sin(Qz),   0,
			sin(Qz),  cos(Qz),   0,
			0      ,        0,   1
				
		   cy*cz         sx*sy*cz -cx*sz     sx*sz + cx*sy*cz
R=RzRyRx=  cy*sz         cx*cz + sx*sy*sz   cx*sy*sz - sx*cz
			-sy      		 sx*cy                  cx*cy
***************************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void rotmatxyz(double  Qx, double Qy, double Qz, double **R)
{

	double sx = sin(Qx);     
	double sy = sin(Qy);      
	double sz = sin(Qz);
	double cx = cos(Qx); 
	double cy = cos(Qy);
	double cz = cos(Qz);

	//// rotation mat
	R[0][0] =  cy*cz ;    R[0][1] =  sx*sy*cz -cx*sz;           R[0][2] =  sx*sz + cx*sy*cz;
	R[1][0] =  cy*sz ;    R[1][1] =  cx*cz  + sx*sy*sz;         R[1][2] =  cx*sy*sz - sx*cz;
	R[2][0] =  -sy   ;    R[2][1] =  sx*cy;                     R[2][2] =  cx*cy;
	
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////// this function is used by the minimizer 
//////// displacement of x, y & z are in v[3-5]  
//////// rotation angles Qx, Qy & Qz are in  v[0-2]
//////// Mah. dist. is computed after the translation and rotation
//////////////////////////////////////////////////////////////////////////////
  void for_minimizer_mah_rot_tran( int n , double *v ,  double *mah_dist,farguments_t *state)
 {
	int i;
	TWO_DIM_ARRAY(R, double, 3, 3);
	rotmatxyz(v[0], v[1], v[2],  R);
	double xn , yn ,zn ;
	double **x = state->x;
	double **xav = state->m;
	ONE_DIM_ARRAY(dx, double, (state->nfit*3) );
	for (i = 0; i < state->nfit; i++)
	{
		xn = v[3] + R[0][0]*x[i][0] + R[0][1]*x[i][1] + R[0][2]*x[i][2] ;
		yn = v[4] + R[1][0]*x[i][0] + R[1][1]*x[i][1] + R[1][2]*x[i][2] ;
		zn = v[5] + R[2][0]*x[i][0] + R[2][1]*x[i][1] + R[2][2]*x[i][2] ;
		////  sub. the mean	
		dx[3*i]    =  xn - xav[i][0] ;
		dx[3*i +1] =  yn - xav[i][1];
		dx[3*i +2] =  zn - xav[i][2];
	}

	*mah_dist= mah_norm( state->ginv, state->nfit*3  , dx  ) ;
	free(dx);
	free(R);
 }
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////////////  fit x to xref using covariance weighted fitting  //////////////
//////////////////////////////////////////////////////////////////////////////
void CW_fitting(rvec 		*x,
				atom_id 	*ifit_x,	// atoms in x for fitting
				rvec		*xref,                   // reference structure used  for fitting (usually the average) 
				atom_id 	*ifit_r,    // atoms in xref for fitting
		 		int 		nfit,
		 		int 		natoms,     // nr atoms in x
		        double      * ginv,		// generalized inverse of the covariance matrix
				int 		verbose     // print the change of MAH dist
				)
{

	
	int  i,d, ndx, ndr;
	double fpls,dbefore ,*dx;
	farguments_t state;
	state.nfit= nfit;
	state.ginv = ginv ;
	ALLO_TWO_DIM_ARRAY(state.m,double,nfit, 3) ;
	ALLO_TWO_DIM_ARRAY(state.x,double,nfit, 3) ;
	///////
	if(ifit_x)
	{
	  for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++) state.x[i][d]= (double) x[ifit_x[i]][d];
	}
	else
	{
	  for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++) state.x[i][d]= (double) x[i][d];
	}
	/////
	if(ifit_r)
	{
	  for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++)  state.m[i][d]= (double) xref[ifit_r[i]][d];
	}
	else
	{
	  for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++)  state.m[i][d]= (double) xref[i][d];
	}
	///////
	if (verbose)
	{
	  ALLO_ONE_DIM_ARRAY(dx,double,(nfit*3)) ;
	  for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++) dx[3*i+d]  = state.x[i][d] -state.m[i][d] ;
	  dbefore = mah_norm( ginv, (nfit*3)  , dx  ) ;
	}

	ONE_DIM_ARRAY( v, double, 6);
	for (i = 0; i < 6; i++) v[i]=0.0;
	//nlmin( (fcn_p)for_minimizer_mah_rot_tran, 6 , v, &state, &fpls,500 , 1,  0);
	nlmin_slim( (fcn_p)for_minimizer_mah_rot_tran, 6 , v, &state, &fpls,500 , verbose);
	rot_trans_frame(x,natoms,v);
	if (verbose)  printf("change of Mahalanobis distance after the\n CW fitting= %f -%f= %f\n", dbefore,fpls ,dbefore-fpls);
	free(v);
	FREE_TWO_DIM_ARRAY(state.m,nfit);
	FREE_TWO_DIM_ARRAY(state.x,nfit);

}

///////////////////////////////////////////////////////////////////////////
 
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*****************************************************************************
Apply covariance weighted fitting of the conformations in a traj.
 to the reference (the average) xav:
 the fitted trajectory are written if out_file is given 
 mat != NULL:computation of covariance matrix and and average after the fitting
 *****************************************************************************/
 
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
						atom_id 	*ifit_xav,      // atoms for fitting, NULL means all the atoms in xav are taken as a reference
						int verbose,
						float first_time,
						float last_time
							)
{
	printf("Computing covariance matrix of traj %s using CW fitting....\n",trxfile );
	XDRFILE *xd;
	if (out_file)
	{
		printf("saving the fitted traj to file %s\n",out_file );
		xd = xdrfile_open(out_file,"w");
		if (NULL == xd)	PRINT_ERR_OPENFILE( out_file);
	}
	int  i,d,j,ndx;
	real t;
	double fpls;
	farguments_t state;
	state.nfit= nfit;
	int ndim= nfit*3;
	ALLO_TWO_DIM_ARRAY(state.m,double,nfit, 3) ;
	double dbefore, dafter, *dx;
	if (verbose >1 )	  ALLO_ONE_DIM_ARRAY(dx,double,(nfit*3)) ;
	for (i = 0; i < nfit; i++)
	{
		if(ifit_xav) ndx= ifit_xav[i];
		else ndx= i;
		state.m[i][0]= (double) xav[ndx][0];
		state.m[i][1]= (double) xav[ndx][1];
		state.m[i][2]= (double) xav[ndx][2];
	}
	ALLO_TWO_DIM_ARRAY(state.x,double,nfit, 3) ;
   	state.ginv = ginv ;
	ONE_DIM_ARRAY( v, double, 6);
	ONE_DIM_ARRAY( xt, double, ndim);
	for (i=0; i< ndim*ndim ; i++) mat[i]=0;
	for(i=0; i < ndim; i++) m[i]=0;
	xtc_t f;
	Read_first_xtc( trxfile,&f, first_time );
	int nframes=0;
	double TMHB=0;
    do
    {
		nframes++;
		if (bPBC) gmx_rmpbc(gpbc, f.natoms, f.box, f.x);
		center_x2ref_diffrenct_atoms(f.x, ifit,xav ,ifit_xav, nfit, f.natoms );
		rotate_B2R(f.x, ifit,xav ,ifit_xav, nfit, f.natoms );

		for (i = 0; i < nfit; i++)
		{
			ndx= ifit[i];
			state.x[i][0]= (double) f.x[ndx][0];
			state.x[i][1]= (double) f.x[ndx][1];
			state.x[i][2]= (double) f.x[ndx][2];
		}

		if (verbose>1)
		{
			for (i = 0; i < nfit; i++) for (d = 0; d < 3; d++) dx[3*i+d]  = state.x[i][d] -state.m[i][d] ;
			dbefore = mah_norm( ginv, (nfit*3)  , dx  ) ;
		}
		for (i = 0; i < 6; i++) v[i]=0.0;
		nlmin_slim( (fcn_p)for_minimizer_mah_rot_tran, 6 , v, &state, &fpls,500 , 0);
		TMHB +=fpls;
		rot_trans_frame(f.x,f.natoms,v);
		if(verbose>1)
		{
			printf("change of Mahalanobis distance after the CW fitting= %f -%f= %f\n", dbefore,fpls ,dbefore-fpls);
		}
		for (i = 0; i < nfit; i++)
		{
			ndx=ifit[i];
			xt[i*3]=(double)f.x[ndx][0] ;
			xt[i*3+1]=(double)f.x[ndx][1] ;
			xt[i*3+2]=(double)f.x[ndx][2] ;
		}
		cblas_daxpy (ndim, 1, xt, 1,m, 1);
		cblas_dsyr( CblasColMajor ,CblasLower, ndim ,1, xt ,1,mat ,ndim );
		if (out_file)	write_xtc(xd,f.natoms,f.step,f.time, f.box,f.x);
	}
    while( Read_next_xtc(&f, last_time));
    reset_frame( &f);
	double inv_nframes = 1.0/ ((double) nframes);
  	cblas_dscal(ndim, inv_nframes, m, 1);
	cblas_dscal(ndim*ndim , inv_nframes, mat, 1);
    cblas_dsyr( CblasColMajor ,CblasLower, ndim ,-1, m ,1,mat ,ndim );
	//// fill the upper part of the sym mat
	for (i=0; i< ndim; i++) for (j=i; j< ndim; j++) mat[j*ndim +i]=mat[i*ndim +j] ;
  	if(verbose) printf("Average Mahalnobis distances= %f\n", inv_nframes * TMHB );
  	if (out_file) xdrfile_close(xd);
	free(v);
	free(xt);
	if (verbose>1 ) free(dx);
}
	
//////////////////////////////////////////////////////////////////////////////////
 
///////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 //////   apply covariance weigthed fitting and write the trajactory to the file out_file
 ////////////////////////////////////////////////////////////////
void CW_fitting_traj(	
						const char 	*trxfile,   // traj file to read
						const char 	*out_file, //  the fitted traj will be written to the file out_file
						rvec  		*ref,      // average structure used as a reference for fitting 
						double      * ginv,   // generalized inverse of the covariance matrix
						int 		nfit,
						atom_id 	*ifit,      // atoms for fitting
						int 		natoms,  // nr atoms to write
						atom_id 	*index, /// atoms to write
						gmx_bool 	bPBC,
						gmx_rmpbc_t     gpbc,
						output_env_t    oenv,
						t_atoms        *atoms,
						atom_id 	*ifit_ref,      // atoms for fitting if reference is different
						float first_time,
						float last_time
					)
{

	XDRFILE *xd= xdrfile_open(out_file,"w");
	if (NULL == xd) PRINT_ERR_OPENFILE( out_file);
	int i, ndx;
	double fpls;
	farguments_t state;
	state.nfit= nfit;
	int ndim= nfit*3;
	ALLO_TWO_DIM_ARRAY(state.m,double,nfit, 3) ;
	for (i = 0; i < nfit; i++)
	{
		if(ifit_ref) ndx= ifit_ref[i];
		else ndx= ifit[i];
		state.m[i][0]= (double) ref[ndx][0];
		state.m[i][1]= (double) ref[ndx][1];
		state.m[i][2]= (double) ref[ndx][2];
	}
	ALLO_TWO_DIM_ARRAY(state.x,double,nfit, 3) ;
   	state.ginv = ginv ;
	ONE_DIM_ARRAY( v, double, 6);
	//ONE_DIM_ARRAY( xt, double, ndim);
	xtc_t f;
	Read_first_xtc( trxfile,&f, first_time );
	int nframes=0;
    do
    {
		nframes++;
		if (bPBC) gmx_rmpbc(gpbc, f.natoms, f.box, f.x);
		if(ref )
		{
			
		center_x2ref_diffrenct_atoms(f.x, ifit,ref ,ifit_ref, nfit, f.natoms );
		rotate_B2R(f.x, ifit,ref ,ifit_ref, nfit, f.natoms );
		}
		for (i = 0; i < nfit; i++)
		{
			ndx= ifit[i];
			state.x[i][0]= (double) f.x[ndx][0];
			state.x[i][1]= (double) f.x[ndx][1];
			state.x[i][2]= (double) f.x[ndx][2];
		}
		for (i = 0; i < 6; i++) v[i]=0.0;
		nlmin_slim( (fcn_p)for_minimizer_mah_rot_tran, 6 , v, &state, &fpls,500 , 0);
		rot_trans_frame(f.x,f.natoms,v);
		write_xtc(xd,f.natoms,f.step,f.time, f.box,f.x);
	}
    while( Read_next_xtc(&f, last_time));
    xdrfile_close(xd);
    reset_frame( &f);
  	free(v);
}
	
//////////////////////////////////////////////////////////////////////////////////
 


/////////////////////////////////////////////////////////////////////////////////
float log_Eval(double *eigval , int n )
{
	float tol =  FLT_EPSILON  * eigval[0] ;
	if (tol < 0.0) tol= FLT_EPSILON ;
	float ll=0 ;
	int i;
	for (i = 0; i < n; i++) 
	{
		if(eigval[i] > tol) ll -=log(eigval[i]); 
	}
	return ll;
}

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


	
		
