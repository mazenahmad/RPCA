#ifndef _FITTING
#define _FITTING

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "vec.h"
#include "rmpbc.h"


//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO the selected 
//// atoms in x is th same cnter as th slected atoms of ref
///// atoms are  assumed to be the same in x and ref
void center_x(rvec x[] ,atom_id ind_cm[],int nfit, int natom, rvec COG);

//////////////////////////////////////////////////////////

void get_COG(rvec ref[] ,atom_id ind_ref[],int nfit, rvec cog );

//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO the selected 
//// atoms in x is th same cnter as th slected atoms of ref
///// atoms are  assumed to be the same in x and ref 
void center_x2ref(rvec x[],rvec ref[] , atom_id ind[],int nfit, int natom );
////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO the selected 
//// atoms in x is th same cnter as th slected atoms of ref
///// atoms are not assumed to be the same in x and ref 
///// all atoms are taken if the coresponding index is NULL 
void center_x2ref_diffrenct_atoms(rvec x[], atom_id ind_x[],rvec ref[] ,atom_id ind_ref[],int nfit, int natom_x );
////////////////////////////////////////
//// print the center of geometry of the selected atoms
void print_cog( rvec x[],atom_id ind[],int nfit);

//////////////////////////////////////////////////////////
	
void jacobi(double **a,int n,double d[],double **v,int *nrot);

/////////////////////////////////////////////////////////
void calc_rotation_R(int ndim,int natoms,rvec *xp,rvec *x,matrix R);

//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//// rotate x based on the  fitting of the slected atoms ind_fit
//// Note: rotation matrix is calculated for the selected while all the atoms are rotated
void rotate_x( int natoms,rvec * xref,rvec *x , int nfit, atom_id ind_fit[]);

/////////////////////////////////////////////////////////
//// fit  x_ b to  ref using the indexed atoms ind_ref in ref and  ind_fit_b in x_b
//// differ from rotate_x that it can be used with diffrent nr of atoms in ref and b
//// Note: rotation matrix is calculated for the selected while all the the atoms of x_b are rotated
///// all atoms are taken if the coresponding index is NULL 

void rotate_B2R(rvec *x_b, atom_id ind_fit_b[], rvec *ref, atom_id ind_ref[] ,int nfit, int natoms_b );




/////////////////////////////////////////////////////////
real calc_rmsd( rvec *x , atom_id ind_fit[] ,rvec * xref,  atom_id ind_fit_ref[] ,int nfit);

////////////////////////////////////////////////////////
void gpa(const char 	*trxfile,   // traj file to read
			rvec  		**xav,       // coor on the average 
			int 		nfit,
			atom_id 	*ifit,
			gmx_bool 	bPBC,
			gmx_rmpbc_t     gpbc,
			int 			TC,
			real 			convergence,
		      int verbose,
			  float first_time,
			  float last_time
						);


			  

/////////////////////////////////////////////////////////			  
void write_fitted_traj(	
						const char 	*trxfile,   // traj file to read
						const char 	*out_file,
						rvec  		*xref,      // refrence coor structure to be used for fitting 
						int 		nfit,
						atom_id 	*ifit,      // atoms for fitting
						int 		natoms,  // nr atoms to write
						atom_id 	*index, /// atoms to write
						gmx_bool 	bPBC,
						gmx_rmpbc_t     gpbc,
						float first_time,
						float last_time
					);			  
			  
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////// write a frame to a one col. data file
////////////////////////////////////////////////////////////////////////
void	write_frame2datafiles ( const char *fileData,rvec *x,int nfit, atom_id 	*ifit,t_atoms *atoms);





#endif	/* _FITTING */
