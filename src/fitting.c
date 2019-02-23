////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
#include <float.h>
#include <math.h>
#include <cblas.h>
#include "fitting.h"
#include "arrays.h"
#include "xdrfile_xtc.h"
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO the selected 
////// group define by the index ind_cm zero
////// THE COG if set to COG if not NULL
void center_x(rvec x[] ,atom_id ind_cm[],int nfit, int natom, rvec COG)
{
	int  i,m;
	real xcm[3]= {0.0, 0.0, 0.0}; 
	real tm=(real) nfit;
	
	if(ind_cm)
	{
		for(i=0; i<nfit; i++) for(m=0; m < 3; m++)  xcm[m]+=x[ind_cm[i]][m];
	}
	else
	{
		for(i=0; i<nfit; i++) for(m=0; m < 3; m++) xcm[m] +=  x[i][m];
	}
	for(m=0; m < 3; m++)  xcm[m]/=tm;
	if(COG) for(m=0; m < 3; m++)  xcm[m] += COG[m];
	for(i=0; i< natom; i++) for(m=0; m < 3; m++) x[i][m] -= xcm[m];
	
}
//////////////////////////////////////////////////////////

void get_COG(rvec ref[] ,atom_id ind_ref[],int nfit, rvec cog )

{
	int  i,m;
	cog[0]=0.0; cog[1]=0.0; cog[2]=0.0;
	if(ind_ref)
	{
		for(i=0; i< nfit; i++)	for(m=0; m < 3; m++)  cog[m]+= ref[ind_ref[i]][m];
	}
	else
	{
		for(i=0; i< nfit; i++)	for(m=0; m < 3; m++)  cog[m]+= ref[i][m];
	}

	for(m=0; m < 3; m++) cog[m]/= (real)(nfit);
}
///////////////////////////////////////////////
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO the selected 
//// atoms in x is th same cnter as th slected atoms of ref
///// atoms are  assumed to be the same in x and ref 
void center_x2ref(rvec x[],rvec ref[] , atom_id ind[],int nfit, int natom )
{
  int  i,m, xi;
  real trans[3] = {0.0, 0.0 , 0.0};
  //rvec trans;
  real tm;
  tm=(real) nfit;
  //clear_rvec(trans);
	if (ind)
	{
		for(i=0; i< nfit; i++) 
		{	
			xi=ind[i];
			for(m=0; m < 3; m++)  trans[m]+= ref[xi][m];
			for(m=0; m < 3; m++)  trans[m]-= x[xi][m];
		}
	}
	else
	{	
		for(i=0; i< nfit; i++) 
		{
			for(m=0; m < 3; m++)  trans[m]+= ref[i][m];
			for(m=0; m < 3; m++)  trans[m]-= x[i][m];
		}
	}
	for(m=0; m < 3; m++) trans[m]/=tm;
	for(i=0; i< natom; i++)  for(m=0; m < 3; m++)  x[i][m] += trans[m];
	////////////////////////////////////
	/***************
	printf("\n COM before= %g  %g  %g \n",trans[0]/tm, trans[1]/tm,trans[2]/tm );
	 clear_rvec(trans);

	  for(i=0; i<nfit; i++) 
  {	
	xi=ind[i];
    for(m=0; m < 3; m++)  trans[m]+=ref[xi][m];
	for(m=0; m < 3; m++)  trans[m]-=x[xi][m];
  }
 	printf("\n COM af= %g  %g  %g \n",trans[0]/tm, trans[1]/tm,trans[2]/tm );
   ****************/
}


////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////// translate x to make the centre of GEO of the selected 
//// atoms in x is the same as th slected atoms of ref
///// atoms are not assumed to be the same in x and ref 
///// all atoms are taken if the coresponding index is NULL 
void center_x2ref_diffrenct_atoms(rvec x[], atom_id ind_x[],rvec ref[] ,atom_id ind_ref[],int nfit, int natom_x )

{
	int  i,m;
	rvec trans;
	clear_rvec(trans);
	
	if (ind_x)
	{
		for(i=0; i< nfit; i++)	for(m=0; m < 3; m++)  trans[m]-= x[ind_x[i]][m];
	}
	else
	{
		for(i=0; i< nfit; i++) for(m=0; m < 3; m++)  trans[m]-= x[i][m];
	}

	//////////////
	real tm = 0.0;
	if(ind_ref)
	{
		for(i=0; i< nfit; i++)
		{
			for(m=0; m < 3; m++)  trans[m]+= ref[ind_ref[i]][m];
			tm +=1;
		}
		
	}
	else
	{
		for(i=0; i< nfit; i++)
		{
			for(m=0; m < 3; m++)  trans[m]+= ref[i][m];
			tm +=1;
		}
	}
	///////////////
	
	for(m=0; m < 3; m++) trans[m]/= tm;
	if (!ind_x && (natom_x != nfit)) PRINT_FATAL_ERR_MES("Error nr of atoms in the frame is diffrent form the nr of fitted atoms\n");
	for(i=0; i< natom_x; i++)  for(m=0; m < 3; m++)  x[i][m] += trans[m];
	/***************
print_cog( x, ind_x, nfit);
print_cog( ref, ind_ref, nfit);
printf("--------------\n");
	///// test

	printf("--------------\n");
	print_cog( x, ind_x, nfit);
	print_cog( ref, ind_ref, nfit);
	printf("--------------\n");
	
	int xi;
    clear_rvec(trans);
    for(i=0; i<nfit; i++) 
  	{     
        	xi=ind_x[i];
    		for(m=0; m < 3; m++)  trans[m]+=ref[i][m];
        	for(m=0; m < 3; m++)  trans[m]-=x[xi][m];
  	}
        printf("\n COM at= %g  %g  %g \n",trans[0], trans[1],trans[2] );
 	 ****************/

}
//////////////////////////////////////////////////////////

void print_cog( rvec x[],atom_id ind[],int nfit)
{
	int i, m;
	real cog[3]={0,0,0}; 
	if(ind) for(i=0; i<nfit; i++) for(m=0; m < 3; m++) cog[m] += x[ind[i]][m];
	else    for(i=0; i<nfit; i++) for(m=0; m < 3; m++) cog[m] += x[i][m];
	 printf("COM at= %g  %g  %g \n",cog[0]/(real)nfit, cog[1]/(real)nfit,cog[2]/(real)nfit );
}
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);
	
void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
	int j,i ,iq,ip;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	ALLO_ONE_DIM_ARRAY(b, double ,n);
	ALLO_ONE_DIM_ARRAY(z, double ,n);
    for (ip=0; ip<n; ip++)
    {
    	for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
    	v[ip][ip]=1.0;
    }
    for (ip=0; ip<n;ip++)
    {
    	b[ip]=d[ip]=a[ip][ip];
    	z[ip]=0.0;
    }
    *nrot=0;
    for (i=1; i<=50; i++)
    {
    	sm=0.0;
    	for (ip=0; ip<n-1; ip++)
    	{
    		for (iq=ip+1; iq<n; iq++)     sm += fabs(a[ip][iq]);
    	}
    	if (sm == 0.0)
    	{
    		free(b);
    		free(z);
    		return ;
    	}
    	if (i < 4)   tresh=0.2*sm/(n*n);
    	else        tresh=0.0;
    	for (ip=0; ip<n-1; ip++)
    	{
    		for (iq=ip+1; iq<n; iq++)
    		{
    			g=100.0*fabs(a[ip][iq]);
    			if (i > 4 && fabs(d[ip])+g == fabs(d[ip]) && fabs(d[iq])+g == fabs(d[iq]))
    					a[ip][iq]=0.0;
    			else if (fabs(a[ip][iq]) > tresh)
    			{
    				h=d[iq]-d[ip];
    				if (fabs(h)+g == fabs(h))	t=(a[ip][iq])/h;
    				else
    				{
    					theta=0.5*h/(a[ip][iq]);
    					t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
    					if (theta < 0.0) t = -t;
    				}
    				c=1.0/sqrt(1+t*t);
    				s=t*c;
    				tau=s/(1.0+c);
    				h=t*a[ip][iq];
    				z[ip] -= h;
    				z[iq] += h;
    				d[ip] -= h;
    				d[iq] += h;
    				a[ip][iq]=0.0;
    				for (j=0; j<ip; j++)
    				{
    					ROTATE(a,j,ip,j,iq)
    				}
    				for (j=ip+1; j<iq; j++)
    				{
    					ROTATE(a,ip,j,j,iq)
    				}
    				for (j=iq+1; j<n; j++)
    				{
    					ROTATE(a,ip,j,iq,j)
    				}
    				for (j=0; j<n; j++)
    				{
    					ROTATE(v,j,ip,j,iq)
    				}
    				++(*nrot);
    			}
    		}
    	}
    	for (ip=0; ip<n; ip++)
    	{
    		b[ip] +=  z[ip];
    		d[ip]  =  b[ip];
    		z[ip]  =  0.0;
    	}
    }
  PRINT_FATAL_ERR_MES("Error: Too many iterations in routine JACOBI");
}




/////////////////////////////////////////////////////////
void calc_rotation_R(int ndim,int natoms,rvec *xp,rvec *x,matrix R)
{
  int    c,r,n,j,m,i,irot,s;
  double **omega,**om;
  double d[2*DIM],xnr,xpc;
  matrix vh,vk,u;
  real   mn;
  int    index;
  real   max_d;

  if (ndim != 3 && ndim != 2)  PRINT_FATAL_ERR_MES("calc_fit_R called with ndim not 3 or 2");
  ALLO_TWO_DIM_ARRAY(omega, double,(2*ndim), (2*ndim));
  ALLO_TWO_DIM_ARRAY(om, double,(2*ndim), (2*ndim));
  for(i=0; i<2*ndim; i++) 
  {
    d[i]=0;
    for(j=0; j<2*ndim; j++) 
	{
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++)
    for(c=0; (c<ndim); c++) {
	xpc=xp[n][c];
	for(r=0; (r<ndim); r++) {
	  xnr=x[n][r];
	  u[c][r]+=xnr*xpc;
	}
      }
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<2*ndim; r++)
    for(c=0; c<=r; c++)
      if (r>=ndim && c<ndim) {
        omega[r][c]=u[r-ndim][c];
        omega[c][r]=u[r-ndim][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }

  /*determine h and k*/
  jacobi(omega,2*ndim,d,om,&irot);
  /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
   *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */
  
 // if (debug && irot==0) fprintf(debug,"IROT=0\n");
  
  index=0; /* For the compiler only */

  /* Copy only the first ndim-1 eigenvectors */  
  for(j=0; j<ndim-1; j++) 
  {
    max_d=-1000;
    for(i=0; i<2*ndim; i++)
      if (d[i]>max_d) 
	  {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<ndim; i++) 
	{
      vh[j][i]=M_SQRT2*om[i][index];
      vk[j][i]=M_SQRT2*om[i+ndim][index];
    }
  }
  if (ndim == 3) 
  {
    /* Calculate the last eigenvector as the outer-product of the first two.
     * This insures that the conformation is not mirrored and
     * prevents problems with completely flat reference structures.
     */  
    cprod(vh[0],vh[1],vh[2]);
    cprod(vk[0],vk[1],vk[2]);
  } 
  else if (ndim == 2) 
  {
    /* Calculate the last eigenvector from the first one */
    vh[1][XX] = -vh[0][YY];
    vh[1][YY] =  vh[0][XX];
    vk[1][XX] = -vk[0][YY];
    vk[1][YY] =  vk[0][XX];
  }

  /* determine R */
  clear_mat(R);
  for(r=0; r<ndim; r++)
    for(c=0; c<ndim; c++)
      for(s=0; s<ndim; s++)
	R[r][c] += vk[s][r]*vh[s][c];
  for(r=ndim; r<DIM; r++)
    R[r][r] = 1;
  FREE_TWO_DIM_ARRAY(omega,(2*ndim));
  FREE_TWO_DIM_ARRAY(om,(2*ndim));

}


//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//// rotate x based on the  fitting of the slected atoms ind_fit
//// Note: rotation matrix is calculated for the selected while all the atoms are rotated
void rotate_x( int natoms,rvec * xref,rvec *x , int nfit, atom_id ind_fit[])
{
  int    i,j,m,r,c, ai;
 
  matrix R;
  rvec   x_old;
  rvec   *fit_x, *fit_xref ;
  ALLO_ONE_DIM_ARRAY(fit_x, rvec ,nfit);
  ALLO_ONE_DIM_ARRAY(fit_xref, rvec ,nfit);
  //snew(fit_x,nfit);
  //snew(fit_xref,nfit);
	if (ind_fit)
	{
		for(i=0; i< nfit; i++) 
		{
			ai=ind_fit[i];
			for(m=0; m<3; m++)
			{
				fit_x[i][m]= x[ai][m];
				fit_xref[i][m]= xref[ai][m];
			}
		}
	}
	else
	{
		for(i=0; i< nfit; i++) 
		{
			for(m=0; m<3; m++)
			{
				fit_x[i][m]= x[i][m];
				fit_xref[i][m]= xref[i][m];
			}
		}
	}
	/* Calculate the rotation matrix R */
	calc_rotation_R(3,nfit,fit_xref,fit_x,R);
  /*rotate X*/
  for(j=0; j<natoms; j++) 
  {
    for(m=0; m<DIM; m++)  x_old[m]=x[j][m];
    for(r=0; r<DIM; r++) 
	{
      x[j][r]=0;
      for(c=0; c<DIM; c++)
      x[j][r]+=R[r][c]*x_old[c];
    }
  }
  //sfree(fit_x);
  //sfree(fit_xref);
  free(fit_x);
  free(fit_xref);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//// fit  x_ b to  ref using the indexed atoms ind_ref in ref and  ind_fit_b in x_b
//// differ from rotate_x that it can be used with diffrent nr of atoms in ref and b
//// Note: rotation matrix is calculated for the selected while all the the atoms of x_b are rotated
///// all atoms are taken if the coresponding index is NULL 

void rotate_B2R(rvec *x_b, atom_id ind_fit_b[], rvec * ref, atom_id ind_ref[] ,int nfit, int natoms_b )
{
	int    i,j,m,r,c, ai, bi;
 	matrix R;
	rvec   x_old;
	//ONE_DIM_ARRAY(fit_x, rvec ,nfit);
	//ONE_DIM_ARRAY(fit_ref, rvec ,nfit);
	rvec  *fit_x , *fit_ref; 
	/////////////////////////////////////
	if(ind_ref)
	{
		ALLO_ONE_DIM_ARRAY(fit_ref, rvec ,nfit);
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  fit_ref[i][m]= ref[ind_ref[i] ][m];
	}
	else
	{
		fit_ref = ref;
		//for(i=0; i< nfit; i++) for(m=0; m<3; m++)  fit_ref[i][m]= ref[i][m];
	}
	////////////////////////////////////
	if(ind_fit_b)
	{
		ALLO_ONE_DIM_ARRAY(fit_x, rvec ,nfit);
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  fit_x[i][m]= x_b [ind_fit_b[i] ][m];
	}
	else
	{
		fit_x = x_b;	
		//for(i=0; i< nfit; i++) for(m=0; m<3; m++)  fit_x[i][m]= x_b [i][m];
	}
	////////////
	/* Calculate the rotation matrix R */
	calc_rotation_R(3,nfit,fit_ref,fit_x,R);
	/*rotate X*/
	for(j=0; j<natoms_b; j++) 
	{
		for(m=0; m<DIM; m++) x_old[m]=x_b[j][m];
		for(r=0; r<DIM; r++) 
		{
			x_b[j][r]=0;
			for(c=0; c<DIM; c++)
			x_b[j][r]+=R[r][c]*x_old[c];
		}
	}

	if(ind_fit_b) free(fit_x);
	if(ind_ref)  free(fit_ref);
}

//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
real calc_rmsd( rvec *x , atom_id ind_fit[] ,rvec * xref,  atom_id ind_fit_ref[] ,int nfit)
{
	int    i,m;
	ONE_DIM_ARRAY(dx, real, (nfit*3) );
	if(ind_fit)
	{
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  dx[i*3 +m ] = x[ind_fit[i]][m];
	}
	else
	{
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  dx[i*3 + m] = x[i][m];
	}
	//////
	if(ind_fit_ref)
	{
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  dx[i*3 +m ] -= xref[ind_fit_ref[i]][m];
	}
	else
	{
		for(i=0; i< nfit; i++) for(m=0; m<3; m++)  dx[i*3 +m ] -= xref[i][m];;
	}
	////////////
	real tt=0;	
	for(i=0; i< 3*nfit ; i++) tt += dx[i] * dx[i] ;
	free(dx);
	return 10*sqrt(tt/(real)nfit);
   
}

////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
void gpa(const char 	*trxfile,   // traj file to read
			rvec  		**xav,       // the output of the function is the coor on the average 
			int 		nfit,
			atom_id 	*ifit,
			gmx_bool 	bPBC,
			gmx_rmpbc_t     gpbc,
			int 	TC,
			real 	convergence,
			int verbose,
			float first_time,
			float last_time
		      )
{
    real            inv_nframes;
    real            t;
    int             nat, nframes0;
    int  i,d;
    xtc_t f, fref;
    Read_first_xtc( trxfile,&fref, first_time );
    nat= fref.natoms;
    int ndim = nat*3;
	ONE_DIM_ARRAY(xref, rvec ,nat);
	for (i = 0; i < nat; i++) for (d = 0; d < DIM; d++) xref[i][d]=fref.x[i][d] ;
	reset_frame(&fref);
	ALLO_ONE_DIM_ARRAY((*xav), rvec ,nat);
	for (i = 0; i < nat; i++) for (d = 0; d < 3; d++) (*xav)[i][d]=0.0;
	printf("\n********************\nStart GPA Iterations\n********************");
	int ci;
	real rmsd,Trmsd=0;
	for (ci = 1; ci < TC; ci++)
	{
		printf("\n************** Iteration Nr...%d *************\n", ci);
		center_x(xref, ifit,  nfit,nat, NULL );
		Read_first_xtc( trxfile,&f, first_time );
		Trmsd=0;
		nframes0 = 0;  
    	do
    	{
    		nframes0++;
	    	if (bPBC) gmx_rmpbc(gpbc, f.natoms, f.box, f.x);
	    	center_x(f.x, ifit,  nfit,f.natoms, NULL );
        	rotate_x(nat,  xref, f.x, nfit, ifit);
			if (verbose) Trmsd += calc_rmsd(f.x, ifit, xref, ifit, nfit);
			#ifdef DOUBLE
			cblas_daxpy (ndim, 1, f.x[0], 1,(*xav)[0], 1);
			#else
			cblas_saxpy (ndim, 1, f.x[0], 1,(*xav)[0], 1);
			#endif
    	}
    	while( Read_next_xtc(&f, last_time));
    	reset_frame( &f);
		inv_nframes = 1.0/nframes0;
		#ifdef DOUBLE
		cblas_dscal(ndim, (double)(1.0/nframes0), (*xav)[0], 1);
        #else
		cblas_sscal(ndim,(float)(1.0/nframes0), (*xav)[0], 1);
		#endif
		if (verbose) printf("------------------\nAverage RMSD of the conformations = %f \n", Trmsd*inv_nframes) ;
		rmsd= calc_rmsd( (*xav), ifit, xref, ifit, nfit);		
		printf("RMSD of NEW/OLD average structure= %.6f\n", rmsd );
    	for (i = 0; i < nat; i++)
    	{ 
			for (d = 0; d < DIM; d++)
			{	
				xref[i][d] = (*xav)[i][d];
				(*xav)[i][d]=0;
			}
    	}

    	if ((ci >1) && ( rmsd < convergence) ) 
    	{
    		printf("#############################################\n");
    		printf("********************************************\n");
			printf("GPA fitting CONVERGED within %d cycles \n using an RMSD threshold of %.6f\n", ci,  convergence);
			printf("********************************************\n");
			break;
    	}
	}
	for (i = 0; i < nat; i++) for (d = 0; d < DIM; d++) (*xav)[i][d]=xref[i][d] ;
	center_x((*xav), ifit,  nfit,nat, NULL );
	free(xref);
}
//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
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
					)
{


		printf("saving the fitted traj to file %s\n",out_file );
		XDRFILE *xd= xdrfile_open(out_file,"w");
		if (NULL == xd) PRINT_ERR_OPENFILE( out_file);
		xtc_t f;
		Read_first_xtc( trxfile,&f, first_time );
		//rvec cog;
		//if (xref) get_COG(xref ,ifit,nfit,  cog );
		center_x(xref, ifit,  nfit,f.natoms, NULL );
    	do
    	{
			if (bPBC) gmx_rmpbc(gpbc, f.natoms, f.box, f.x);
        	if (xref)
			{
        		center_x(f.x, ifit,  nfit,f.natoms, NULL );
				rotate_x(f.natoms,  xref, f.x, nfit, ifit);
			}
        	write_xtc(xd,f.natoms,f.step,f.time, f.box,f.x);
    	}
    	while( Read_next_xtc(&f, last_time));
    	reset_frame( &f);
    	xdrfile_close(xd);
}

////////////////////////////////////////////////////////////////////////
////// write a frame to a one col. data file
////////////////////////////////////////////////////////////////////////
void	write_frame2datafiles ( const char *fileData,rvec *x,int nfit, atom_id 	*ifit,t_atoms *atoms)
{
	int i ,ai;
	FILE           *out;
	out = fopen(fileData, "w");
	fprintf(out, "##### atomnr_atomname coor of the average of structure of %s (x1y1z1x2y2z2...) \n", fileData);
	for (i = 0; i < nfit; i++)
	{ 
		ai=ifit[i];
		fprintf(out, "%d_%s_x  %f \n",ai+1 ,*(atoms->atomname[ai]), x[ai][0]);
		fprintf(out, "%d_%s_y  %f \n",ai+1 ,*(atoms->atomname[ai]), x[ai][1]);
		fprintf(out, "%d_%s_z  %f \n",ai+1 ,*(atoms->atomname[ai]), x[ai][2]);
	}
	fclose(out);
}

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


	
		
