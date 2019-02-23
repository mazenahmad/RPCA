
#ifndef UNCMIN_H
#define UNCMIN_H


/* type of pointer to the target and gradient functions */
typedef void (*fcn_p)(int, double *, double *, void *);

/* type of pointer to the hessian functions */
typedef void (*d2fcn_p)(int, int, double *, double *, void *);


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////		


void nlmin( fcn_p fcn , int n , double *x, void *state, double *fpls ,int itnlim , int method,int print_code );

void nlmin_slim( fcn_p fcn , int n , double *x, void *state, double *fpls ,int itnlim ,  int verbose );
		
		
		
#endif
