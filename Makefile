#modify for the correct path to gromacs headers and library
GROMACS_FOLDER = /usr/local/gromacs
INCLUDE      = $(GROMACS_FOLDER)/include/gromacs/
# path to gromacs library 
GMXLIB       = $(GROMACS_FOLDER)/lib 
### modify for the correct libraries of Lapack abd Blas
LIBS         = -lgmx -lm -ldl -llapack  -lblas -lcblas
LOCAL        = ./
LOCAL_INC    =  ./include
 
LDFLAGS      = -L$(GMXLIB) 
CFLAGS	     = -I$(LOCAL) -I$(LOCAL_INC)  -I$(INCLUDE) 
 
CC           = gcc
LD           = $(CC)
 
 
g_RPCA:	src/gmx_RPCA.c
	 $(LD) -fopenmp $(CFLAGS) $(LDFLAGS) src/gmx_RPCA.c src/fitting.c src/covariance.c src/linear_algebra.c  src/strTools.c src/slim_uncmin.c src/CWfitting.c src/xdrfile.c  src/xdrfile_xtc.c -o $@ $(LIBS)
g_sdiag: src/SIMDAIG.c
	 $(LD) -fopenmp $(CFLAGS) $(LDFLAGS) src/SIMDAIG.c  src/covariance.c src/linear_algebra.c src/xdrfile.c src/xdrfile_xtc.c src/fitting.c src/strTools.c -o $@ $(LIBS)	
g_rpcana: src/RPCA_ANA.c
	$(LD) -fopenmp $(CFLAGS) $(LDFLAGS) src/RPCA_ANA.c src/fitting.c src/covariance.c src/linear_algebra.c  src/strTools.c  src/xdrfile.c  src/xdrfile_xtc.c -o $@ $(LIBS)
 
g_GPA:  src/gmx_GPA.c
	$(LD) -fopenmp $(CFLAGS) $(LDFLAGS) src/gmx_GPA.c src/fitting.c src/covariance.c src/linear_algebra.c   src/strTools.c   src/xdrfile.c src/xdrfile_xtc.c -o $@ $(LIBS) 
 
 
all:   g_RPCA  g_GPA g_sdiag g_rpcana


