# Relative Principle Components Analysis (RPCA)

This is the source code of the computational tools to perform RPCA of the conformational changes upon biomolecular interactions.

Ahmad, M.; Helms, V.; Kalinina, O. V.; Lengauer, T. [Relative Principle Components Analysis: Application to Analyzing Biomolecular Conformational Changes.](https://doi.org/10.1021/acs.jctc.8b01074)  J. Chem. Theory Comput. 2019 , doi:10.1021/acs.jctc.8b01074.
### Installing 
#### Prerequisites:
**GROMACS-4.6.5**: The tools use few functions from Gromacs-4.6.5 and have to be linked to its library. See the [installation guide of Gromacs-4.6](http://www.gromacs.org/Documentation/Installation_Instructions_4.6).

**LAPACK and BLAS** 

The tools g_RPCA  g_GPA g_sdiag g_rpcana will be compiled using the command
```
make all
```
Please make sure to modify the Makefile to include the path to the header files of Gromacs and the header file of cblas.

### License
This code is distributed without a particular license. You can freely use it for your scientific research. 
      
### Acknowledgments
* The source code of Gromacs was used as a base for a part of the source code.
* The nonlinear minimizer was written and modified after the [source code](https://github.com/SurajGupta/r-source/blob/master/src/appl/uncmin.c) of nlm () function from R, which in-turn is a C translation of the FORTRAN code of the nonlinear minimization algorithm by [Dennis and Schnabel](https://www.amazon.com/Numerical-Unconstrained-Optimization-Nonlinear-Mathematics/dp/0898713641)

## Coming soon: official distribution of RPCA independent from Gromacs with documentation and a tutorial

