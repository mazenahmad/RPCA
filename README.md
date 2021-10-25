# Relative Principle Components Analysis (RPCA)

This is the source code of the computational tools to perform RPCA of the conformational changes upon biomolecular interactions.

Ahmad, M.; Helms, V.; Kalinina, O. V.; Lengauer, T. [Relative Principle Components Analysis: Application to Analyzing Biomolecular Conformational Changes.](https://doi.org/10.1021/acs.jctc.8b01074)  J. Chem. Theory Comput. 2019 , doi:10.1021/acs.jctc.8b01074.
### Installing 
#### Prerequisites:
**GROMACS-4.6.5**: The tools use few functions from Gromacs-4.6.5 and have to be linked to its library. See the [installation guide of Gromacs-4.6](https://www.gromacs.org/Documentation_of_outdated_versions/Installation_Instructions_4.6).
To install GROMACS-4.6.5 in the folder /usr/local:
```
cd /usr/local && mkdir gromacs && wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.5.tar.gz \
    && tar xfz gromacs-4.6.5.tar.gz && rm gromacs-4.6.5.tar.gz  && cd gromacs-4.6.5 && mkdir build  && cd build \
    && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs \
    && make && make install
```

**LAPACK and BLAS** 

Please make sure to modify the Makefile to include the path to the header files of Gromacs and the header file of cblas.
A simple way to install BLAS and LAPACK libraries on Ubuntu:
```
sudo apt-get update -y && sudo apt-get install libatlas-base-dev 
```
The tools g_RPCA  g_GPA g_sdiag g_rpcana will be compiled using the command
```
make all
```


### License
This code is distributed without a particular license in the hope that it will be useful. You can freely use it for your scientific research BUT WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
### Acknowledgments
* The source code of Gromacs was used as a base for a part of the source code.
* The nonlinear minimizer was written and modified after the [source code](https://github.com/SurajGupta/r-source/blob/master/src/appl/uncmin.c) of nlm () function from R, which in-turn is a C translation of the FORTRAN code of the nonlinear minimization algorithm by [Dennis and Schnabel](https://www.amazon.com/Numerical-Unconstrained-Optimization-Nonlinear-Mathematics/dp/0898713641)

## Coming soon: Documentation and a tutorial

