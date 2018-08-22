# SORDlatest
SORD latest version maintained by Yongfei Wang

## New branch
Since Zheqiang Shi, Brittany Erickson, Qian Yao and Yongfei Wang joined Steven Day's rupture dynamics group, a new branch of SORD has been developed. Currently, Yongfei Wang is the chief maintainer who use it as a simulator of a physically plausible earthquake to study earthquake parameters and its evolution. This github repository stores this branch.

## Publications
-Shi, Z. Q., and Day, S. M. (2013), Rupture dynamics and ground motion from 3-D rough-fault simulations, J. Geophys. Res., 118(3), 1122-1141, doi:10.1002/jgrb.50094.
-Erickson, B. A., and Day, S. M. (2016), Bimaterial effects in an earthquake cycle model using rate-and-state friction, J. Geophys. Res., 121(4), 2480-2506, doi:10.1002/2015jb012470.
-Wang, Y., and Day, S. M. (2017), Seismic source spectral properties of crack-like and pulse-like modes of dynamic rupture, J. Geophys. Res., 122(8), 6657-6684, doi:10.1002/2017jb014454.

## original SORD
The Support Operator Rupture Dynamics (SORD) code simulates spontaneous rupture within a 3D isotropic viscoelastic solid. Wave
motions are computed on a logically rectangular hexahedral mesh, using the generalized finite difference method of support operators. Stiffness and viscous hourglass corrections are employed to suppress suppress zero-energy grid oscillation modes. The fault surface is modeled by coupled double nodes, where the strength of the coupling is determined by a linear slip-weakening friction law. External boundaries may be reflective or absorbing, where absorbing boundaries are handled using the method of perfectly matched layers (PML). The hexahedral mesh can accommodate non-planar ruptures and surface topography SORD simulations are configured with Python scripts. Underlying computations are coded in Fortran 95 and parallelized for multi-processor execution using Message Passing Interface (MPI). The code is portable and tested with a variety of Fortran 95 compilers, MPI implementations, and operating systems (Linux, Mac OS X, IBM AIX, SUN Solaris). The source code is stored [Ely's github](https://elygeo.net/coseis/index.html#sord)

## Publications
-Ely, G. P., Day, S. M., and Minster, J.-B. (2008), A support-operator method for viscoelastic wave modelling in 3-D heterogeneous media, Geophys. J. Int., 172(1), 331-344, doi:10.1111/j.1365-246X.2007.03633.x.
-Ely, G. P., Day, S. M., and Minster, J.-B. (2009), A support-operator method for 3-D rupture dynamics, Geophys. J. Int., 177(3), 1140-1150, doi:10.1111/j.1365-246X.2009.04117.x.
-Ely, G. P., Day, S. M., and Minster, J.-B. (2010), Dynamic Rupture Models for the Southern San Andreas Fault, Bull. Seismol. Soc. Am., 100(1), 131-150, doi:10.1785/0120090187.


## Development logs
- Aug 21, 2018 Add frequency-dependent attenuation (to be checked by AWP benchmark)
