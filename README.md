# SORDlatest
## SORD latest version maintained by Yongfei Wang

# SORD
## The Support Operator Rupture Dynamics (SORD) code simulates spontaneous rupture within a 3D isotropic viscoelastic solid. Wave
motions are computed on a logically rectangular hexahedral mesh, using the generalized finite difference method of support operators. Stiffness and viscous hourglass corrections are employed to suppress suppress zero-energy grid oscillation modes. The fault surface is modeled by coupled double nodes, where the strength of the coupling is determined by a linear slip-weakening friction law. External boundaries may be reflective or absorbing, where absorbing boundaries are handled using the method of perfectly matched layers (PML). The hexahedral mesh can accommodate non-planar ruptures and surface topography SORD simulations are configured with Python scripts. Underlying
computations are coded in Fortran 95 and parallelized for multi-processor execution using Message Passing Interface (MPI). The code is portable and tested with a variety of Fortran 95 compilers, MPI implementations, and operating systems (Linux, Mac OS X, IBM AIX, SUN Solaris). 

BACKGROUND
The formulation, numerical algorithm, and verification of the SORD method are described by Ely, Day, and Minster (2008) for wave propagation, and Ely, Day, and Minster (2009) for spontaneous rupture. Ely, Day, and Minster (2010) present an application to simulating earthquakes in southern California.



- Aug 21, 2018 Add frequency-dependent attenuation (to be checked by AWP benchmark)
