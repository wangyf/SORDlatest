# SORDlatest
SORD latest version maintained by Yongfei Wang

## New branch
Since Zheqiang Shi, Brittany Erickson, Qian Yao and Yongfei Wang joined Steven Day's rupture dynamics group, a new branch of SORD has been developed. Currently, Yongfei Wang is the chief maintainer who use it as a simulator of a physically plausible earthquake to study earthquake parameters and its evolution. This github repository stores this branch.

## Publications
1. [Shi, Z. Q., and Day, S. M. (2013), Rupture dynamics and ground motion from 3-D rough-fault simulations, J. Geophys. Res., 118(3), 1122-1141](doi.org/10.1002/jgrb.50094)
2. [Wang, Y., and Day, S. M. (2017), Seismic source spectral properties of crack-like and pulse-like modes of dynamic rupture, J. Geophys. Res., 122(8), 6657-6684](doi.org/10.1002/2017jb014454)
3. Wang, Y., Day, S. M. and Denolle, M. A. (2018), Geometric Controls on Pulse-like Rupture in a Dynamic Model of the 2015 Gorkha Earthquake, J. Geophys. Res., *submitted*


## Original SORD
The Support Operator Rupture Dynamics (SORD) code simulates spontaneous rupture within a 3D isotropic viscoelastic solid. Wave
motions are computed on a logically rectangular hexahedral mesh, using the generalized finite difference method of support operators. Stiffness and viscous hourglass corrections are employed to suppress suppress zero-energy grid oscillation modes. The fault surface is modeled by coupled double nodes, where the strength of the coupling is determined by a linear slip-weakening friction law. External boundaries may be reflective or absorbing, where absorbing boundaries are handled using the method of perfectly matched layers (PML). The hexahedral mesh can accommodate non-planar ruptures and surface topography SORD simulations are configured with Python scripts. Underlying computations are coded in Fortran 95 and parallelized for multi-processor execution using Message Passing Interface (MPI). The code is portable and tested with a variety of Fortran 95 compilers, MPI implementations, and operating systems (Linux, Mac OS X, IBM AIX, SUN Solaris). The source code is stored in [Ely's github](https://elygeo.net/coseis/index.html#sord)

## Publications
1. [Ely, G. P., Day, S. M., and Minster, J.-B. (2008), A support-operator method for viscoelastic wave modelling in 3-D heterogeneous media, Geophys. J. Int., 172(1), 331-344](doi.org/10.1111/j.1365-246X.2007.03633.x)
2. [Ely, G. P., Day, S. M., and Minster, J.-B. (2009), A support-operator method for 3-D rupture dynamics, Geophys. J. Int., 177(3), 1140-1150](doi.org/10.1111/j.1365-246X.2009.04117.x)
3. [Ely, G. P., Day, S. M., and Minster, J.-B. (2010), Dynamic Rupture Models for the Southern San Andreas Fault, Bull. Seismol. Soc. Am., 100(1), 131-150](doi.org/10.1785/0120090187)

## Publications in which SORD is used
- [Harris, R. A., et al. (2009), The SCEC/USGS dynamic earthquake rupture code verification exercise, Seismol. Res. Lett., 80(1), 119-126.](doi.org/10.1785/gssrl.80.1.119)
- [Ben-Zion, Y., Rockwell, T. K., Shi, Z. Q., and Xu, S. Q. (2012), Reversed-polarity secondary deformation structures near fault stepovers, J. Appl. Mech., 79(3), 031025](doi.org/10.1115/1.4006154)
- [Song, S. G., Dalguer, L. A., and Mai, P. M. (2013), Pseudo-dynamic source modelling with 1-point and 2-point statistics of earthquake source parameters, Geophys. J. Int., 196(3), 1770-1786](doi.org/10.1093/gji/ggt479)
- [Fan, W. Y., Shearer, P. M., and Gerstoft, P. (2014), Kinematic earthquake rupture inversion in the frequency domain, Geophys. J. Int., 199(2), 1138-1160](doi.org/10.1093/gji/ggu319)
- [Baumann, C., and Dalguer, L. A. (2014), Evaluating the compatibility of dynamic rupture-based synthetic ground motion with empirical ground-motion prediction equation, Bull. Seismol. Soc. Am., 104(2), 634-652](doi.org/10.1785/0120130077)
- [Song, S. G. (2015), The effect of fracture energy on earthquake source correlation statistics, Bull. Seismol. Soc. Am., 105(2a), 1042-1048](doi.org/10.1785/0120140207)
- [Vyas, J. C., Mai, P. M., and Galis, M. (2016), Distance and azimuthal dependence of ground-motion variability for unilateral strike-slip ruptures, Bull. Seismol. Soc. Am., 106(4), 1584-1599](doi.org/10.1785/0120150298)
- [Mai, P., Galis, M., Thingbaijam, K., Vyas, J., and Dunham, E. (2017), Accounting for Fault Roughness in Pseudo-Dynamic Ground-Motion Simulations, Pure. Appl. Geophys., 174(9), 3419-3450](doi.org/10.1007/s00024-017-1536-8)
- [Song, S. G., and Dalguer, L. A. (2017), Synthetic Source Inversion Tests with the Full Complexity of Earthquake Source Processes, Including Both Supershear Rupture and Slip Reactivation, Pure. Appl. Geophys., 174(9), 3393-3418](doi.org/10.1007/s00024-017-1514-1)
- [Passone, L., and Mai, P. M. (2017), Kinematic Earthquake Ground-Motion Simulations on Listric Normal Faults, Bull. Seismol. Soc. Am., 107(6), 2980-2993](doi.org/10.1785/0120170111)
- [Vyas, J. C., Mai, P. M., Galis, M., Dunham, E. M., and Imperatori, W. (2018), Mach wave properties in the presence of source and medium heterogeneity, Geophys. J. Int., 214(3), 2035-2052](doi.org/10.1093/gji/ggy219)
- [Harris, R. A., et al. (2018), A Suite of Exercises for Verifying Dynamic Earthquake Rupture Codes, Seismol. Res. Lett., 89(3), 1146-1162](doi.org/10.1785/0220170222)

## Development logs
- Aug 21, 2018 Add frequency-dependent attenuation (to be checked by AWP benchmark)
- Aug 23, 2018 Add a function to exert a point source by setting strike, dip, rake and moment.

## Bug found and fixed
- Bug in computing radiated energy, negative value

*If you found a bug in this program, you are so welcome to report it.*
