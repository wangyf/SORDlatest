# -*- coding: utf-8 -*-
# @Author: yow004
# @Date:   2018-08-18 01:41:21
# @Last Modified by:   yow004
# @Last Modified time: 2018-08-23 15:29:21
#!/usr/bin/env python
"""
Default simulation parameters

Spatial difference operator level:

  0: Auto pick 2 or 6
  1: Mesh with constant spacing dx
  2: Rectangular mesh
  3: Parallelepiped mesh
  4: One-point quadrature
  5: Exactly integrated elements
  6: Saved operators, nearly as fast as 2, but doubles the memory usage
"""

# I/O and code execution parameters
np3 = 1, 1, 1			# number of processors in j k l
mpin = 1			# input:  0=separate files, 1=MPI-IO, -1=non-collective MPI-IO
mpout = 1			# output: 0=separate files, 1=MPI-IO, -1=non-collective MPI-IO
itstats = 10			# interval for calculating statistics
itio = 50			# interval for writing i/o buffers
itcheck = 0			# interval for check-pointing (0=off)
itstop = 0			# for testing check-pointing, simulates code crash
debug = 0			# >0 verbose, >1 sync, >2 mpi vars, >3 I/O

# Wave model parameters
nn = 41, 41, 42			# number of nodes in j, k, l (double nodes counted)
nt = 41				# number of time steps
dx = 100.0, 100.0, 100.0	# spatial step length
dt = 0.0075			# time step length
tm0 = 0.0			# initial time
affine = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0) # grid transformation
gridnoise = 0.0			# random noise added to mesh, assumes planar fault
oplevel = 0			# spatial difference operator level
vdamp = -1.0			# Vs dependent damping
hourglass = 1.0, 1.0		# hourglass stiffness (1) and viscosity (2)
fieldio = [			# field I/O
    ('=', 'rho', [], 2670.0),	# density
    ('=', 'vp',  [], 6000.0),	# P-wave speed
    ('=', 'vs',  [], 3464.0),	# S-wave speed
    ('=', 'gam', [],    0.0),	# viscosity
    ('=', 'mco', [],    0.0),   # material cohesion
]
rho1 = -1.0			# min density
rho2 = -1.0			# max density
vp1 = -1.0			# min P-wave speed
vp2 = -1.0			# max P-wave speed
vs1 = -1.0			# min S-wave speed
vs2 = -1.0			# max S-wave speed
gam1 = -1.0			# min viscosity
gam2 = 0.8			# max viscosity
npml = 10			# number of PML damping nodes
ppml = 2			# PML exponend, 1-4. Generally 2 is best.
vpml = -1.0			# PML damping velocity, <0 default to min, max V_s harmonic mean
bc1 = 0, 0, 0			# boundary condition - near side
bc2 = 0, 0, 0			# boundary condition - far side
ihypo = 0, 0, 0			# hypocenter indices (with fractional values), 0 = center
rexpand = 1.06			# grid expansion ratio
n1expand = 0, 0, 0		# number of grid expansion nodes - near side
n2expand = 0, 0, 0		# number of grid expansion nodes - far side

# Dynamic rupture parameters
faultnormal = 0			# normal direction to fault plane (0=no fault)
faultopening = 0		# 0=not allowed, 1=allowed
slipvector = 1.0, 0.0, 0.0	# shear traction direction for ts1
svtol = 0.001			# slip velocity considered rupturing

# nucleation by gradually overstressing
vrup = -1.0                     # nucleation rupture velocity, negative = no nucleation
rcrit = 0.0                     # nucleation critical radius
trelax = -1.0                   # nucleation relaxation time (0.07sec)
rrelax = 100					# prescribed cohesive zone size (100 m by default)
tslope = -1.0                   # nucleation relaxation slope (choose trelax or tslope) (5MPa/km)
rnucl = 0.0                     # nucleation radius/half-width for Ts perturbation
tmnucl = 1.0                    # duration time for Ts pertubation. if negative, perturb at t=0
delts = 0.0                     # peak amplitude increase for Ts perturbation
psidelts = -1.0                 # peak amplitude reduction ratio for psi perturbation
nstage = 1						# used in forced rupture. Numbers of stages of different rupture velocity
vrupstage = 3300				# rupture velocity in each segment vrup = 3000, 3500, 5000
sizestage = 2000				# each segment size

# friction type
friction = 'slipweakening'      # or 'rateandstate' 'frictionless (need input sa)' 'forced (need rrelax)'
intype = 'sa'			  # or 'sv' or 'sa' this defines kinematic input type slip rate or slip acc (time history)
					  # fsu is input static final slip to calculate static stress change (tiny density)
					  # period is used for artificial period static slip input

# Prakash-Clifton gradual response to change of normal traction
pcdep = 'no'                    # or 'yes'  

# elastic or plastic
eplasticity = 'elastic'         # or 'elastoplastic' 
tv = -1.0                       # viscoplastic relaxation time Tv = dx/Vs = 50./3464.

# frequency-dependent attenuation
attenuation = 'elastic'         # or 'anelastic'
fac = 1.0			# factor to shift bandwith for constant Q only (f0 in paper)
q0 = 150.0			# Q0 value in Q0*f^n relation (not directly used here)
ex = 0.0 			# Exp of powerlaw (use 0.1:0.1:0.9 for Q(f),0.0 for constant
fp = 1.0			# reference frequency (velocity of input media)

# Skempton's coefficient in pore elastic medium
#skepb = 0.0                     # skepb in (0,1) 0.0 means no pore pressure change, 
                                # 1.0 means effective normal traction keeps constant

# initial volume stress input 
ivols = 'no'

# Finite source parameters
source = 'none'		# 'moment', 'potency', 'force','fault', or 'none'
nsource = 0			# number of sub-faults

# Point source parameters
source1 = 0.0, 0.0, 0.0		# normal components
source2 = 0.0, 0.0, 0.0		# shear components
strike = 90.
dip = 90.
rake = 0.
m0 = 1.e16                # N.m (Mw4.6)
timefunction = 'none'		# time function, see util.f90 for details.
period = 10 * dt		# dominant period

# Placeholders
i1pml = None
i2pml = None

