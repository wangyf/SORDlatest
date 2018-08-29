#!/usr/bin/env python

"""
Source Spectral Study (Asperity model)
"""
import numpy, SORDlatest, sys
#debug = 3
#itstats=1
rundir = 'mdl3'

dx = 25, 25, 25

np3 = 8,8,16

# dimensions
L = 21000, 21000., 21000.

T = 10. 
dt = dx[0] / 12500.0 
nt = int( T / dt + 1.5 )
nn = (
   int( L[0] / dx[0] + 1.5 ),
   int( L[1] / dx[1] + 1.5 ),
   int( L[2] / dx[2] + 2.5 ),
)

ihypo = int(L[0]/2./dx[0]+1.5), int(L[1]/2./dx[1]+1.5), int(L[2]/2./dx[2]+1.5)
#ihypo = 181, 106, 76.5
#ihypo = (383-1)*4+1,(109-1)*4+1, 145
#ihypo = 393,118, 38
print ihypo
npml = 10
faultnormal = 3

# boundary conditions
bc1 = 10, 10, 10 
bc2 = 10, 10, 10

# material properties
hourglass = 1.0, 2.0
#eplasticity = 'plastic'
_beta=3464
fieldio = [
   ( '=', 'rho', [], 2670.0  ),
   ( '=', 'vp',  [], 6000.0  ),
   ( '=', 'vs',  [], _beta  ),
   ( '=', 'gam', [], 0.1    ),
#   ( '=', 'mco', [], 5e6),
#   ( '=', 'phi', [], 0.75),
]  

# initial volume stress input
#ivols = 'yes'
#_scale = 1.5
fieldio += [
   ( '=', 'tn',  [], -10e6 ),
#   ( '=r','ts',  [], 'stress.bin'),
   ( '=', 'ts', [],   4e6  ),
]  
slipvector = (1.0, 0.0, 0.0)

# forced rupture
#friction = 'slipweakening'
friction = 'forced'
#nstage = 10
#vrupstage = 0.9*_beta,0.4*_beta,0.9*_beta,0.4*_beta,0.9*_beta,0.4*_beta,0.9*_beta,0.4*_beta,0.9*_beta,0.4*_beta
#sizestage = 2000,2000,2000,2000,2000,2000,2000,2000,2000,2000
nstage = 1
vrupstage = 0.9*_beta
sizestage = 10000
rrelax = 125

fieldio += [
    ( '=', 'dc',  [],  1e10),
    ( '=', 'mus', [],  1e10),
    ( '=r', 'mud', [],  'mud.bin'),

]

# Read grid
fieldio += [
#   ( '=r', 'x1',  [],  'x.bin'  ),
#   ( '=r', 'x2',  [],  'y.bin'  ),
#   ( '=r', 'x3',  [],  'z.bin'  ),
]

svtol = 0.001
# nucleation
rcrit = 10000
#trelax = 0.08
vrup = 0.9*_beta



fieldio += [
    ( '=w', 'trup', [(),(),1,-1], 'trup'),
    ( '=w', 'tarr', [(),(),1,-1], 'tarr'),
    ( '=w', 'sum', [(),(),1,-1], 'slip'),
    ( '=w', 'tsm', [(),(),1,1],'ts0'),
    ( '=w', 'tsm', [(),(),1,-1],'tse'),
    ( '=w', 'x1',  [(),(),ihypo[2],1],'faultx'),
    ( '=w', 'x2',  [(),(),ihypo[2],1],'faulty'),
    ( '=w', 'x3',  [(),(),ihypo[2],1],'faultz'),

#    ( '=w', 'mr11',[], 'mr11'),
#    ( '=w', 'mr22',[], 'mr22'),
#    ( '=w', 'mr33',[], 'mr33'),
#    ( '=w', 'mr12',[], 'mr12'),
    ( '=w', 'mr31',[], 'mr31'),
#    ( '=w', 'mr23',[], 'mr23'),

    ( '=w', 'svm', [(),(),1,(1,-1,25)], 'sliprate'),
    ( '=w', 'tsm', [(),(),1,(1,-1,25)], 'shearstress'),
]

SORDlatest.run( locals() )

