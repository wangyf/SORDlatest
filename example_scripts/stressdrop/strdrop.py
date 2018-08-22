#!/usr/bin/env python

"""
Simulation: 3D static stress drop
"""
import numpy, SORDlatest, sys
#itstats = 1
#debug = 3

rundir = 'strdrop'

dx = 200, 200, 200

np3 = 2,2,1

# dimensions
L = 28000., 28000., 28000

T = 25. 
#dt = dx[0]/12500.
dt = 0.02
nt = int( T / dt + 1.5 )
#nt = 4
nn = (
   int( L[0] / dx[0] + 1.5 ),
   int( L[1] / dx[1] + 1.5 ),
   int( L[2] / dx[2] + 2.5 ),
)

# material properties and initial conditions
ihypo = (nn[0] + 1) / 2+0.5, (nn[1] + 1) / 2+0.5, (nn[2] + 1) / 2
print ihypo
npml = 10
#ppml = 2
#vpml = 10
faultnormal = 3

# boundary conditions
bc1 = 10, 10, 10 
bc2 = 10, 10, 10

# material properties
hourglass = 1.0, 2.0
_lam = 3.3e10
_mu = 3.3e10
_rho = 2700.
print _rho,numpy.sqrt(_lam+2*_mu/_rho),numpy.sqrt(_mu/_rho)
fieldio = [
    ( '=', 'rho', [], _rho  ),
    ( '=', 'vp',  [], numpy.sqrt((_lam+2*_mu)/_rho)  ),
    ( '=', 'vs',  [], numpy.sqrt(_mu/_rho)  ),
#    ( '=', 'gam', [], 0.1),
] 

# initial volume stress input
#ivols = 'yes'
#eplasticity = 'plastic'

fieldio += [
	( '=', 'ts', [],   60.e6 ),
	( '=', 'tn', [],   -120.e6 ),
]

#fieldio += [
#    ( '=', 'a11', [], -100.e6 ),
#    ( '=', 'a22', [], -100.e6 ),
#    ( '=', 'a33', [], -100.e6 ),    
#    ( '=', 'a12', [],   45.e6 ),
#    ( '=', 'mco', [], 0),
#    ( '=', 'mco', [], 0e6),
#    ( '=', 'phi', [], 0.8),
#]

friction = 'frictionless'
intype = 'fsu'
period = 0.1*T
fieldio += [
	( '=r', 'su1', [(21,-1-20),(21,-1-20),(),1],   'slip.bin' ),
]

#k = 25,nn[1]-26
#fieldio += [
#    ( '=', 'dc',  [], 0.4     ),
#    ( '=', 'mus', [], 1e10     ),
#    ( '=', 'mus', [k,(),()], 0.6     ),
#    ( '=', 'mud', [], 0.4     ),
#]
#friction = 'frictionless'
#intype = 'sa'
# Read grid
fieldio += [
#   ( '=r', 'x1',  [],  'x.bin'  ),
#   ( '=r', 'x2',  [],  'y.bin'  ),
#   ( '=r', 'x3',  [],  'z.bin'  ),
]


# nucleation
#rnucl = 1000.
#delts = 0.0
#tmnucl = -1. 

fieldio += [

     ( '=w', 'tsm',  [(21,-1-20),(21,-1-20),1,1], 'ts0'),
     ( '=w', 'tnm',  [(21,-1-20),(21,-1-20),1,1], 'tn0'),
     ( '=w', 'tsm',  [(21,-1-20),(21,-1-20),1,-1], 'tse'),
     ( '=w', 'tnm',  [(21,-1-20),(21,-1-20),1,-1], 'tne'),
     ( '=w', 'sum',  [(21,-1-20),(21,-1-20),1,-1], 'sue'),
     ( '=w', 'sum',  [(21,-1-20),(21,-1-20),(),()], 'suline'),
     ( '=w', 'svm',  [(21,-1-20),(21,-1-20),(),()], 'svline'),
     ( '=w', 'sam',  [(21,-1-20),(21,-1-20),(),()], 'saline'),
     ( '=w', 'tsm',  [(21,-1-20),(21,-1-20),(),()], 'tsline'),
     

]

SORDlatest.run( locals() )

