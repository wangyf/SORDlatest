#!/usr/bin/env python
"""
2-D Simulation: rough fault, RSD friction
"""
import numpy, sord, sys

rundir = 'run'

dx = 100.0, 100.0, 100.0

np3 = 8, 1, 1

# dimensions
L = 60000.0, 30000.0

T = 200 
dt = dx[0] / 12500.0
nt = int( T / dt + 1.5 )

nn = (
    int( L[0] / dx[0] + 1.5 ),
    int( L[1] / dx[1] + 2.5 ),
    2,
)

npml = 10
faultnormal = 2

# boundary conditions
bc1 = 10, 10, 1
bc2 = 10, 10, 1

# material properties
hourglass = 1.0, 2.0
fieldio = [
    ( '=', 'rho', [], 2670.0  ),
    ( '=', 'vp',  [], 6000.0  ),
    ( '=', 'vs',  [], 3464.0  ),
]  

# initial volume stress input 
ivols = 'yes'
fieldio += [
    ( '=', 'a11', [], -100.e6 ),
    ( '=', 'a22', [], -100.e6 ),
    ( '=', 'a33', [], -100.e6 ),
    ( '=', 'a12', [],   30.e6 ),
]

# elasticity or elastoplasticity
eplasticity = 'plastic'
fieldio += [
    ( '=', 'phi',  [],  0.75  ),
]

# friciton type and parameters
friction = 'rateandstate'
k1_ = 1, nn[1] / 2
k2_ = nn[1] / 2 + 1, nn[1]
fieldio += [
    ( '=', 'af',  [],   0.01  ),
    ( '=', 'bf',  [],  0.014  ),
    ( '=', 'vw',  [],    0.1  ),
    ( '=', 'v0',  [],  1.e-6  ),
    ( '=', 'f0',  [],    0.6  ),
    ( '=', 'll',  [],    0.4  ),
    ( '=', 'fw',  [],    0.2  ),
    ( '=', 'v1',  [(), k1_, (), 1],  0.5e-16 ),
    ( '=', 'v1',  [(), k2_, (), 1], -0.5e-16 ),
]

# read grid file 
# fieldio += [
#    ( '=R', 'x2', [(), (), 1, ()], 'y_coarse'),
#]

fieldio += [
    ( '=w', 'x1', [(), (), 1, ()], 'xx' ),
    ( '=w', 'x2', [(), (), 1, ()], 'yy' ),
]

# nucleation
delts = 2
rnucl = 2000

# output 
_dnt = int( 0.04 / dt )
ihypo = 30000./dx[0], (nn[1] + 1) / 2, 1.5
k = ihypo[1]

fieldio += [
    ( '=w', 'sl',  [(),k,1,(1,-1,_dnt)], 'sl'   ),
    ( '=w', 'tsm', [(),k,1,(1,-1,_dnt)], 'tsm'  ),
    ( '=w', 'svm', [(),k,1,(1,-1,_dnt)], 'svm'  ),
    ( '=w', 'tnm', [(),k,1,(1,-1,_dnt)], 'tnm'  ),
    ( '=w', 'psi', [(),k,1,(1,-1,_dnt)], 'psi'  ),
    ( '=w', 'a1', [(),(),1,(1,-1,_dnt)], 'ax'   ),
    ( '=w', 'a2', [(),(),1,(1,-1,_dnt)], 'ay'   ),
    ( '=w', 'epm', [(),(),1,(1,-1,_dnt)], 'epm' ),
]

sord.run( locals() )

