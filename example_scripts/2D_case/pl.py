#!/usr/bin/env python
"""
2-D Simulation: straight fault, SW friction
"""
import  SORDlatest

rundir = 'pl'

dx = 100.0, 100.0, 100.0
#dx = 50.0, 50.0, 50.0

np3 = 1, 1, 1

# dimensions
L = 40000.0, 30000.0

T = 5.0 
dt = dx[0] / 12500.0
nt = int( T / dt + 1.5 )
#nt = 10

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

# material properties and initial conditions
ihypo = (nn[0] + 1) / 2, (nn[1] + 1) / 2, 1.5

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
#    ( '=', 'a33', [], -100.e6 ),
    ( '=', 'a12', [],   45.e6 ),
]

# elasticity or elastoplasticity
eplasticity = 'plastic'
fieldio += [
    ( '=', 'phi',  [],  0.7  ),
]

# friction type and parameters
friction = 'slipweakening'
fieldio += [
    ( '=', 'dc',  [], 0.4     ),
    ( '=', 'mus', [], 0.6     ),
    ( '=', 'mud', [], 0.4     ),
]

fieldio += [
#    ( '=', 'gam', [], 0.02    ),
]


fieldio += [
    ( '=w', 'x1', [(), (), 1, ()], 'xx' ),
    ( '=w', 'x2', [(), (), 1, ()], 'yy' ),
]

# nucleation
delts = 0.38 

# output 
_dnt1 = 1 # int( 0.016 / dt )
_dnt2 = 5 # int( 0.08 / dt )
ihypo = (nn[0] + 1) / 2, (nn[1] + 1) / 2, 1.5
k = ihypo[1]

fieldio += [
    ( '=w', 'sl',  [(),k,1,(1,-1,_dnt1)],  'sl'  ),
    ( '=w', 'tsm', [(),k,1,(1,-1,_dnt1)],  'tsm' ),
    ( '=w', 'tnm', [(),k,1,(1,-1,_dnt1)],  'tnm' ),
    ( '=w', 'svm', [(),k,1,(1,-1,_dnt1)],  'svm' ),
    ( '=w', 'w22', [(),(),1,(1,-1,_dnt2)], 's22' ),
    ( '=w', 'w12', [(),(),1,(1,-1,_dnt2)], 's12' ),
    ( '=w', 'w11', [(),(),1,(1,-1,_dnt2)], 's11' ),
    ( '=w', 'epm', [(),(),1,(1,-1,_dnt2)], 'epm' ),
    ( '=w', 'a1',  [(),(),1,(1,-1,_dnt2)], 'ax'  ),
    ( '=w', 'a2',  [(),(),1,(1,-1,_dnt2)], 'ay'  ),
    ( '=w', 'v1',  [(),(),1,(1,-1,_dnt2)], 'vx'  ),
    ( '=w', 'v2',  [(),(),1,(1,-1,_dnt2)], 'vy'  ), 
]

SORDlatest.run( locals() )
