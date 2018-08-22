#!/usr/bin/env python
"""
2-D Simulation: rough fault, RSD friction
"""
import numpy, sord, sys

rundir = 'elastic'

dx = 50.0, 50.0, 50.0

np3 = 64, 64, 1

# dimensions
L = 60000.0, 30000.0

T = 8 
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
    ( '=', 'a12', [],   45.e6 ),
]

# elasticity or elastoplasticity
#eplasticity = 'plastic'
#fieldio += [
#    ( '=', 'phi',  [],  0.75  ),
#]
8
# friction type and parameters
friction = 'slipweakening'
fieldio += [
    ( '=', 'dc',  [], 0.04     ),
    ( '=', 'mus', [], 0.6     ),
    ( '=', 'mud', [], 0.4     ),
    ( '=', 'gam', [], 0.02    ),
]


# friciton type and parameters
#friction = 'rateandstate'
#k1_ = 1, nn[1] / 2
#k2_ = nn[1] / 2 + 1, nn[1]
#fieldio += [
#    ( '=', 'af',  [],   0.01  ),
#    ( '=', 'bf',  [],  0.014  ),
#    ( '=', 'vw',  [],    0.1  ),
#    ( '=', 'v0',  [],  1.e-6  ),
#    ( '=', 'f0',  [],    0.6  ),
#    ( '=', 'll',  [],    0.4  ),
#    ( '=', 'fw',  [],    0.2  ),
#    ( '=', 'v1',  [(), k1_, (), 1],  0.5e-16 ),
#    ( '=', 'v1',  [(), k2_, (), 1], -0.5e-16 ),
#]

# read grid file 
#fieldio += [
#    ( '=R', 'x2', [(), (), 1, ()], 'y_coarse'),
#]

fieldio += [
    ( '=w', 'x1', [(), (), 1, ()], 'xx' ),
    ( '=w', 'x2', [(), (), 1, ()], 'yy' ),
]

# nucleation
delts = .38
rnucl = 1000

# output 
_dnt = int( 0.04 / dt )
_ddnt = int( 0.4 / dt)
ihypo = 30000./dx[0], (nn[1] + 1) / 2, 1.5
k = ihypo[1]

fieldio += [
    ( '=w', 'sl',  [(),k,1,(1,-1,_dnt)], 'sl'   ),
    ( '=w', 'tsm', [(),k,1,(1,-1,_dnt)], 'tsm'  ),
    ( '=w', 'svm', [(),k,1,(1,-1,_dnt)], 'svm'  ),
    ( '=w', 'tnm', [(),k,1,(1,-1,_dnt)], 'tnm'  ),
    ( '=w', 'psi', [(),k,1,(1,-1,_dnt)], 'psi'  ),
    ( '=w', 'a1', [(),(),1,(1,-1,_ddnt)], 'ax'   ),
    ( '=w', 'a2', [(),(),1,(1,-1,_ddnt)], 'ay'   ),
    ( '=w', 'epm', [(),(),1,(1,-1,_ddnt)], 'epm' ),
    ( '=w', 'v1', [(),(),1,(1,-1,_ddnt)], 'vx'   ),
    ( '=w', 'v2', [(),(),1,(1,-1,_ddnt)], 'vy'   ),
    ( '=w', 'u1', [(),(),1,(1,-1,_ddnt)], 'ux'   ),
    ( '=w', 'u2', [(),(),1,(1,-1,_ddnt)], 'uy'   ),
    ( '=w', 'w11', [(),(),1,(1,-1,_ddnt)], 's11'   ),
    ( '=w', 'w12', [(),(),1,(1,-1,_ddnt)], 's12'   ),
    ( '=w', 'w22', [(),(),1,(1,-1,_ddnt)], 's22'   ),
]



# Write slip, slip velocity, acceleration and shear traction time histories
p1_ = ihypo[0] + 0.0 / dx[0], ihypo[1], ihypo[2], () # mode II point indices
p2_ = ihypo[0] + 10000.0 / dx[0], ihypo[1], ihypo[2], () # mode II point indices
p3_ = ihypo[0] - 10000.0 / dx[0], ihypo[1], ihypo[2], () # mode II point indices

p4_ = ihypo[0], ihypo[1] + 3000.0 / dx[1], ihypo[2], () # mode III point indices
p5_ = ihypo[0], ihypo[1] - 3000.0 / dx[1], ihypo[2], () # mode III point indices
p6_ = ihypo[0] + 10000.0 / dx[0], ihypo[1] + 3000.0 / dx[1], ihypo[2], () # mode III point indices
p7_ = ihypo[0] + 10000.0 / dx[0], ihypo[1] - 3000.0 / dx[1], ihypo[2], () # mode III point indices
p8_ = ihypo[0] - 10000.0 / dx[0], ihypo[1] + 3000.0 / dx[1], ihypo[2], () # mode III point indices
p9_ = ihypo[0] - 10000.0 / dx[0], ihypo[1] - 3000.0 / dx[1], ihypo[2], () # mode III point indices


for f in 'u1', 'u2', 'v1', 'v2', 't1', 't2', 'a1', 'a2':
    fieldio += [
        ( '=w', f, p1_, 'P1-' + f ),        # mode II point
        ( '=w', f, p2_, 'P2-' + f ),        # mode II point
        ( '=w', f, p3_, 'P3-' + f ),        # mode II point
        ( '=w', f, p4_, 'P4-' + f ),        # mode III/Users/britte/Desktop/short_course_ex.pdf  point
        ( '=w', f, p5_, 'P5-' + f ),        # mode III point
        ( '=w', f, p6_, 'P6-' + f ),        # mode III point
        ( '=w', f, p7_, 'P7-' + f ),        # mode III point
        ( '=w', f, p8_, 'P8-' + f ),        # mode III point
        ( '=w', f, p9_, 'P9-' + f ),        # mode III point
    ]




sord.run( locals() )

