# -*- coding: utf-8 -*-
# @Author: yow004
# @Date:   2018-08-23 19:24:14
# @Last Modified by:   yow004
# @Last Modified time: 2018-08-25 22:46:37
#!/usr/bin/env python

"""
Simulation: 2D In-plane rupture 
"""
import numpy, SORDlatest, sys
#itstats = 1
#debug = 3

rundir = 'debug'

dx = 100, 100, 400

np3 = 1,1,1

# dimensions
L = 4000., 4000.

T = 10. 
dt = dx[0] / 12500.0
nt = int( T / dt + 1.5 )
nt = 401
nn = (
   int( L[0] / dx[0] + 1.5 ),
   int( L[1] / dx[1] + 1.5 ),
   2,
)

# material properties and initial conditions
ihypo = (nn[0] + 1) / 2, (nn[1] + 1) / 2, 1.5

npml = 2
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
    ( '=', 'a12', [],   45.e6 ),
]

friction = 'slipweakening'
fieldio += [
    ( '=', 'dc',  [], 0.4     ),
    ( '=', 'mus', [], 0.6     ),
    ( '=', 'mud', [], 1e6     ),
]
#friction = 'frictionless'

# Read grid
fieldio += [
#   ( '=r', 'x1',  [],  'x.bin'  ),
#   ( '=r', 'x2',  [],  'y.bin'  ),
#   ( '=r', 'x3',  [],  'z.bin'  ),
]


# nucleation
rnucl = 100.
delts = 0.38 
tmnucl = -1. 

k = ihypo[1]
fieldio += [
#     ( '=w', 'a1',  [(),(k,k+1),1,(1,-1)], 'a1'),
#     ( '=w', 'v1',  [(),(k,k+1),1,(1,-1)], 'v1'),
#     ( '=w', 'u1',  [(),(k,k+1),1,(1,-1)], 'u1'),
#     ( '=w', 'a2',  [(),(k,k+1),1,(1,-1)], 'a2'),
#     ( '=w', 'v2',  [(),(k,k+1),1,(1,-1)], 'v2'),
#     ( '=w', 'u2',  [(),(k,k+1),1,(1,-1)], 'u2'),
     
#     ( '=w', 'a3',  [(),(k,k+1),1,(1,-1)], 'a3'),
#     ( '=w', 'v3',  [(),(k,k+1),1,(1,-1)], 'v3'),
#     ( '=w', 'u3',  [(),(k,k+1),1,(1,-1)], 'u3'), 
#     ( '=r', 'sv1',  [], 'sv1'),
#     ( '=r', 'sv2',  [], 'sv2'),         
#     ( '=w', 'sv1',  [], 'sv1'),
#     ( '=w', 'sa1',  [], 'sa1'),
#     ( '=w', 'sa2',  [], 'sa2'),
#     ( '=w', 'sv2',  [], 'sv2'),
#     ( '=w', 'ts1',  [], 'ts1'),
#     ( '=w', 'tnm',  [], 'tnm'),
     
     
#    ( '=w', 'sv3',  [], 'slipr_z'  ),
#    ( '=w', 'su2',  [(),(),801,-1], 'slip_y'),
#    ( '=w', 'su3',  [(),(),801,-1], 'slip_z'),
#    ( '=w', 'x1',  [], 'xx'  ),
#    ( '=w', 'x2',  [], 'yy'  ),
#    ( '=w', 'x3',  [], 'zz'  ), # final horizontal slip
#    ( '=w', 'mus', [], 'mustatic'),
#    ( '=w', 'nhat1',[ihypo[0],(),1,1], 'n_x'),
#    ( '=w', 'nhat2',[ihypo[0],(),1,1], 'n_y'),
#    ( '=w', 'nhat3',[ihypo[0],(),1,1], 'n_z'),
#    ( '=w', 'a2',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'a2'),
#    ( '=w', 'a3',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'a3'),
#    ( '=w', 'v2',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'v2'),
#    ( '=w', 'v3',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'v3'),
#    ( '=w', 'u2',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'u2'),
#    ( '=w', 'u3',  [ihypo[0],(),(ihypo[2],ihypo[2]+1),(1,-1,50)],'u3'),
#    ( '=w', 'af',  [], 'af'),
#    ( '=w', 'vw',  [], 'vw'),
#    ( '=w', 'ts1',  [], 'ts1'),
#    ( '=w', 'ts2',  [], 'ts2'),
#    ( '=w', 'ts3',  [], 'ts3'),
#    ( '=w', 'su1',  [], 'su1'),
#    ( '=w', 'su2',  [], 'su2'),
#    ( '=w', 'su3',  [], 'su3'),
#    ( '=w', 'sv1',  [], 'sv1'),
#    ( '=w', 'sv2',  [], 'sv2'),
#    ( '=w', 'sv3',  [], 'sv3'),

]

SORDlatest.run( locals() )

