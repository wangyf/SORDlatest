#!/usr/bin/env python
import sord

rundir = 'tpv3'

L = 34000, 19000, 20000
dx = 100, 100, 100
dt = dx[0] / 12500.0                           # time step size

nn = (
    int( L[0] / dx[0] + 1.5 ),
    int( L[1] / dx[1] + 1.5 ),
    int( L[2] / dx[2] + 2.5 ),
)

np3 = 8, 1, 1
#careful about nodes/processor

# Near side boundary conditions:
bc1 = 10, 10, 10

# Far side boundary conditions:
# PML absorbing boundaries for the x, y and z boundaries
bc2 = 10, 10, 10

# Material properties
fieldio = [
   ( '=', 'rho', [],     2670.0  ),        # density
   ( '=', 'vp',  [],     6000.0  ),        # P-wave speed
   ( '=', 'vs',  [],     3464.0  ),        # S-wave speed
   ( '=', 'gam', [],       0.02  ),        # high viscosity
]
hourglass = 1.0, 2.0

# Fault parameters
faultnormal = 3                             # fault plane of constant z
ihypo = (nn[0] + 1) / 2, (nn[1] + 1) / 2, nn[2]/2 + 0.5                  # hypocenter indices

j = ihypo[0] - 15000.0 / dx[0], ihypo[0] + 15000.0 / dx[0]     # X fault extent
k = ihypo[1] -  7500.0 / dx[1], ihypo[1] +  7500.0 / dx[1]     # X fault extent
l = ihypo[2]

j1_ =  ihypo[0]-1500.0 / dx[0], ihypo[0] + 1500.0/ dx[0]                # nucleation patch x extent
k1_ =  ihypo[1]-1500.0 / dx[1], ihypo[1] + 1500.0/ dx[1]                # nucleation patch y extent

fieldio += [
    ( '=', 'dc',  [],           0.4   ),    # slip weakening distance
    ( '=', 'mud', [],           0.525 ),    # dynamic friction
    ( '=', 'mus', [],             1e4 ),    # static friction - locked section
    ( '=', 'mus', [j,k,(),()],  0.677 ),    # static friction - slipping section
    ( '=', 'tn',  [],          -120e6 ),    # normal traction
    ( '=', 'ts',  [],            70e6 ),    # shear traction
    ( '=', 'ts',  [j1_,k1_,(),()], 81.6e6 ),    # shear traction - nucleation patch
]

 # Write fault plane output
fieldio += [
    ( '=w', 'su1',  [j,k,(),-1], 'su1'  ),  # final horizontal slip
    ( '=w', 'su2',  [j,k,(),-1], 'su2'  ),  # final vertical slip
    ( '=w', 'psv',  [j,k,(),-1], 'psv'  ),  # peak slip velocity
    ( '=w', 'trup', [j,k,(),-1], 'trup' ),  # rupture time
]

# Write slip, slip velocity, and shear traction time histories
p1_ = ihypo[0] + 10000.0 / dx[0], ihypo[1], l, () # mode II point indices
p2_ = ihypo[0], ihypo[1] + 5000.0 / dx[1], l, () # mode III point indices
for f in 'su1', 'su2', 'sv1', 'sv2', 'ts1', 'ts2':
    fieldio += [
        ( '=w', f, p1_, 'P1-' + f ),        # mode II point
        ( '=w', f, p2_, 'P2-' + f ),        # mode III point
    ]

# Launch SORD code
sord.run( locals() )

