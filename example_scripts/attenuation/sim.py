# -*- coding: utf-8 -*-
# @Author: yow004
# @Date:   2018-08-20 03:30:02
# @Last Modified by:   yow004
# @Last Modified time: 2018-08-23 15:43:24
#!/usr/bin/env python
import SORDlatest                         # import the sord module
#debug = 1
rundir = 'run'                            # directory location for output
dx = 100.0, 100.0, 100.0                  # spatial step length in x, y, and z
dt = 0.0075                               # time step length
nn = 61, 61, 61                           # number of mesh nodes in x, y, and z
nt = 60                                   # number of time steps
np3 = 2,1,1
fieldio = [                               # field variable input and output
    ( '=',  'rho', [], 2670.0 ),          # material density
    ( '=',  'vp',  [], 6000.0 ),          # material P-wave velocity
    ( '=',  'vs',  [], 3464.0 ),          # material S-wave velocity
    ( '=',  'gam', [], 0.3    ),          # material viscosity
    ( '=',  'qp', [], 100    ),          # P wave attenuation
    ( '=',  'qs', [], 50    ),		  # S wave attenuation
    ( '=w', 'v1',  [0,0,31,-1], 'vx' ),   # write X velocity slice output
    ( '=w', 'v2',  [0,0,31,-1], 'vy' ),   # write Y velocity slice output
    ( '=w', 'attr11',[0,0,31,0],'attr11'),
    ( '=w', 'attr22',[0,0,31,0],'attr22'),
    ( '=w', 'attr33',[0,0,31,0],'attr33'),
    ( '=w', 'attr23',[0,0,31,0],'attr23'),
    ( '=w', 'attr31',[0,0,31,0],'attr31'),
    ( '=w', 'attr12',[0,0,31,0],'attr12'),
]                                 
ihypo = 31.0, 31.0, 31.0                  # source location
source = 'fault'                        # source type
strike = 60
dip = 60
rake = 30
m0 = 1e5
#source1 = 1e6, 1e6, 1e6                   # source normal components
#source2 = 0.0, 0.0, 0.0                   # source shear components
timefunction = 'cos'                    # source time function
period = 6 * dt                           # source dominant period
attenuation = 'anelastic'                 # anelastic setup
fac = 1									  # transition frequency
ex = 0.0								  # power in power law attenuation
fp = 1									  # reference frequency
q0 = 150								  # reference attenuation
SORDlatest.run( locals() )                      # launch SORD job
