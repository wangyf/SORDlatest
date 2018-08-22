notes = """
IGPP laptop

ssh -fNL localhost:8022:pisco.sdsu.edu:22 sciences.sdsu.edu
ssh -p 8022 localhost

Use MPICH instead of OpenMPI:
export PATH="/opt/mpich2/gnu/bin:${PATH}"
"""
login = 'igpp13grad9-2.local'
hosts = 'igpp13grad9-2.local',
maxnodes = 1
maxcores = 2
maxram = 8000000
fortran_serial = 'gfortran',
fortran_mpi = 'mpif90',
#_ = '-fimplicit-none', '-Wall', '-Wtabs', '-Wl,-no_compact_unwind', '-o'
_ = '-fimplicit-none', '-o' 
fortran_flags = {
    'g': ('-fbounds-check', '-ffpe-trap=invalid,zero,overflow', '-g') + _,
    't': ('-fbounds-check', '-ffpe-trap=invalid,zero,overflow') + _,
    'p': ('-O', '-pg') + _,
    'O': ('-O3',) + _,
}

