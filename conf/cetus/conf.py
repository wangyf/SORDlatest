notes = """
account: jmz
"""
login = 'cetus.alcf.anl.gov'
hosts = 'wangyf',
maxcores = 16
maxram = 65536 
minnodes = 1
maxnodes = 22640 
maxtime = 24, 00
rate = 1e6
fortran_serial = 'gfortran',
fortran_mpi = 'mpif90',

#_ = '-fimplicit-none', '-Wall', '-Wtabs','-fopenmp', '-o'
#_ = '-Wall', '-Wtabs','-qsmp=omp','-o',
_ = '-o',
fortran_flags = {
    'g': _,
    't': _,
    'p': ('-O', '-p',) + _,
    'O': ( '-O3','-g') + _,
}


