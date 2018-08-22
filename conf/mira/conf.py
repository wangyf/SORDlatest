notes = """
account: jmz
"""
login = 'mira.alcf.anl.gov'
hosts = 'wangyf',
maxcores = 16
maxram = 65536 
minnodes = 1
maxnodes = 22640 
maxtime = 24, 00
rate = 1e6
fortran_serial = 'xlf90',
fortran_mpi = 'mpixlf90',
#_ = '-fimplicit-none', '-Wall', '-o'
_ = '-o',
fortran_flags = {
    'g': _,
    't': _,
    'p': ('-O', '-p',) + _,
    'O': ( '-O3','-g') + _,
}


