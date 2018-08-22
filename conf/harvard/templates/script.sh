#!/bin/bash

#SBATCH -n %(np)s                    # Number of cores
#SBATCH -N %(nodes)s                    # Ensure that all cores are on one machine
#SBATCH -t %(walltime)s              # Runtime in D-HH:MM
#SBATCH -p general       # Partition to submit to
#SBATCH --contiguous     #Ensure that all of the cores are on the same Infiniband network
#SBATCH --mem-per-cpu=1024               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o stdout      # File to which STDOUT will be written
#SBATCH -e stderr      # File to which STDERR will be written
#SBATCH --mail-type=BEGIN,END,FAIL,ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yow004@ucsd.edu # Email to which notifications will be sent

mpiexec -np %(np)s %(bin)s
