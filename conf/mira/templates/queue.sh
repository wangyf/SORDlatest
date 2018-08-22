#!/bin/bash -e

cd %(rundir)r

qsub -A GMSeismicSim -t %(walltime)s --proccount %(np)s --mode c%(ppn)s -n %(nodes)s -O %(name)s %(bin)s

