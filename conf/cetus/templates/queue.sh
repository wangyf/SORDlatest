#!/bin/bash -e

cd %(rundir)r

qsub -A GMSeismicSim -t 10 --proccount %(np)s --mode c%(ppn)s -n %(nodes)s -O %(name)s %(bin)s
#qsub --env BGLOCKLESSMPIO_F_TYPE=0x47504653 -A SeismicHazard -t 60 --proccount %(np)s --mode c%(ppn)s -n %(nodes)s -O %(name)s %(bin)s
