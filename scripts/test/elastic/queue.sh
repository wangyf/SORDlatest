#!/bin/bash -e

cd '/gpfs/mira-home/wangyf/local/sord/scripts/test/elastic'

qsub -A GMSeismicSim -t 0:10:00 --proccount 3904 --mode c4 -n 976 -O elastic ./sord-mO

