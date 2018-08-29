#!/bin/bash -e

mode='m'
cd '/Users/yow004/Google_Drive_UCSD/code_under_development/SORDlatest/scripts/muljob_script/mdl001'
[ "$mode" = m ] && ./mpd.sh

echo "$( date ): mdl001 started" >> log

case "$mode${1:--i}" in
    s-i)   time ./sord-mO ;;
    s-g)   gdb  ./sord-mO ;;
    s-ddd) ddd  ./sord-mO ;;
    m-i)   mpiexec -np 1024 time ./sord-mO ;;
    m-g)   mpiexec -gdb -np 1024 ./sord-mO ;;
    m-ddd) mpiexec -np 1024 ddd  ./sord-mO ;;
esac

echo "$( date ): mdl001 finished" >> log

