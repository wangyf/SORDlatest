#!/bin/bash -e

mode='s'
cd '/Users/yow004/Desktop/SORDlatest/script/asperity'
[ "$mode" = m ] && ./mpd.sh

echo "$( date ): asperity started" >> log

case "$mode${1:--i}" in
    s-i)   time ./sord-sO ;;
    s-g)   gdb  ./sord-sO ;;
    s-ddd) ddd  ./sord-sO ;;
    m-i)   mpiexec -np 1 time ./sord-sO ;;
    m-g)   mpiexec -gdb -np 1 ./sord-sO ;;
    m-ddd) mpiexec -np 1 ddd  ./sord-sO ;;
esac

echo "$( date ): asperity finished" >> log

