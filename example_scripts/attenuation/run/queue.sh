#!/bin/bash -e

cd '/Users/yow004/Google_Drive_UCSD/code_under_development/SORDlatest/scripts/attenuation/run'
[ "$mode" = m ] && ./mpd.sh

nice nohup ./run.sh > out.log &
pid=$!
echo "$( date ): PID: $pid" >> log
echo "run started with PID: $pid"

