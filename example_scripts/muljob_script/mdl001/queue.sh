#!/bin/bash -e

cd '/Users/yow004/Google_Drive_UCSD/code_under_development/SORDlatest/scripts/muljob_script/mdl001'
[ "$mode" = m ] && ./mpd.sh

nice nohup ./run.sh > out.log &
pid=$!
echo "$( date ): PID: $pid" >> log
echo "mdl001 started with PID: $pid"

