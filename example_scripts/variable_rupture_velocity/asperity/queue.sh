#!/bin/bash -e

cd '/Users/yow004/Desktop/SORDlatest/script/asperity'
[ "$mode" = m ] && ./mpd.sh

nice nohup ./run.sh > out.log &
pid=$!
echo "$( date ): PID: $pid" >> log
echo "asperity started with PID: $pid"

