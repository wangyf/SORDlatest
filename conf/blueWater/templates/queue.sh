#!/bin/bash -e

cd %(rundir)r
if [ $( /bin/pwd | grep -v scratch ) ]; then
    echo "Error: jobs must be run from /scratch"
    exit
fi

echo "$( date ): %(name)s queued with ID: $( qsub script.sh )" >> log

