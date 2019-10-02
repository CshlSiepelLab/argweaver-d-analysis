#!/bin/bash

## CRF code was supplied by Sriram Sankararaman, in an executable called
## "caller"

## The script in ../generate/genRecentSims.sh creates input files to
## the caller program, so just need to call it for each simulated data set


dir=recentSims
wd=`pwd`
for i in `seq 1 100`; do
    if [[ -e $dir/$i/crf/par.caller.test ]]; then
	cd $dir/$i/crf
	caller -p par.caller.test
	cd $wd
    fi
done
