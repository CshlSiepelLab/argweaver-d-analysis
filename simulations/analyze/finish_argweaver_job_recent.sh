#!/bin/bash
set -e
## similar to finish_condor_job.sh but simpler because we are dealing with phased genomes

bo=$1
migstats=$2
migfile=$3
if [[ ! -e $bo.bed.gz || $bo.stats -nt $bo.bed.gz ]]; then
    if [[ -e $bo.2000.smc.gz ]]; then
       smc2bed-all $bo
    fi
fi
if [[ ! -e $bo.migStats.bed.gz || $bo.bed.gz -nt $bo.migStats.bed.gz ]]; then
arg-summarize --burnin 500 --mean --mig-file $migstats --log-file $bo.log -a $bo.bed.gz | merge_bed.sh | bgzip > $bo.migStats.bed.gz
fi

f=$bo.migStatsHap.bed.gz
if [[ ! -e $f ]]; then
     echo "running arg-summarize on $bo" >> /dev/stderr
     arg-summarize --log-file $bo.log -a $bo.bed.gz --burnin 500 \
         --hap-mig-file $migfile --burnin 500 --mean | merge_bed.sh | bgzip > $f
fi
