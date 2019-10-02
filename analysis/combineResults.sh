#!/bin/bash

## This script was used to take results computed in many separate regions in the genome
## and combine them into a single set of predictions.

## model and subset are the same as given in runDeep.sh or runOOA.sh

model=$1
subset=$2
migfile=migs_${model}_${subset}
lastmig=`awk '{print $1}' $migfile | tail -n 1 | sed 's/_[12]$//g'`
outdir=$model/regions
migs=`cut -f 1 -d '.' $migfile | sort | uniq`

## This just makes sure that we only have results for finished runs
for r in $regions; do
    bo=$model/$r/out
    if [[ ! -e $bo.bed.gz ]]; then continue; fi
    maxrep=`zcat $bo.bed.gz | awk -v maxrep=0 -v startCoord=-1 '$0 !~ /^#/ {if (startCoord < 0) {startCoord=$2};
if ($4 > maxrep) {maxrep=$4};
if ($2 > startCoord) {print maxrep; exit}}'`
    if [[ $maxrep -lt 2000 ]]; then continue; fi
    if [[ -e $bo.$lastmig.either.bed.starch ]]; then
        nl=`unstarch $bo.$lastmig.either.bed.starch | wc -l`
        if [[ $nl == 0 ]]; then
            rm -f $bo.*.starch $bo.migStatsHap.bed.gz
        else
            continue
        fi
    fi
done

mkdir -p $outdir
migs=`cut -f 1 -d '.' $migfile | sort | uniq`

## need to parse through these files and get fractions of homo,het
## for each mig+ind
for mig in $migs; do
    inds=`awk -v mig=$mig -F "[. ]" '$1==mig {print $2}' $migfile | sed 's/_[12]$//g' | uniq`
    for ind in $inds; do
        for t in het hom either; do
            outf=$outdir/$mig.$ind.$t.score.bed.gz
            if [[ -e $outf ]]; then
                files=`find $model -name out.$mig.$ind.$t.bed.starch -newer $outf`
                if [[ -z $files ]]; then continue; fi
            fi
            echo "making score file for $model $mig $ind $t"
            ( for r in $regions; do
                  tempf=$model/$r/out.$mig.$ind.$t.bed.starch
                  if [[ ! -e $tempf ]]; then
                      continue
                  fi
                  unstarch $tempf
              done ) | merge_bed.sh | bgzip > $outdir/$mig.$ind.$t.score.bed.gz
        done
    done
done

## Create file showing which regions have finished results
needCovered=1
if [[ -e $outdir/covered.bed ]]; then
    files=`find $model -name out.$lastmig.either.bed.starch -newer $outdir/covered.bed`
    if [[ -z $files ]]; then
        needCovered=
    fi
fi
if [[ -n $needCovered ]]; then
    echo "making $outdir/covered.bed" >> /dev/stderr
    ( for r in $regions; do
          pf=$model/$r/out.$lastmig.either.bed.starch
          if [[ -e $pf ]]; then
              unstarch $pf | bedops -m -
          fi
      done ) | bedops -m - > $outdir/covered.bed
fi

## after running the above, we should have a file for each mig, ind, and het/hom
## giving the fraction of runs (after burnin) that have a migration of that type in that ind
## make file for each potential cutoff
for cutoff in 0.5 0.9; do
    ## first homozygous is easy
    for t in hom het either; do
        for f in $outdir/*.$t.score.bed.gz; do
            bn=`basename $f .score.bed.gz`
            zcat $f | awk -v cutoff=$cutoff -v OFS="\t" '$4 >= cutoff {print $1,$2,$3}' |
                bedops -m - > $outdir/$bn.$cutoff.bed
        done
    done
done

awk '$1 != "X"' $outdir/covered.bed > $outdir/covered.autosome.bed
for i in `seq 1 22` X; do
    awk -v chr=$i '$1==chr' $outdir/covered.bed > $outdir/covered.$i.bed
done

## report coverage, element lengths in X and autosomes
for chr in X autosome; do
    cov=$outdir/covered.$chr.bed
    totcov=`sorted-bed-sum $cov`
    if [[ $totcov == 0 ]]; then continue; fi
    echo "$model run $chr total coverage=$totcov"
    for mig in $migs; do
        cov5=`bedops -u $outdir/$mig.*.either.0.5.bed | bedops -i - $cov | sorted-bed-sum`
        cov9=`bedops -u $outdir/$mig.*.either.0.9.bed | bedops -i - $cov | sorted-bed-sum`
        numel5=`bedops -m $outdir/$mig.*.either.0.5.bed | bedops -i - $cov | wc -l`
        numel9=`bedops -m $outdir/$mig.*.either.0.9.bed | bedops -i - $cov | wc -l`
        if [[ $numel5 -gt 0 ]]; then
            meanlen5=`echo  $cov5/$numel5 | bc`
        else meanlen5=NA; fi
        if [[ $numel9 -gt 0 ]]; then
            meanlen9=`echo $cov9/$numel9 | bc`
        else meanlen9=0; fi
        cov5=`echo $cov5/$totcov | bc`
        cov9=`echo $cov9/$totcov | bc`
        echo -e "$mig\t$cov5\t$cov9\t$numel5\t$numel9\t$meanlen5\t$meanlen9"
        inds=`awk -v mig=$mig -F "[. ]" '$1==mig {print $2}' $migfile | sed 's/_[12]$//g' | uniq`
        for ind in $inds; do
            cov5=`bedops -u $outdir/$mig.$ind.either.0.5.bed | bedops -i - $cov | sorted-bed-sum`
            cov9=`bedops -u $outdir/$mig.$ind.either.0.9.bed | bedops -i - $cov | sorted-bed-sum`
            numel5=`bedops -m $outdir/$mig.$ind.either.0.5.bed | bedops -i - $cov | wc -l`
            numel9=`bedops -m $outdir/$mig.$ind.either.0.9.bed | bedops -i - $cov | wc -l`
            if [[ $numel5 -gt 0 ]]; then
                meanlen5=`echo $cov5/$numel5 | bc`
            else meanlen5=NA; fi
            if [[ $numel9 -gt 0 ]]; then
                meanlen9=`echo $cov9/$numel9 | bc`
            else meanlen9=NA; fi
            cov5=`echo $cov5/$totcov | bc`
            cov9=`echo $cov9/$totcov | bc`
            echo -e "\t$ind\t$cov5\t$cov9\t$numel5\t$numel9\t$meanlen5\t$meanlen9"
        done
    done
done

