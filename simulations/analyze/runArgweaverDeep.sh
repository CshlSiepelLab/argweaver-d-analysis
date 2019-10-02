#!/bin/bash
set -e
rep=$1
run1=$2
modelFile=$3
subset=$4
recomb=$5
repl=$6

## This script takes up to 6 arguments:
## 1) rep number (which rep of the simulation to analyze)
## 2) run1: the main directory containing all simulations
## 3) modelFile: the model file passed to the --pop-tree-file argument of arg-sample.
##      Assumes the model name is given by `basename $modelFile .txt`
## 4) subset: a string specifying the subset of individuals to use. Usually this is "2afr", though
##      we experimented with 4afr and 8afr. The script will expect a file named
##      subsites_${subset}.txt containing a list of haplotypes from the sites file to keep.
## 5) recomb: should be either "real" (to use real recomb map) or a number (e.g., 5e-9) specifing
##      a constant recombination rate to use in analysis
## 6) repl: an optional string to append to output directory, in case this is a replicate of
##      a previous run


chrX=""
if [[ $run1 =~ '_X' ]]; then
   chrX=1
fi

if [[ $subset == "2afr_vinOnly" ]]; then
  ageFile=sample_ages_vinOnly.txt
else
  ageFile=sample_ages.txt
fi
if [[ -z $recomb ]]; then
  # recomb should be real or number ie 5e-9
  recomb=real
fi
timesfile=times2.txt
popsizefile=pop_sizes2.txt
if [[ -n $chrX ]]; then
  popsizefile=pop_sizes2_chrX.txt
fi
model=`basename $modelFile .txt`
outdir=$run1/$rep/$model.$subset.${recomb}_recomb
mkdir -p $outdir
bo=$outdir/out$repl
resume="--overwrite"
if [[ -e $bo.0.smc.gz ]]; then
  resume="--resume"
fi
recombArg=""
if [[ $recomb == "real" ]]; then
  recombArg="--recombmap $run1/$rep/recomb_map.bed.gz"
elif [[ -n $recomb ]]; then
  recombArg="--recombrate $recomb"
fi
if [[ ! -e $run1/$rep/ind_filters.txt ]]; then
  awk -F "\t" -v d="$run1/$rep" '{print $1,d"/"$2}' ind_filters.txt > $run1/$rep/ind_filters.txt
fi
if [[ ! -e $bo.2000.smc.gz ]]; then
  arg-sample --unphased-file ${subset}_pairs.txt --subsites subsites_${subset}.txt -c 10 --times-file $timesfile --resample-window-iters 3 $recombArg --mutmap $run1/$rep/mu_map.bed.gz -o $bo $resume --start-mig 100 --pop-file pop_assignments_$subset.txt --pop-tree-file $modelFile --popsize-file $popsizefile --age-file $ageFile -n 2000 --sites $run1/$rep/sim.sites.gz --randomize-phase 1.0 --maskmap $run1/$rep/filters/filter.bed.gz --ind-maskmap $run1/$rep/ind_filters.txt --smc-prime
fi
bash finish_argweaver_job.sh $bo migs_$model.txt migs_${model}_${subset}_haps.txt
