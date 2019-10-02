#!/bin/bash
set -e
indir=$1
modelFile=$2
recomb=$3
subset=nea2afr2eur
if [[ -z $modelFile ]]; then
    ## This is model used in main paper;
    ## the supplement also looks at eurAfr_simpLe_model2.txt
    modelFile=eurAfr_simple_model.txt
fi

if [[ -z $recomb ]]; then
    ## This value was used for analysis in paper
    recomb=5e-9
fi

model=`basename $modelFile .txt`

ageFile=sample_ages_noden.txt
timesfile=times2.txt
popsizefile=pop_sizes_eurafr2.txt
outdir=$indir/$model/${recomb}_recomb_constMut
mkdir -p $outdir
bo=$outdir/out$repl
resume="--overwrite"
if [[ -e $bo.0.smc.gz ]]; then
  resume="--resume"
fi
recombArg=""
if [[ $recomb == "real" ]]; then
  recombArg="--recombmap $indir/recomb_map.bed.gz"
elif [[ -n $recomb ]]; then
  recombArg="--recombrate $recomb"
fi
if [[ ! -e $bo.2000.smc.gz ]]; then
arg-sample --subsites ${subset}.txt -c 10 --times-file $timesfile --resample-window-iters 3 $recombArg -o $bo $resume --start-mig 100 --pop-file pop_assignments_$subset.txt --pop-tree-file $modelFile --popsize-file $popsizefile --age-file $ageFile -n 2000 --sites $indir/sim.sites.gz --smc-prime --mutrate 1.5e-8
fi
bash finish_argweaver_job_recent.sh $bo migs_$model.txt migs_${model}_haps.txt
