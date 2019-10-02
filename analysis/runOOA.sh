#!/bin/bash

## This is how each region was analyzed with ARGweaver to get introgressed
## regions from Nea/Den into modern non-African humans.
## The model used was "ooaModel_nosuper" (ooa for Out-of-Africa;
## nosuper refers to the lack of super-archaic introgression in the model).
## All runs included the Vindija and Altai Neanderthal and the Denisovan.
## The subsets were 1pap2afr (1 papuan and 2 africans), and 1fra2afr (1 french 2 afr)

## regionName taken from 4th column of windows/filtered_windows.bed, i.e., chr10.4
regionName=$1
model=$2
subset=$3

if [[ -z $regionName ]]; then
    echo "error; no regionName given"
    exit 1
fi

export PATH=~/argweaver.github/bin:~/bin:$PATH
export PYTHONPATH=~/argweaver.github:$PYTHONPATH

winfile=windows/filtered_windows.bed
chrom=`awk -v region=$regionName '$4==region {print $1}' $winfile`
chromStart=`awk -v region=$regionName '$4==region {print $2+1}' $winfile`
chromEnd=`awk -v region=$regionName '$4==region {print $3}' $winfile`
region=${chrom}:${chromStart}-${chromEnd}
echo "regionName=$regionName region=$region"

popAssignFile=pop_assignments.txt
if [[ $chrom == "X" ]]; then
  popSizeFile=pop_sizes2_chrX.txt
else
  popSizeFile=pop_sizes2.txt
fi
timesFile=times2.txt
indMaskFile=ind_masks.txt
ageFile=sample_ages.txt
modelFile=$model.txt
popAssignFile=pop_assignments_nosuper.txt
if [[ $chrom == "X" ]]; then
  popSizeFile=pop_sizes_nosuper_chrX.txt
else
  popSizeFile=pop_sizes_nosuper.txt
fi
migFile=migs_${model}_${subset}.txt
outdir=${model}_$subset/$regionName
vcfFiles=vcf_files_$subset.txt
nameMap=name_map_$subset.txt

bo=out

# do not resume; just quite if a run is started. Will have to delete this later
#if [[ -e $outdir/$bo.0.smc.gz ]]; then
#  exit 0
#fi

if [[ $model == "panmictic" ]]; then
  modelArgs="--popsize 10000"
else
  modelArgs="--start-mig 100 --pop-file $popAssignFile --pop-tree-file $modelFile --popsize-file $popSizeFile"
fi

mkdir -p $outdir

if [[ -e $outdir/$bo.2000.smc.gz ]]; then
  echo "already done $outdir/$bo"
  exit 0
fi

if [[ -e $outdir/$bo.0.smc.gz ]]; then
  overwrite="--resume"
else
  overwrite="--overwrite"
fi

recombArg="--recombrate 5e-9"

arg-sample --vcf-files $vcfFiles -c 10 \
       --region $region \
       --times-file $timesFile \
       --resample-window-iters 3 \
       --rename-seqs $nameMap \
       --ind-maskmap $indMaskFile --vcf-genotype-filter "DP<20;DP>80" \
       --vcf-min-qual 20 --mutmap rate_files/subst_rate.bed.gz \
       $recombArg -o $outdir/$bo $overwrite \
       $modelArgs \
       --age-file $ageFile \
       --maskmap filters/filter.bed.gz \
       -n 2000 $optArg

smc2bed-all $outdir/$bo
bash summarizeOneRun.sh $outdir/$bo $migFile $regionName
