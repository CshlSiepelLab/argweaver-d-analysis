#!/bin/bash
## This is how each region was analyzed with ARGweaver to look for "deep" introgression
## between Africans, Neanderhtals, Denisovans, and super-archaic hominins.
## The main models used in the paper are mig250_div1000 and mig150_div1500
## The subset used was "2afr" (which includes Mandenka and Khomani_San Africans,
## as well as all the archaics).

## The regions are listed in windows/filtered_windows.bed; the fourth column
## gives all the region name, which should be the first argument to this script

regionName=$1
modelFile=$2
subset=$3

model=`basename $modelFile .txt`

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
hapMigFile=migs_${model}_${subset}_haps.txt
modelFile=${model}Model.txt
if [[ -z $subset ]]; then
  outdir=$model/$regionName
  vcfFiles=vcf_files.txt
  nameMap=name_map.txt
else
  outdir=${model}_$subset/$regionName
  vcfFiles=vcf_files_$subset.txt
  nameMap=name_map_$subset.txt
fi
bo=out

modelArgs="--start-mig 100 --pop-file $popAssignFile --pop-tree-file $modelFile --popsize-file $popSizeFile"

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

arg-sample --vcf-files $vcfFiles -c $compress \
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
       -n 2000 --smc-prime

smc2bed-all $outdir/$bo
bash summarizeOneRun.sh $outdir/$bo $hapMigFile $regionName

