#!/bin/bash

## This is helper script used with genDeepSims.sh
## After generating simulated data using some true mutation and recombination
## rate maps from hg19, it creates a map file to be used with the simulated
## data by shifting the coordinates to the range of the simulated data.

## dir should have a file called used_rates.bed
## will output files in dir with mask.bed.gz, Altai.bed.gz, Vindija.bed.gz,
## Denisova.bed.gz, panTro4.bed.gz

maskdir=simFilters
masks="filter.bed.starch Altai_mask.bed.starch Denisova_mask.bed.starch
Vindija_mask.bed.starch panTro4_filter.bed.starch"


dir=$1
chr=`zcat $dir/sim.sites.gz | head -n 2 | grep REGION | awk '{print $2}'`
inbed=$dir/used_rates.bed

outdir=$dir/filters
mkdir -p $outdir

if [[ ! -e $inbed ]]; then
    echo "Error: Do not see file $inbed"
    exit 0
fi

for mask in $masks; do
    maskfile=$maskdir/$mask
    if [[ ! -e $maskfile ]]; then
	echo "Error: Could not find $maskfile"
	exit 0
    fi
    outfile=$outdir/`basename $mask .starch`.gz
    echo "Making $outfile"
    bedops -m $inbed |
	bedmap --delim "\t" --multidelim "\t" --echo --echo-map - $maskfile |
	gawk -v OFS="\t" -v chr=$chr -v startPos=0 '{
offset=$2-startPos
for (i=1; i<NF; i+=3) {
  $i=chr
  $(i+1) = $(i+1)-offset
  $(i+2) = $(i+2)-offset
}
startPos = startPos + $3-$2
for (i=4; i<NF; i+=3) {
  if ($(i+1) < $2) {$(i+1)=$2};
  if ($(i+2) > $3) {$(i+2)=$3};
  print $i,$(i+1),$(i+2)
}
}' | bgzip > $outfile
done

