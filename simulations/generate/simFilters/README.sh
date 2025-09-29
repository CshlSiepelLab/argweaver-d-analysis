#!/bin/bash

## these files are same as those in ~/arghmm/SGDP/filters,
## except they are starched and have "chr" in chrom name

masks="filter.bed.gz Altai_mask.bed.gz Denisova_mask.bed.gz
Vindija_mask.bed.gz panTro4_filter.bed.gz"
maskdir=~/arghmm/SGDP/filters

mkdir -p filters
for mask in $masks; do
  echo "making mask $mask"
  maskfile=$maskdir/$mask
  if [[ ! -e $maskfile ]]; then
    echo "cannot find $maskfile"
    exit 0
  fi
  outfile=filters/`basename $mask .gz`.starch
  zcat $maskfile | awk '{$1="chr"$1; print $0}' |
     sort-bed - | starch - > $outfile
done
