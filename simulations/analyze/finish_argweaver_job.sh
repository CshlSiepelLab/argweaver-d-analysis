#!/bin/bash
set -e

## This script parses the results from arg-sample and creates files that summarize
## the introgression calls across the genomic region. It is called by runArgweaverDeep.sh

## The first argument should be the "out-root" (bo = baseout) argument given to arg-sample,
##   so there should exist files $bo.stats, $bo.log, $bo.*.smc.gz, etc
## The second argument is the "migfile" file which summarizes the migration events in
##   the model used by arg-sample, and is passed to the "--mig-file" option in arg-summarize
## The third argument is the "hapmigfile" which summarizes which haplotypes can take which
##   migrations in the model used by arg-sample, and is passed to the "--hap-mig-file" option
##   in arg-summarize


bo=$1
migfile=$2
hapmigfile=$3


## First, call smc2bed-all; this combines all smc files into one big sorted bed file
smc2bed-all $bo

## This command produces a file giving the probability of each type of migration at each site
arg-summarize --burnin 500 --mean --mig-file $migfile --log-file $bo.log -a $bo.bed.gz |
    merge_bed.sh |
    bgzip > $bo.migStats.bed.gz

## This command is the same as above, but does not average over MCMC replicates. Instead it
##   indicates whether there is each type of migration at each site for each MCMC rep
arg-summarize --mig-file $migfile --log-file $bo.log -a $bo.bed.gz |
    bgzip > $bo.migstats.bed.gz


## The rest of the file deals with whether specific haplotypes are introgressed.
## First get a list of migration events. They are in the first column of $hapmigfile
migs=`cut -f 1 -d '.' $hapmigfile | sort | uniq`


f=$bo.migStatsHap.bed.gz
echo "running arg-summarize on $bo" >> /dev/stderr
arg-summarize --log-file $bo.log -a $bo.bed.gz --burnin 500 \
     --hap-mig-file $hapmigfile | bgzip > $f

pf=$bo.migStatsHap.part.bed.starch
echo "making $pf" >> /dev/stderr
zcat $bo.migStatsHap.bed.gz | sort-bed - | bedops --partition - |
  starch - > $pf

for mig in $migs; do
    inds=`awk -v mig=$mig -F "[. ]" '$1==mig {print $2}' $hapmigfile | sed 's/_[12]$//g' | uniq`
    for ind in $inds; do
        for t in het hom either; do
             tempf=$bo.$mig.$ind.$t.bed.starch
             if [[ ! -e $tempf ]]; then
                  echo "making score file for $region $mig $ind $t" >> /dev/stderr
                  zcat $f | awk -v ind=$ind -v mig=$mig -v OFS="\t" -v t=$t '
         NR==3 {col1=-1; col2=-1;
            for (i=5; i<=NF; i++) {if ($i==mig"."ind"_1") {col1=i}; if ($i==mig"."ind"_2") {col2=i}};
            if (col1== -1 || col2==-1) {print "Error finding columns" >> "/dev/stderr"; exit}
          };
          NR > 3 {
            score=0
            if (t=="hom" && $col1==1 && $col2==1) {score=1};
            if (t=="het" && xor($col1, $col2)) {score=1};
            if (t=="either" && ($col1==1 || $col2==1)) {score=1};
            print $1,$2,$3,$4,score}' | sort-bed - |
                    bedmap --echo --mean --delim "\t" $pf - | merge_bed.sh |
                   starch - > $tempf
             fi
        done
    done
    cmd='bedops -p $bo.$mig.*.either.bed.starch'
    for ind in $inds; do
        cmd="$cmd | bedmap --delim \"\t\" --echo --echo-map - $bo.$mig.$ind.either.bed.starch"
    done
    eval $cmd | awk -v OFS="\t" '{maxscore=$7;
                                 for (i=11; i<=NF; i+=4) {
                                    if ($i > maxscore) {maxscore=$i}};
                                 print $1,$2,$3,maxscore}' |
       merge_bed.sh | starch - > $bo.$mig.any.bed.starch
done


## make summary file with all migs
files=`ls $bo.*.any.bed.starch`
cmd='bedops -p $files'
header="#chrom\tchromStart\tchromEnd"
cols=1,2,3
nextcol=7
for mig in $migs; do
    cmd="$cmd | bedmap --delim \"\t\" --echo --echo-map - $bo.$mig.any.bed.starch"
    header="$header\t$mig"
    cols="$cols,$nextcol"
    nextcol=$(($nextcol+4))
done
echo $cmd
( echo -e $header
eval $cmd | cut -f $cols ) | bgzip > $bo.migStats2.bed.gz
