#!/bin/bash
bo=$1
migfile=$2
region=$3
migs=`cut -f 1 -d '.' $migfile | sort | uniq`


winFile=windows/filtered_windows_nonoverlapping.bed
if [[ ! -e $winFile ]]; then
    echo "$winFile not found" >> /dev/stderr
    exit 0
fi

if [[ ! -e $bo.bed.gz ]]; then
    smc2bed-all $bo
fi
f=$bo.migStatsHap.bed.gz
if [[ ! -e $f ]]; then
     echo "running arg-summarize on $bo" >> /dev/stderr
     arg-summarize --log-file $bo.log -a $bo.bed.gz --burnin 500 \
         --hap-mig-file $migfile | bgzip > $f
fi
pf=$bo.migStatsHap.part.bed.starch
if [[ ! -e $pf ]]; then
    echo "making $pf" >> /dev/stderr

    startCoord=`awk -v r=$region '$4==r {print $2}' $winFile`
    endCoord=`awk -v r=$region '$4==r {print $3}' $winFile`
    zcat $bo.migStatsHap.bed.gz |
      awk -v startCoord=$startCoord -v endCoord=$endCoord -v OFS="\t" '
         NR > 3 && $2 < endCoord && $3 > startCoord {
            if ($2 < startCoord) {$2=startCoord};
            if ($3 > endCoord) {$3=endCoord};
            print $1,$2,$3}' | sort-bed - | bedops --partition - |
      starch - > $pf
fi
for mig in $migs; do
    inds=`awk -v mig=$mig -F "[. ]" '$1==mig {print $2}' $migfile | sed 's/_[12]$//g' | uniq`
    for ind in $inds; do
        for t in het hom either; do
             tempf=$bo.$mig.$ind.$t.bed.starch
             if [[ -e $tempf ]]; then
                 if [[ $pf -nt $tempf ]]; then
                     rm -f $tempf
                 fi
             fi
             if [[ ! -e $tempf ]]; then
                  echo "making score file for $region $mig $ind $t" >> /dev/stderr
                  zcat $f | awk -v ind=$ind -v mig=$mig -v OFS="\t" -v t=$t '
         NR==3 {col1=-1; col2=-1;
            for (i=5; i<=NF; i++) {if ($i==mig"."ind"_1") {col1=i}; if ($i==mig"."ind"_2") {col2=i}};
            if (col1== -1 || col2==-1) {print col1,col2,ind,mig,"Error finding columns" >> "/dev/stderr"; exit}
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
done
