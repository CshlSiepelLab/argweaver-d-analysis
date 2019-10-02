#!/bin/bash

## Call msprime to generate ARGs
## Then use SeqGen to generate mutations conditional on simulated ARGs
## Use mutation rates, recombination rates, and missing data patterns
##  drawn from random autosomal regions of autosomal hg19

#export PYTHONPATH=~/msprime:$PYTHONPATH

## test for seq-gen
tmp=`which seq-gen`
if [[ -z $tmp ]]; then
    export PATH="$HOME/Seq-Gen.v1.3.3/source:$PATH"
    tmp=`which seq-gen`
    if [[ -z $tmp ]]; then
	echo "Error: need to install seq-gen and/or add to PATH"
	exit 1
    fi
fi

## Also assumed installed:
## R
## fa2sites (this comes with argWeaver)
## bedops/bedmap command
## bgzip



migtime=250
divtime=1000
startrep=1
endrep=100
chrX=""
meanMutRate=1.45e-8
numAfr=4
numAlt=2
numVin=2
numDen=2
numArc=1
numSup=1
numChimp=1
supPopSize=500
nomig=""
recomb=""
mainoutdir=""

while getopts ":hnm:d:a:s:e:o:x" opt; do
    case ${opt} in
	h )
	    echo "usage: bash genDeepSims.sh"
	    echo " -h: help (print this message)"
	    echo " -m: migration time (in kya; default 250)"
	    echo " -d: super-archaic divergence time (in kya; default 1000)"
	    echo " -a: number of haploid african samples (default: 4)"
	    echo " -x: chrX simulation (multiply popsizes by 3/4)"
	    echo " -n: set all migration rates to zero"
	    echo " -r: give constant recombination rate to use (otherwise, will use a section of the human recombination map)"
	    echo " -o: specify output directory (default: migMIG_divDIV with _X appended if -x option used and _nomig appended if -n used)"
	    exit 0
	    ;;
	m )
	    migtime=$OPTARG
	    ;;
	d )
	    divtime=$OPTARG
	    ;;
	a )
	    numAfr=$OPTARG
	    ;;
	s )
	    startrep=$OPTARG
	    ;;
	e )
	    endrep=$OPTARG
	    ;;
	n )
	    nomig="--no-mig"
	    ;;
	x)
	    chrX="--chrX"
	    ;;
	r)
	    recomb="--recomb $OPTARG"
	    ;;
	o)
	    mainoutdir=$OPTARG
	    ;;
	\? )
	    echo "Invalid option $OPTARG"
	    exit 1
	    ;;
	: )
	    echo "Invalid option: $OPTARG requires an argument"
	    exit 1
    esac
done

echo "migtime = $migtime kya"
echo "divtime = $divtime kya"
echo "numAfr  = $numAfr"
echo "start   = $startrep"
echo "end     = $endrep"
if [[ -n $chrX ]]; then echo "chrX turned on" else echo "chrX turned off"; fi
if [[ -n $nomig ]]; then echo "nomig turned on"; else echo "nomig turned off"; fi

if [[ -z $mainoutdir ]]; then
    mainoutdir=mig${migtime}_div${divtime}
    if [[ -n $chrX ]]; then mainoutdir=${mainoutdir}_X; fi
    if [[ -n $nomig ]]; then mainoutdir=${mainoutdir}_nomig; fi
fi

numHap=$(($numAfr + $numAlt + $numVin + $numDen + $numSup + $numArc + $numChimp ))

## output will be in order african, altai, vindija, denisova, super-archaic, chimp
## this gets the indicies of first and last samples of each type
afrStart=1
afrEnd=$numAfr
altStart=$(($afrEnd+1))
altEnd=$(($altStart + $numAlt - 1))
vinStart=$(($altEnd+1))
vinEnd=$(($vinStart+$numVin-1))
denStart=$(($vinEnd+1))
denEnd=$(($denStart+$numDen-1))
supStart=$(($denEnd+1))
supEnd=$(($supStart+$numSup-1))
arcStart=$(($supStart+1))
arcEnd=$(($arcStart+$numArc-1))
chimpStart=$(($arcEnd+1))
chimpEnd=$(($chimpStart+$numChimp-1))

for rep in `seq $startrep $endrep`; do
    outdir=$mainoutdir/$rep
    mkdir -p $outdir

    seedfile=$outdir/seeds.txt
    echo 'write(file="'$seedfile'",
                round(runif(3, min=1, max=1000000)), sep="\n")' | R --vanilla
    seed1=`awk 'NR==1' $seedfile`
    seed2=`awk 'NR==2' $seedfile`
    seed3=`awk 'NR==3' $seedfile`

    python simDeep.py --migtime $migtime --divtime $divtime \
	   --afr $numAfr --altai $numAlt --vindija $numVin \
           --denisova $numDen --super $numSup --chimp $numChimp --arc $numArc \
           -o $outdir -c chr$rep -s $seed1 --super-size $supPopSize $chrX $nomig $recomb \
       > $outdir/msprime_script_out.txt
    bgzip -f $outdir/recomb_map.bed
    bgzip -f $outdir/mu_map.bed

    ## use SeqGen to add mutations to generated local trees
    ## the ugly awk command translates between sample numbers and produces more intutitive
    ## haploid sample names (alt1, alt2, vin1, vin2, afr1, afr2, etc)

    ## seqgen requires the number of non-recombining blocks; this is nump
    nump=`awk '$0 ~ /^\[/ {c++}; END{print c}' $outdir/trees.txt`
    ## currMu is the total mutation rate across all sites
    currMu=`zcat $outdir/mu_map.bed | awk '{len=$3-$2; totlen+=len; s+=$4*len}; END{print s/totlen}'`
  seq-gen -q -z $seed3 -p $nump -mHKY -t3.0 \
    -f0.3,0.2,0.2,0.3 -l2000000 -s $currMu \
   < $outdir/trees.txt |
   awk 'NR > 1' |
   awk -v afrStart=$afrStart -v afrEnd=$afrEnd \
       -v altStart=$altStart -v altEnd=$altEnd \
       -v vinStart=$vinStart -v vinEnd=$vinEnd \
       -v denStart=$denStart -v denEnd=$denEnd \
       -v supStart=$supStart -v supEnd=$supEnd \
       -v arcStart=$arcStart -v arcEnd=$arcEnd \
       -v chimpStart=$chimpStart -v chimpEnd=$chimpEnd '
      $1 >= afrStart && $1 <= afrEnd {ind=sprintf("afr%i", $i-afrStart+1)};
      $1 >= altStart && $1 <= altEnd {ind=sprintf("alt%i", $i-altStart+1)};
      $1 >= vinStart && $1 <= vinEnd {ind=sprintf("vin%i", $i-vinStart+1)};
      $1 >= denStart && $1 <= denEnd {ind=sprintf("den%i", $i-denStart+1)};
      $1 >= supStart && $1 <= supEnd {ind=sprintf("sup%i", $i-supStart+1)};
      $1 >= arcStart && $1 <= arcEnd {ind=sprintf("arc%i", $i-arcStart+1)};
      $1 >= chimpStart && $1 <= chimpEnd {
         if (chimpStart==chimpEnd) {ind="chimp"} else {ind=sprintf("chimp%i", $i-chimpStart+1)}}
    {print ">"ind};
    {print $2}' | fa2sites -c chr$rep --start 1 -f /dev/stdin \
         -o $outdir/sim.sites
    gzip -f $outdir/sim.sites

   # get true trees:
    echo $rep
    cp $outdir/trees.txt $outdir/trees_named.txt
    msfile=$outdir/trees_named.txt
    idx=1
    for i in `seq 1 $numHap`; do
       if [[ $i == $afrStart && $numAfr -gt 0 ]]; then str="afr"; idx=1; fi
       if [[ $i == $altStart && $numAlt -gt 0 ]]; then str="alt"; idx=1; fi
       if [[ $i == $vinStart && $numVin -gt 0 ]]; then str="vin"; idx=1; fi
       if [[ $i == $denStart && $numDen -gt 0 ]]; then str="den"; idx=1; fi
       if [[ $i == $supStart && $numSup -gt 0 ]]; then str="sup"; idx=1; fi
       if [[ $i == $arcStart && $numArc -gt 0 ]]; then str="arc"; idx=1; fi
       if [[ $i == $chimpStart && $numChimp -gt 0 ]]; then str="chimp"; idx=""; fi
    sed -i "s/\([(,]\)$i:/\1${str}$idx:/g" $msfile
    idx=$(($idx+1))
    done
    awk -F "[][]"  \
    -v chr=chr$rep -v start=0 -v rep=0 -v OFS="\t" \
    'NF==3 && $2 != 0 {print chr,start,start+$2,rep,$3; start=start+$2}' $msfile |
	bgzip > $outdir/trueArg.bed.gz
    bash used_rates_to_mask.sh $outdir
    rm -f $outdir/trees.txt $outdir/trees_named.txt
done

