#!/bin/bash

awk -v OFS="\t" -v prevchr="" '
  BEGIN{prevval[4]=""};
  $0 ~ /^#/ {print $0};
  $0 !~ /^#/ { if (prevchr==$1 && prevend==$2) {
                 matches=1;
                 for (j=4; j<=NF; j++) {
                     if (prevval[j] != $j) {matches=0};
                 }
               } else {matches=0}
               if (matches==1) {
                 prevend=$3;
               } else {
                 if (prevval[4] != "") {
                   printf("%s\t%s\t%s\t%s", prevchr, prevstart, prevend, prevval[4]);
                   for (j=5; j<=NF; j++) {printf("\t%s", prevval[j])};
                   printf("\n");             
                 }
                 prevchr=$1; prevstart=$2; prevend=$3;
                 for (j=4; j<=NF; j++) {prevval[j]=$j};
               }
             }
             END {if (prevval[4] != "") {
                   printf("%s\t%s\t%s\t%s", prevchr, prevstart, prevend, prevval[4]);
                   for (j=5; j<=NF; j++) {printf("\t%s", prevval[j])};
                   printf("\n");             
}}'
