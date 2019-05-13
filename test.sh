#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data/nanonome
cell=MCF7
bam=$root/pooled/bam/${cell}_nanoNOMe.pooled.bam
cpg=$root/pooled/mbed/${cell}_nanoNOMe.pooled.cpg.meth.bed.gz
gpc=$root/pooled/mbed/${cell}_nanoNOMe.pooled.gpc.meth.bed.gz
t=10
fai=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa.fai
out=$root/../../igv/${cell}_nanoNOMe.methylation.bam
log=$root/../../igv/${cell}_nanoNOMe.methylation.log

python convert_bam_for_methylation.py --verbose -t $t \
  -b $bam -c $cpg -g $gpc -a $fai -o $out &> $log
