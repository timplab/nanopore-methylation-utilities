#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data/nanonome/pooled
cell=MCF7
bam=$(find $root/bam -name "*$cell*bam")
cpg=$(find $root/mbed -name "*$cell*cpg*bed.gz") 
gpc=$(find $root/mbed -name "*$cell*gpc*bed.gz")
t=10
outdir=$root/igv
out=$outdir/${cell}_nanoNOMe.pooled.methylation.bam
log=$outdir/${cell}_nanoNOMe.pooled.bamconvert.log

python ../convert_bam_for_methylation.py --verbose -t $t \
  -b $bam -c $cpg -g $gpc --all #|\
#  samtools sort -o $out

samtools index $out
