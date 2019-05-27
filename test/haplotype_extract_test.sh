#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data/hap/all
cell=GM12878
reg=chrX
t=10
outdir=$root/mbed
[ -e $outdir ]||mkdir $outdir

for hap in 1 2; do
  bam=$(find $root/ -name "*$cell*$reg*hap$hap.bam")
  if [ ! -e $bam.bai ]; then
    echo "indexing bam"
    samtools index $bam 
  fi
  for mod in cpg gpc; do
    mbed=$(find $root/../../nanonome/pooled/mbed -name "*$cell*$mod*bed.gz") 
    out=$outdir/${cell}_nanoNOMe.$reg.$mod.hap$hap.meth.bed.gz
    log=$outdir/${cell}_nanoNOMe.$reg.$mod.hap$hap.extract.log
    python ../extract_mbed_by_qname.py --verbose -t $t \
      -b $bam -m $mbed |\
      sort -T $outdir -k1,1 -k2,2n | bgzip > $out
    tabix -p bed $out
  done
done
