#!/bin/bash
root=/kyber/Data/Nanopore/Analysis/gilfunk/190508_recalledDATA_gup_flip/
#cell=MCF7
bam=$root/onTarg_mm2/min_gup3_prim_onTarg.bam
cpg=$root/190510_methyLATION/isac_calls/min_gm12878_SORT_isac.mcalls.bed.gz
#gpc=$root/pooled/mbed/${cell}_nanoNOMe.pooled.gpc.meth.bed.gz
t=10
#fai=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa.fai
bed=$root/meth_bed/methyL_2kb_flank.bed
#out=$root/../../igv/${cell}_nanoNOMe.methylation.bam
#log=$root/../../igv/${cell}_nanoNOMe.methylation.log
out=$root/meth_plots/min_gm12878_isacHACK_isactry.bam
ref=/mithril/Data/NGS/Reference/human38/GRCH38.fa

python convert_bam_for_methylation.py --verbose -t $t \
  -f $ref -b $bam -c $cpg -r $bed -o $out

