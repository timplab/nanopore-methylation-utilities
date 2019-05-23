#!/bin/bash
root=/kyber/Data/Nanopore/Analysis/gilfunk/190508_recalledDATA_gup_flip/
#cell=MCF7
bam=$root/onTarg_mm2/mcf10a_prim_onTarg.bam
#bam=$root/rawbams_ngmlr/minion_gup3_sorted.bam
cpg=$root/190510_methyLATION/isac_calls/mcf10a_SORT_isac.mcalls.bed.gz
#gpc=$root/pooled/mbed/${cell}_nanoNOMe.pooled.gpc.meth.bed.gz
t=10
#fai=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa.fai
bed=$root/bedz/methyL_2kb_flank.bed
#out=$root/../../igv/${cell}_nanoNOMe.methylation.bam
#log=$root/../../igv/${cell}_nanoNOMe.methylation.log
out=$root/meth_plots/mcf10a_isacHACK_isactry_mm2.bam
ref=/mithril/Data/NGS/Reference/human38/GRCH38.fa

python convert_bam_for_methylation.py --verbose -t $t \
  -f $ref -b $bam -c $cpg -r $bed --all |\
  samtools sort -o $out

samtools index $out

