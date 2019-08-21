#!/bin/bash

# for version 4 mapping that came in cram format

pwdir=$(readlink -f ../)
dir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly
outdir=$dir/isac
cd $dir
[ -e $outdir ]||mkdir $outdir

if [ "$1" == "meth" ]; then
	python3 $pwdir/mtsv2bedGraph.py -i methylation.tsv |\
		sort -k1,1 -k2,2n | bgzip > methylation.bed.gz
	tabix -p bed methylation.bed.gz

fi

if [ "$1" == "bam" ]; then
  outbam=$outdir/rel2_to_v0.4_DXZ4_isac_nopoor.bam 
  python3 $pwdir/convert_bam_for_methylation.py --debug --verbose -b DXZ4.bam \
    -c methylation.bed.gz -f chm13.draft_v0.4.fasta \
    -w chrX_fixedBionanoSV_centromereV3:112000000-114000000 \
    --remove_poor |\
    samtools sort -o $outbam
  samtools index $outbam
fi
