# nanopore-methylation-utilities
Set of utilities for analyzing nanopore methylation data

bed-style format methylation file
------
I convert the nanopolish methylation calling output into bed-style format, such that each line is

|Contig |Start  |End  |Read name  |Methylation call string  |Log-likelihood ratios  |Motif context  |
|-------|-------|-----|-----------|-------------------------|-----------------------|---------------|

where Methylation call string is arranged such that 
- numbers are separated by methylation calls
- each number is cumulative distance from the "start"
- methylation call corresponds to the motif at position preceding the letter
- "m" means methylated, "u" means unmethylated, and "x" means uncalled (not confident)

The resulting bed-style file is sorted, [bgzipped](http://www.htslib.org/doc/bgzip.html), and [tabix](http://www.htslib.org/doc/tabix.html) indexed for easy manipulation.  
```
./mtsv2bedGraph.py -i [path/to/nanopolish/methylation.tsv] |\
  sort -k1,1 -k2,2n | bgzip > [methylation.bed.gz]
tabix -p bed [methylation.bed.gz]
```

converting bam for igv
------
Using the converted bed-style methylation file, the original bam file can be "bisulfite converted _in silico_" for easy visualization on IGV via their bisulfite mode.
There are three options for specifying the region to convert:
- `-a,--fai` : for the entire genome, supply the fasta fai index
- `-r,--regions` : for multiple regions, supply the bed file
- `-w,--window` : for one region, supply the coordinate (chr:start-end)

For minimap2 alignments : the default output does not have MD tags, and MD tags are necessary for using pysam to get the reference sequence. To get around this, the fasta of reference genome must be supplied via `-f,--fasta`.
This also makes the process slower and memory intensive than using ngmlr reads or aligning minimap2 with `--MD` option

```
./convert_bam_for_methylation.py -b [path/to/sorted.bam] \
  -c [path/to/cpg.methylation.bed.gz] -a [path/to/fasta.fai] |\
  samtools sort -o [path/to/converted.bam]
samtools index [path/to/cnverted.bam]
```
