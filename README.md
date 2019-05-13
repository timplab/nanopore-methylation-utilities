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
- "m" means methylated, "u" means unmethylated, and "x" means uncalled (not confident)
- methylation call corresponds to the motif at position preceding the letter

The resulting bed-style file is sorted, bgzipped, and tabix indexed for easy manipulation.

```
./mtsv2bedGraph.py -i [path/to/nanopolish/methylation.tsv] |\
  sort -k1,1 -k2,2n | bgzip > [methylation.bed.gz]
tabix -p bed [methylation.bed.gz]
```

converting bam for igv
------
Using the original bam file and the converted bed-style methylation file, methylation motifs can be "bisulfite converted _in silico_" for easy visualization on IGV via their bisulfite mode.
There are three options for specifying the region to perform the conversion:
- `-a,--fai` : for the entire genome, supply the fasta fai index
- `-r,--regions` : for multiple regions, supply the bed file
- `-w,--window` : for one region, supply the coordinate (chr:start-end)

```
./convert_bam_for_methylation.py -b [path/to/sorted.bam] -c [path/to/cpg.methylation.bed.gz] \
  -a [path/to/fasta.fai] -o [path/to/converted.bam]
```
