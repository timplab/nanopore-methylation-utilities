## nanopore-methylation-utilities
Set of utilities for analyzing nanopore methylation data from the Timp Lab

- [Bed-Style format methylation File](#bedstyle)
- [BAM conversion for Methylation Viewing in IGV](#igv)
- [Citation](#cite)

# <a name="bedstyle"></a> bed-style format methylation file


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

# <a name="igv"></a>converting bam for igv


Using the converted bed-style methylation file, the original bam file can be "bisulfite converted _in silico_" for easy visualization on IGV via their bisulfite mode.
There are three options for specifying the region to convert:
- `-r,--regions` : for multiple regions, supply the bed file
- `-w,--window` : for one region, supply the coordinate (chr:start-end)
- without either of the above options, all reads will be converted

```
./convert_bam_for_methylation.py -b [path/to/sorted.bam] \
  -c [path/to/cpg.methylation.bed.gz] -f [path/to/reference.fasta ] |\
  samtools sort -o [path/to/converted.bam]
samtools index [path/to/cnverted.bam]
```
#### For minimap2 alignments 

Using `--MD` option during alignment is recommended.

The default output does not have MD tags, and MD tags are necessary for using pysam to get the reference sequence. To get around this, the fasta of reference genome must be supplied via `-f,--fasta`.

# <a name="cite"></a> Citation

If you use this package in your work, please cite:

> Lee, I. et al. (2019). Simultaneous profiling of chromatin accessibility and methylation on human cell lines with nanopore sequencing.
> *bioRxiv*. [doi:10.1101/504993v2][doi]


[doi]: https://doi.org/10.1101/504993
