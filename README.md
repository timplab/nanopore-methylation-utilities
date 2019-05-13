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


