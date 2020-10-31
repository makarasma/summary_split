# summary_split
A Python script that splits variants in a GWAS summary statistic file into two lists: genome-wide significant variants and insignificant variants.

This script was used as a part of my workflow in the master's thesis project where I separated the GWAS variants into two lists.
The genome-wide significant variants are clumped to adjust for linkage disequilibrium. Then the remaining variants are used to clump the non-significant ones. 

Use -s (--summary) flag to select a GWAS summary statistic file.
Use -g (--genotype_data) flag to select genotype data stored in plink format (file trio - .bed, .bim, .fam)

The script requires installed PLINK 1.9.
