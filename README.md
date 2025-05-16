# RBCLIB

Utility functions to manipulated BCF/VCF files in R usig `htslib` and `bctools` C libraries and a wrapper of the `bcftools` command line interface.

The main motivating gaol is to ultimately reproduce the fantastic [bcftools munge](https://github.com/freeseek/score) utility as a faster alternative to the great [{MungeSumstats}](https://github.com/Al-Murphy/MungeSumstats). 

TODO

- [ ] Provide full wrapping of BCFTOOLS command line tool from R
  
  The strategy is copied from [pysam](https://github.com/pysam-developers/pysam), especially the pysam.h and pysam.c file. The adaptation here of is static linking.
  
  - [ ] Non Interactive mode with Rscript
 
  - [ ] Interactive mode : issue here is how to manage sdtout
 
  - [ ] Allow user interupts : this may be done by running bcftools on a separate thread or background R session. 

- [ ] Basic VCF/BCF manipulations and loading in R
  
  Our goal here is to provide a more dataframe like access to BCF/VCF files. But we will start with basic scanning functions and/or streaming.

- [ ] faidx functions for fasta manipulation

- [ ] kfunctions for intervals

- [ ] Tabix Tabular imports to BCF/VCF

