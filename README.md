# SCE Analysis
This repository contains the R script used for processing Strand-seq data using the R package `BreakpointR` and collecting bed file output of putative SCE breakpoints

### Usage
* Before beginning, generate your [BreakpointR](https://github.com/daewoooo/BreakPointR) output by running this package on your processed Strand-seq BAM files.

* Download this repo and run this command from your terminal to run R script: `Rscript sce_analysis.R /path/to/BreakpointR/data`

* OPTIONAL: generate custom `centromerers2.txt` file using coordinates specific to your cell line
