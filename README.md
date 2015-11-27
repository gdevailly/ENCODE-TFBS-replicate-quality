# ENCODE-TFBS-replicate-quality
R code for "Variable reproducibility in genome-scale public data: A case study using ENCODE ChIP sequencing resource" article, FEBS letter, by Guillaume Devailly, Anna Mantsoki, Tom Michoel and Anagha Joshi. [Link to the article](http://www.sciencedirect.com/science/article/pii/S0014579315010315)

## How to use
This repository contains R scripts (.R), various tab-delimited text files (.txt and .bed) as well as R binary data (.RData). It DOES NOT contains .bed and .bam files downloaded from ENCODE websites. Howerver, instruction to download those is explained in comments in various R scripts.
* s0_sessionInfo.R contains version of R and packages that were used through the analysis.
* sa_*.R scripts are scripts to generate the various .RData files from ENCODE .bed and .bam files.
* sp_*.R scripts are scripts used to generate the various figure of the article using .RData file. One does not need to download ENCODE .bed and .bam files to generate the plots, but only needs the RData and .txt files. The only exception to this is sp5_CorHeatmapByTF_sf4.R, that require a *heavy* RData object generated for another project, missing in this repository.

## What this repository is for
* in depth details on how the analyses described in the article were done
* attempt to reproduce analyses described in the article

## What this repository is NOT
* It is not a *download and run* project. You will need to properly set up various directory structure, to download files from ENCODE, to manually edit path in R, to install various R packages, and to addapt mc.cores in mclapply() functions according to your computer
* It is not a pedagogical ressource to learn R or bioinformatics. The R code is arguably a rather complete collection of bad R code practices. Many steps could have been done in cleaner or faster ways. It could be noted that I made significant progress in R since the begining of this project.
* It is not aimed to be further updated or improved as it should reflect analyses describe in the article. Exception to this are correction of typos, errors and upload of any missing file.

## Contact
Mail: guillaume.devailly _at_ roslin.ed.ac.uk
