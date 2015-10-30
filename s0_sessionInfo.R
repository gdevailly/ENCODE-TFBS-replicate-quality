library(GenomicRanges)
library(Repitools)
library(readr)
library(magrittr)
library(parallel)
library(motifStack)
library(lsr)
library(gplots)
library(RColorBrewer)
library(cba)
library(SDMTools)

sessionInfo()
# R version 3.2.0 (2015-04-16)
# Platform: x86_64-unknown-linux-gnu (64-bit)
# Running under: Scientific Linux release 6.5 (Carbon)
# 
# locale:
#     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
# 
# attached base packages:
#     [1] grid      stats4    parallel  stats     graphics  grDevices utils
# [8] datasets  methods   base
# 
# other attached packages:
#     [1] SDMTools_1.1-221     cba_0.2-15           proxy_0.4-15
# [4] RColorBrewer_1.1-2   gplots_2.17.0        lsr_0.5
# [7] motifStack_1.12.0    Biostrings_2.36.3    XVector_0.8.0
# [10] ade4_1.7-2           MotIV_1.24.0         grImport_0.9-0
# [13] XML_3.98-1.3         magrittr_1.5         readr_0.1.1
# [16] Repitools_1.14.0     GenomicRanges_1.20.5 GenomeInfoDb_1.4.1
# [19] IRanges_2.2.5        S4Vectors_0.6.2      BiocGenerics_0.14.0
# 
# loaded via a namespace (and not attached):
#     [1] Biobase_2.28.0               vsn_3.36.0
# [3] edgeR_3.10.2                 splines_3.2.0
# [5] R.utils_2.1.0                gtools_3.5.0
# [7] Ringo_1.32.0                 R.devices_2.13.0
# [9] affy_1.46.1                  BSgenome_1.36.3
# [11] Rsamtools_1.20.4             rGADEM_2.16.0
# [13] RSQLite_1.0.0                lattice_0.20-33
# [15] limma_3.24.15                digest_0.6.8
# [17] preprocessCore_1.30.0        Matrix_1.2-2
# [19] R.oo_1.19.0                  aroma.affymetrix_2.13.2-9000
# [21] genefilter_1.50.0            zlibbioc_1.14.0
# [23] xtable_1.7-4                 gdata_2.17.0
# [25] aroma.apd_0.6.0              affyio_1.36.0
# [27] BiocParallel_1.2.20          annotate_1.46.1
# [29] aroma.core_2.13.1            survival_2.38-3
# [31] R.methodsS3_1.7.0            R.cache_0.11.0
# [33] MASS_7.3-43                  truncnorm_1.0-7
# [35] R.rsp_0.20.0                 Cairo_1.5-8
# [37] BiocInstaller_1.18.5         tools_3.2.0
# [39] matrixStats_0.14.2           cluster_2.0.3
# [41] AnnotationDbi_1.30.1         lambda.r_1.1.7
# [43] caTools_1.17.1               futile.logger_1.4.1
# [45] RCurl_1.95-4.7               R.huge_0.9.0
# [47] Rsolnp_1.15                  bitops_1.0-6
# [49] base64enc_0.1-3              DNAcopy_1.42.0
# [51] gsmoothr_0.1.7               DBI_0.3.1
# [53] R.filesets_2.8.0             PSCBS_0.44.0
# [55] GenomicAlignments_1.4.1      seqLogo_1.34.0
# [57] rtracklayer_1.28.6           futile.options_1.0.0
# [59] KernSmooth_2.23-15           Rcpp_0.12.1

# linux module loaded

module load apps/gcc/R/3.2.0
module load apps/gcc/BEDTools/2.23.0
module load apps/gcc/homer/2014-02-03


