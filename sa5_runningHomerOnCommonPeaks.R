setwd("myDirectory/")
library(GenomicRanges)
library(magrittr)
load("Rdata/binaryRepMatrix_57.RData")
discrepencies <- c("K562_None_ATF3", "H1-hESC_None_CHD1", "H1-hESC_None_c-Myc",  "K562_None_HDAC2",
                 "K562_None_p300", "GM12878_None_PAX5", "K562_None_Pol2", "K562_None_Pol2_phosphoS2")
names(discrepencies) <- discrepencies

writeCommonPeaksIntBedFile <- function(repName) {
    bed <- expMatrix[[repName]]
    dim(bed)
    bed <- subset(bed, rowSums(bed[5:ncol(bed)]) == ncol(bed) - 4)
    dim(bed)
    write.table(bed[,1:3], file = paste0("TFtables/bedReps/",repName, "_common.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
empty <- lapply(discrepencies, writeCommonPeaksIntBedFile)

prefixInput <- "TFtables/bedReps/"
prefixOutput <- "homerMotifOutput/repCommonPeaks/"
makeHomerCommands <- function(repName) {
    paste0("findMotifsGenome.pl ", prefixInput, repName,"_common.bed hg19 ", prefixOutput, repName,
           "/ -size given -p 20 -S 20 -preparsedDir tempHomer")
}
myCommands <- sapply(discrepencies, makeHomerCommands)
write.table(as.data.frame(myCommands), file = "homerOnCommonPeaksRepplicates.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)




