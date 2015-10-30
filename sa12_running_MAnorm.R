setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")
library(Repitools)
library(magrittr)
library(readr)
testRep <- c("HepG2_None_JunD", "K562_None_Rad21", "GM12878_None_YY1",
             "HeLa-S3_None_c-Myc", "GM12878_None_PAX5", "GM12891_None_Pol2")

prepareForMAnorm <- function(repName) {
    metaD <- listOfAllRep[[repName]]
    
    myBeds <- metaD$name
    myBeds <- lapply(myBeds, function(x) read_tsv(
        paste0("TFtables/sortedTFtables/", x, ".narrowPeak.s"),
        col_names = FALSE))
    names(myBeds) <- metaD$name
    
    myBams <- metaD$name
    myBams <- lapply(myBams, function(x) BAM2GRanges(
        paste0("TFtables/bamReps/", x, ".bam")
        ))
    myBams <- lapply(myBams, as.data.frame)
    names(myBams) <- metaD$name
    
    system(paste0("mkdir TFtables/MAnorm/", repName))
    
    empty <- lapply(as.character(metaD$name), function(x) {
        write.table(
            myBeds[[x]][,1:3],
            file = paste0("TFtables/MAnorm/", repName, "/", x, ".peaks.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
            )
        write.table(
            myBams[[x]][,c(1:3, 5)],
            file = paste0("TFtables/MAnorm/", repName, "/", x, ".reads.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
        )
    })
    
    system(paste0("cp TFtables/MAnorm/MAnorm.sh TFtables/MAnorm/", repName, "/MAnorm.sh"))
    system(paste0("cp TFtables/MAnorm/MAnorm.r TFtables/MAnorm/", repName, "/MAnorm.r"))
    
    message(paste(repName, "done !"))
    
    return(paste0(
        "./MAnorm.sh ",
        metaD$name[1], ".peaks.bed ", metaD$name[2], ".peaks.bed ",
        metaD$name[1], ".reads.bed ", metaD$name[2], ".reads.bed ",
        "150 150"
        ))
}

manormCommands <- lapply(testRep, prepareForMAnorm)

###
# execute these


./MAnorm.sh wgEncodeAwgTfbsHaibHepg2JundPcr1xUniPk.peaks.bed wgEncodeAwgTfbsSydhHepg2JundIggrabUniPk.peaks.bed wgEncodeAwgTfbsHaibHepg2JundPcr1xUniPk.reads.bed wgEncodeAwgTfbsSydhHepg2JundIggrabUniPk.reads.bed 150 150

./MAnorm.sh wgEncodeAwgTfbsHaibK562Rad21V0416102UniPk.peaks.bed wgEncodeAwgTfbsSydhK562Rad21UniPk.peaks.bed wgEncodeAwgTfbsHaibK562Rad21V0416102UniPk.reads.bed wgEncodeAwgTfbsSydhK562Rad21UniPk.reads.bed 150 150

./MAnorm.sh wgEncodeAwgTfbsHaibGm12878Yy1sc281Pcr1xUniPk.peaks.bed wgEncodeAwgTfbsSydhGm12878Yy1UniPk.peaks.bed wgEncodeAwgTfbsHaibGm12878Yy1sc281Pcr1xUniPk.reads.bed wgEncodeAwgTfbsSydhGm12878Yy1UniPk.reads.bed 150 150

./MAnorm.sh wgEncodeAwgTfbsSydhHelas3CmycUniPk.peaks.bed wgEncodeAwgTfbsUtaHelas3CmycUniPk.peaks.bed wgEncodeAwgTfbsSydhHelas3CmycUniPk.reads.bed wgEncodeAwgTfbsUtaHelas3CmycUniPk.reads.bed 150 150

./MAnorm.sh wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk.peaks.bed wgEncodeAwgTfbsHaibGm12878Pax5n19Pcr1xUniPk.peaks.bed wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk.reads.bed wgEncodeAwgTfbsHaibGm12878Pax5n19Pcr1xUniPk.reads.bed 150 150

./MAnorm.sh wgEncodeAwgTfbsHaibGm12891Pol2Pcr1xUniPk.peaks.bed wgEncodeAwgTfbsSydhGm12891Pol2IggmusUniPk.peaks.bed wgEncodeAwgTfbsHaibGm12891Pol2Pcr1xUniPk.reads.bed wgEncodeAwgTfbsSydhGm12891Pol2IggmusUniPk.reads.bed 150 150
