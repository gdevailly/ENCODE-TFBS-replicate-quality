setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
library(readr)
library(magrittr)
load("Rdata/listOfAllRep_57.RData")
load("Rdata/binaryRepMatrix_57.RData")
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

testRep <- c("HepG2_None_JunD", "K562_None_Rad21", "GM12878_None_YY1",
             "HeLa-S3_None_c-Myc", "GM12878_None_PAX5", "GM12891_None_Pol2")

formatDataFor <- function(repName){
    MAnormResults <- read_tsv(paste0("TFtables/MAnorm/", repName, "/MAnorm_result_commonPeak_merged.xls"))
    colnames(MAnormResults) <- c(colnames(MAnormResults)[-9], "minus_log10_pval")
    sorted_pval <- MAnormResults$minus_log10_pval[order(MAnormResults$minus_log10_pval)]
    
    red_sorted_pval <- rep(0, 200)
    red_sorted_pval[(200 - round(length(which(sorted_pval >= 2))*200/length(sorted_pval))):200] <- 1
    red_sorted_pval[(200 - round(length(which(sorted_pval >= 3))*200/length(sorted_pval))):200] <- 2
    red_sorted_pval[(200 - round(length(which(sorted_pval >= 5))*200/length(sorted_pval))):200] <- 3
    
    labels <- c(0, length(which(sorted_pval < 2)), length(which(sorted_pval < 3)),
                length(which(sorted_pval < 5)), length(sorted_pval))
    at <- c(0, (which(red_sorted_pval == 1)[1] -1) /200, (which(red_sorted_pval == 2)[1] -1)/200,
            (which(red_sorted_pval == 3)[1] -1)/200, 1)
    
    return(list("red_sorted_pval" = red_sorted_pval, "labels" = labels, "at" = at))
}

MAnorm_results <- lapply(testRep, formatDataFor)
names(MAnorm_results) <- testRep
save(MAnorm_results, file = "Rdata/MAnorm_results.RData")
