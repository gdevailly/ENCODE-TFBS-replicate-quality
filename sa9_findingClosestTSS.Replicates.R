# first we write bed tools comands to find closest TSS for each peak file
setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")
setwd("TFtables/mergedDuplicates/")
files <- list.files()
merged <- files[grep(".merged.bed$", files)]

TSStable <- "../../gencode.transcriptTSS.hg19.s.bed"
command2 <- vector(mode="character")
for(i in c(6, 18, 42)) {
    command2[i] <- paste0("bedtools closest -a ", merged[i], " -b ", TSStable, " -D b -t first ",
                        "> ", merged[i], ".closestTSS" )
}
command2 <- as.data.frame(command2[c(6, 18, 42)])

write.table(command2[,], file = "findClosestTSS_57.sh", sep = "\t", quote=F, row.names=F, col.names=F)
# we run those command and we start a new R session
###

setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")

loadAnnot <- function(name) {
    tempT <- read.table(paste0("TFtables/mergedDuplicates/", name, ".merged.bed.closestTSS"))
    colnames(tempT) <- c("chr", "start", "end", "score","TSS_chr", "TSS_start", "TSS_end", "ENSG", "whatever", "TSS_strand", "ENST", "type",
                       "name", "yet_another_name", "yet_another_score", "dist")
    return(tempT)
}

a <- Sys.time()
closestTSS <- lapply(names(listOfAllRep), loadAnnot)
Sys.time() - a
names(closestTSS) <- names(listOfAllRep)
save(closestTSS, file = "Rdata/closestTSS.replicates_57.RData")
# this is used just bellow
# we start a new session
###

setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")
load("Rdata/binaryRepMatrix_57.RData")
defs <- c("distal" = 10000, "proximal" = 1000)
load("Rdata/closestTSS.replicates_57.RData")
categ <- read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

# we format data for easy ploting
prepareData <- function(name){
    Nrep <- nrow(listOfAllRep[[name]])
    binMat <- expMatrix[[name]][,5:(Nrep+4)]
    closest <- closestTSS[[name]]
    counts <- matrix(nrow = 7, ncol = Nrep+1)
    rownames(counts) <- c("far_up", "mid_up", "close_up", "at", "close_down", "mid_down", "far_down")
    colnames(counts) <- c("common_peaks", as.character(listOfAllRep[[name]][,"name"]))
    
    commons <- which(rowSums(binMat) == Nrep)
    counts["far_up", "common_peaks"] <- length(which(closest[commons, "dist"] > defs["distal"]))
    counts["mid_up", "common_peaks"] <- length(which(closest[commons, "dist"] <= defs["distal"] &
                                                       closest[commons, "dist"] > defs["proximal"]))
    counts["close_up", "common_peaks"] <- length(which(closest[commons, "dist"] <= defs["proximal"] &
                                                         closest[commons, "dist"] > 0))
    counts["at", "common_peaks"] <- length(which(closest[commons, "dist"] == 0))
    counts["close_down", "common_peaks"] <- length(which(closest[commons, "dist"] >= -defs["proximal"] &
                                                           closest[commons, "dist"] < 0))
    counts["mid_down", "common_peaks"] <- length(which(closest[commons, "dist"] >= -defs["distal"] &
                                                         closest[commons, "dist"] < -defs["proximal"]))
    counts["far_down", "common_peaks"] <- length(which(closest[commons, "dist"] < -defs["distal"]))
    
    for (i in 1:Nrep) {
        peaks <- which(binMat[,i] == 1)
        counts["far_up", i+1] <- length(which(closest[peaks, "dist"] > defs["distal"]))
        counts["mid_up", i+1] <- length(which(closest[peaks, "dist"] <= defs["distal"] &
                                                closest[peaks, "dist"] > defs["proximal"]))
        counts["close_up", i+1] <- length(which(closest[peaks, "dist"] <= defs["proximal"] &
                                                  closest[peaks, "dist"] > 0))
        counts["at", i+1] <- length(which(closest[peaks, "dist"] == 0))
        counts["close_down", i+1] <- length(which(closest[peaks, "dist"] >= -defs["proximal"] &
                                                    closest[peaks, "dist"] < 0))
        counts["mid_down", i+1] <- length(which(closest[peaks, "dist"] >= -defs["distal"] &
                                                  closest[peaks, "dist"] < -defs["proximal"]))
        counts["far_down", i+1] <- length(which(closest[peaks, "dist"] < -defs["distal"]))
    }
    normCounts <- counts
    for (i in 1:(Nrep+1)) {
        normCounts[,i] <- counts[,i]/sum(counts[,i])
    }
    return(normCounts)
}

a <- Sys.time()
listOfPeakPostions <- lapply(names(listOfAllRep), prepareData)
Sys.time() - a # 3 s
names(listOfPeakPostions) <- names(listOfAllRep)

save(listOfPeakPostions, file = "Rdata/TSSproximity.replicates_57.RData")




