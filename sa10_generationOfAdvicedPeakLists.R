setwd("myDirectory/")
advices <- read.table("CurratedPeakListAdvices.txt", header = T, sep = "\t") # manually generated
load("Rdata/listOfAllRep_57.RData")

# merge
merge <- subset(advices, advices$Adviced.peak.list == "Merge peak lists")[, 1]

writeNeatMergeCommands <- function(repName){
    md <- listOfAllRep[[repName]][,"name"]
    filenames <- paste0("sortedTFtables/", md, sep = ".narrowPeak.s")
    filenames <- paste0(filenames, collapse = " ")
    commands <- rbind(
        paste0("cat ", filenames, " > advicedPeakLists/", repName, ".cat.bed"),
        paste0("sortBed -i advicedPeakLists/", repName, ".cat.bed", " > advicedPeakLists/", repName, ".cat.sorted.bed"),
        paste0("bedtools merge -c 7 -o sum -i advicedPeakLists/", repName, ".cat.sorted.bed > advicedPeakLists/", repName, ".bed")
    )
    return(commands)
}

mergeCommands <- do.call(rbind, lapply(merge, writeNeatMergeCommands))
write.table(mergeCommands, file = "mergeReplicatesForFEBS.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
# run this

# intersect
# all intersect are only two files, thanks god
intersect <- subset(advices, advices$Adviced.peak.list == "Intersect peak lists")[,1]

setwd("TFtables/")
library(GenomicRanges)

writeNiceIntersectedFiles <- function(repName){
    md <- listOfAllRep[[repName]][,"name"]
    filenames <- paste0("sortedTFtables/", md, sep = ".narrowPeak.s")
    beds<-list(
        read.table(filenames[1]),
        read.table(filenames[2])
    )
    GR <- lapply(beds, function(x) GRanges(x[,1], IRanges(x[,2], x[,3])))
    mtch <- as.matrix(findOverlaps(GR[[1]], GR[[2]]))
    # Be wary of index nightmare and nested expressions
    tab<-cbind(as.character(beds[[1]][mtch[,1],1]),
               do.call(pmax, as.data.frame(cbind(beds[[1]][mtch[,1],2], beds[[2]][mtch[,2],2]))),
               do.call(pmin, as.data.frame(cbind(beds[[1]][mtch[,1],3], beds[[2]][mtch[,2],3]))),
               rowMeans(cbind(beds[[1]][mtch[,1],7], beds[[2]][mtch[,2],7]))
               )
    write.table(tab, file = paste0("advicedPeakLists/", repName, ".bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    message(paste(repName, "done!"))
}

empty <- lapply(intersect, writeNiceIntersectedFiles)

# copy one file
copies <- c("SydhH1hescCmycIggrab", "HaibK562Hdac2sc6296V0416102", "SydhK562P300Iggrab", "SydhK562Pol2s2", "HaibSknshNrsfV0416101","SydhK562Atf3")
names(copies) <- c("H1-hESC_None_c-Myc","K562_None_HDAC2","K562_None_p300","K562_None_Pol2_phosphoS2","SK-N-SH_None_NRSF","K562_None_ATF3")

writeNiceCopiesOfFiles <- function(i){
    sfn <- copies[i]
    table <- read.table(paste0("sortedTFtables/wgEncodeAwgTfbs", sfn, "UniPk.narrowPeak.s"))
    write.table(table[,c(1,2,3,7)], file = paste0("advicedPeakLists/", names(copies)[i], ".bed"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    message(names(copies)[i])
}
empty <- lapply(1:length(copies), writeNiceCopiesOfFiles)

# K562_None_Pol2
repName <- "K562_None_Pol2"
filenames <- as.character(listOfAllRep[["K562_None_Pol2"]][,"name"])[c(1,3,4)]
filenames <- paste0("sortedTFtables/", filenames, sep = ".narrowPeak.s")
filenames <- paste0(filenames, collapse = " ")
paste0("cat ", filenames, " > advicedPeakLists/", repName, ".cat.bed")
paste0("sortBed -i advicedPeakLists/", repName, ".cat.bed", " > advicedPeakLists/", repName, ".cat.sorted.bed")
paste0("bedtools merge -c 7 -o sum -i advicedPeakLists/", repName, ".cat.sorted.bed > advicedPeakLists/", repName, ".bed")

# GM12878_None_p300
repName <- "GM12878_None_p300"
filenames <- as.character(listOfAllRep[[repName]][,"name"])

NiceIntersectedFiles <- function(pairs){
    md <- pairs
    filenames <- paste0("sortedTFtables/", md, sep = ".narrowPeak.s")
    beds <- list(
        read.table(filenames[1]),
        read.table(filenames[2])
    )
    GR <- lapply(beds, function(x) GRanges(x[,1], IRanges(x[,2], x[,3])))
    mtch <- as.matrix(findOverlaps(GR[[1]], GR[[2]]))
    # Be wary of index nightmare and nested expressions
    tab <- cbind(as.character(beds[[1]][mtch[,1],1]),
               do.call(pmax, as.data.frame(cbind(beds[[1]][mtch[,1],2], beds[[2]][mtch[,2],2]))),
               do.call(pmin, as.data.frame(cbind(beds[[1]][mtch[,1],3], beds[[2]][mtch[,2],3]))),
               rowMeans(cbind(beds[[1]][mtch[,1],7], beds[[2]][mtch[,2],7]))
    )
    return(tab)
}
tab1 <- NiceIntersectedFiles(filenames[1:2])
tab2 <- NiceIntersectedFiles(filenames[2:3])
tab3 <- NiceIntersectedFiles(filenames[c(1,3)])

tab <- rbind(tab1, tab2, tab3)

write.table(tab, file=paste0("advicedPeakLists/", repName, ".unmerged.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
paste0("sortBed -i advicedPeakLists/", repName, ".unmerged.bed", " > advicedPeakLists/", repName, ".cat.sorted.bed")
paste0("bedtools merge -c 4 -o sum -i advicedPeakLists/", repName, ".cat.sorted.bed > advicedPeakLists/", repName, ".bed")
