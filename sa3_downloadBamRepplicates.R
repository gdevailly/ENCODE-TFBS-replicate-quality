# This step describe the process of downloading and merging the .bam files
setwd("myDirectory/")
# A(nother) manualy build table.
bamMD <- read.table("replicatesTable_bam_handCurrated_new.txt", header = TRUE)
sum(bamMD$bamSize) # 192 Go

summary(bamMD$bamSize)

nrow(bamMD)

commands <- as.data.frame(with(bamMD, paste0("wget ", bam)))
write.table(commands, file = "TFtables/bamReps/download.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
# execute this and wait...

nrow(bamMD)
bamMD <- bamMD[grep("wgEncode", bamMD$name),]
nrow(bamMD)

bamMD$bamFileName <- with(bamMD, substr(as.character(bam), nchar(as.character(bam)) - 14, nchar(as.character(bam))))

length(levels(factor(as.character(bamMD$allRep))))
length(levels(factor(as.character(bamMD$name))))
bamMD$allRep <- as.character(bamMD$allRep)
bamMD[which(bamMD$allRep == "K562_None_Pol2(phosphoS2)"), "allRep"] <- "K562_None_Pol2_phosphoS2"
subrep <- levels(factor(bamMD$name))
length(subrep)

writeMergeBamCommands <- function(expName){
    tempT <- subset(bamMD, bamMD$name == expName)
    bamN <- as.character(tempT[1, "bamFileName"])
    if(nrow(tempT) > 1) {
        for (i in 2:nrow(tempT)) {
            bamN <- paste(bamN, tempT[i, "bamFileName"])
        }
    }
    commands <- as.data.frame(c(
        paste0("samtools merge ", expName, ".bam ", bamN),
        paste0("rm ", bamN)
    ))
    return(commands)
}

listOfCommands <- lapply(subrep, writeMergeBamCommands)

tableOfCommands <- listOfCommands[[1]]
for(i in 2:length(listOfCommands)) {
    tableOfCommands <- rbind(tableOfCommands, listOfCommands[[i]])
}

write.table(tableOfCommands, file = "TFtables/bamReps/mergeBam.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
# run this, require samtools
