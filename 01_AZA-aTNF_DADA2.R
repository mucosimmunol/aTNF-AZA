# AZA Mikrobiom - DADA2 Workflow v1 - on EC2 t3.large instance #


library(dada2)
library(stringr)
library(ggplot2)

packageVersion("dada2")
outpath <- "dada2/"

######  ---------------- Filter --------------------------
#Filename parsing
path <- "input/"
filtpath <- file.path(path, "filtered")
fns <- list.files(path, pattern="fastq.gz")
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names
snames <- sapply(str_split(basename(fnFs), "_"), `[`, 1)
sample.names <- sapply(str_split(snames, "-",n = 2 ), '[', 2)

#Inspect read quality
qpF <- plotQualityProfile(fnFs[1:2])
qpR <- plotQualityProfile(fnRs[1:2])

ggsave(paste0(outpath, "qualprofileF.png"), qpF, device="png")
ggsave(paste0(outpath, "qualprofileR.png"), qpR, device="png")


#Filter and Trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(290, 270), maxN=0, maxEE = c(2,3),
                     truncQ = 2, compress=TRUE, multithread = FALSE, verbose=TRUE)

#Discard low depth samples, because they can cause errors in the merge step 
keep <- out[,"reads.out"] > 20 # Or other cutoff
filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq"))[keep]
filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq"))[keep]

saveRDS(out, "dada2/out.rds")



#####----------- Error Rates --------------------
# Learn error rates
set.seed(100)
errF <- learnErrors(filtFs, nbases = 1e8, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE)

## save error models
saveRDS(errF, paste0(outpath, "errF.RDS"))
saveRDS(errR, paste0(outpath, "errR.RDS"))
## Plot errors
perrf <- plotErrors(errF, nominalQ=TRUE)
perrR <- plotErrors(errR, nominalQ=TRUE)
ggsave(filename = paste0(outpath,"errorsF.png"), plot = perrf, device = "png")
ggsave(filename = paste0(outpath,"errorsR.png"), plot = perrR, device = "png")

####---------- Sample inference and merger of paired-end reads -----
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in seq_along(sample.names)) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
x <- dim(seqtab)
print(x)
dist <- table(nchar(getSequences(seqtab)))
print(dist)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(488,539)]
saveRDS(seqtab2, paste0(outpath, "seqtab_cut.rds"))
saveRDS(seqtab, paste0(outpath, "seqtab_all.rds"))


####----- Remove chimera and assign taxonomy-----
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE)
saveRDS(seqtab.nochim, paste0(outpath, "seqtab_nochim.rds"))
dim(seqtab.nochim)
x<-sum(seqtab.nochim)/sum(seqtab)
# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, "silva132/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(tax, paste0(outpath, "tax_final.rds"))


#Track reads through pipeline
getN <- function(x) sum(getUniques(x))

mergers_track <- cbind(1,getN(mergers[[1]]))
for (i in seq(2,185)) {
  sample<-i
  abundance<-getN(mergers[[i]])
  mergers_track <- rbind(mergers_track, cbind(sample,abundance))
}

track <- cbind(out, mergers_track[,2], rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "tabled", "nonchim")
track <- as.data.frame(track)
track$chimrate <- track$nonchim/track$tabled*100
rownames(track) <- sample.names
head(track)
write.table(track, file = "track_reads.csv", append=F, quote=T, sep=";", row.names = T, col.names = T)  


saveRDS(ddF, paste0(outpath, "ddF.rds"))
saveRDS(ddR, paste0(outpath, "ddR.rds"))
saveRDS(mergers, paste0(outpath, "mergers.rds"))
rm(list = ls())

## For species assignment, change instance type to r5.2xlarge, restart instance and run the following snippet
library(dada2)
setwd("18077/")
outpath <- "dada2/"
tax <- readRDS(paste0(outpath, "tax_final.rds"))
spec <- addSpecies(tax, "silva132/silva_species_assignment_v132.fa.gz", verbose=TRUE)
saveRDS(spec, paste0(outpath, "species.rds"))