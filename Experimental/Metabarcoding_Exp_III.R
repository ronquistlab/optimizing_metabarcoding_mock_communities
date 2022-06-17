library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
install.packages("seqRFLP")
library("seqRFLP")

setwd('~/Documents/PROJECTS/Exp_III')
path <- "~/Documents/PROJECTS/Exp_III" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Make list of your files:
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
fnFs
# TRIM PRIMERS - COI - BF3, BR2:
# Input your primer sequences:
FWD <- "CCHGAYATRGCHTTYCCHCG"  
REV <- "TCDGGRTGNCCRAARAAYCA"  

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#Remove sequences with Ns, cause they cause problems in matching sequences:
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
# now fitered files are in fnFs.filtN and fnRs.filtN directories.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[5]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[5]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[5]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[5]]))

cutadapt <- "/Users/elzbiwas/Library/Python/3.8/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#Set up arguments for cutadapt:
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 


#Run cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "--discard-untrimmed", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


#######Proceeeding with trimmed reads:
# Reading in files again. Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])  


filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filt=filtFs, cutRs, filt.rev=filtRs,
                     minLen=180, 
                     maxN=0, maxEE=2,
                     compress=TRUE, verbose=TRUE)
head(out)

derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)

#Learn Error rates
errF <- learnErrors(derepF1, multithread=FALSE)
errR <- learnErrors(derepR1, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Infer sample composition - denoising
dadaF1 <- dada(derepF1, err=errF, multithread=FALSE)
dadaR1 <- dada(derepR1, err=errR, multithread=FALSE)
print(dadaF1)
print(dadaR1)

#Merge forward and reverse reads
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose=TRUE)

#Remove chimeras:
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)

seqtab <- makeSequenceTable(merger1.nochim)
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
#Inspect dimensions of the asv table and distribution of seqs
#dim(seqtab)
#table(nchar(getSequences(seqtab)))
#getSequences(seqtab)
#nchar(getSequences(seqtab))

dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
t.seqtab.nonchim <- t(seqtab.nochim)
asv <-rownames(t.seqtab.nonchim)
rownames(t.seqtab.nonchim) <- paste("ASV", 1:length(asv), sep="")