#Experiment II
#DADA2 anlysis and visualization.

#loading packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
a
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
install.packages("seqRFLP")
library("seqRFLP")

#Set working directory:
setwd('~/Desktop/Desktop/1.0.Mock_Re-Seq/Re-seq_Mock1')
path <- "~/Desktop/Desktop/1.0.Mock_Re-Seq/Re-seq_Mock1" 
list.files(path)

#Make list of your files:
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#COI primers - find and count occurences:
      #COI-mito - BF3-BR2
      FWD <- "CCHGAYATRGCHTTYCCHCG"  #BF3
      REV <- "TCDGGRTGNCCRAARAAYCA"  #BR2
      
      Orients <- function(primer) {
        # Create all orientations of the input sequence
        require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                     RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
      }
      FWD.orients <- Orients(FWD)
      REV.orients <- Orients(REV)
      
      primerHits <- function(primer, fn) {
        # Counts number of reads in which the primer is found
        nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
      }
      
      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[22]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[22]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[22]]), 
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[22]]))

#Trim reads with cutadapt program (called from R but is an external feature):
      cutadapt <- "/Users/elzbiwas/miniconda3/bin/cutadapt" 
      system2(cutadapt, args = "--version") # Run shell commands from R
      
      path.cut <- file.path(path, "cutadapt2")
      if(!dir.exists(path.cut)) dir.create(path.cut)
      fnFs.cut <- file.path(path.cut, basename(fnFs))
      fnRs.cut <- file.path(path.cut, basename(fnRs))
      fnFs.cut
      
      R1.flags <- paste("-g", FWD)
      R2.flags <- paste("-G", REV)
      
      #Run cutadapt
      for(i in seq_along(fnFs)) {
        system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "--discard-untrimmed", # -n 2 required to remove FWD and REV from reads
                                   "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                                   fnFs[i], fnRs[i])) # input files
      }
      
      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[22]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[22]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[22]]), 
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[22]]))


#######Proceeding with trimmed reads:
# Reading in files again. Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[4:5])
plotQualityProfile(cutRs[4:5])  

filtFs <- file.path(path.cut, "filtered2", basename(cutFs))
filtRs <- file.path(path.cut, "filtered2", basename(cutRs))

out <- filterAndTrim(cutFs, filt=filtFs, cutRs, filt.rev=filtRs,
                     truncLen=c(225, 210), minLen=180, 
                     maxN=0, maxEE=2,
                     compress=TRUE, verbose=TRUE)
head(out)
out<-as.data.frame(out)
write.csv(out, "COI_counts_perSamples.csv")
out$num<- 1:length(out$reads.in)
out
filtFs1<-filtFs[-105:-114]
filtRs1<-filtRs[-105:-114]

#Dereplicate:
derepF1 <- derepFastq(filtFs1, verbose=TRUE)
derepR1 <- derepFastq(filtRs1, verbose=TRUE)

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
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
t.seqtab.nonchim <- t(seqtab.nochim)

count.df <- as.data.frame(t.seqtab.nonchim)
count.df$seq.len <- nchar(getSequences(seqtab.nochim))
count.df$seq.len

COI <- subset(count.df, count.df$seq.len==418 | count.df$seq.len==415)


### Creating FASTA file with all 1173 COI (length 415 and 418) sequences output of dada2:
Seqs<-rownames(COI)
ASV<-paste("ASV", 1:length(Seqs), sep="")
d<-data.frame(ASV, Seqs)
d.fasta = dataframe2fas(d, file="COI_ASVs_Sequences.fasta")
write(d.fasta, "COI_ASVs_Sequences.fasta")



####__________________________________________________________________________###
###Visualisation Experiment II

#Load packages:
library("reshape2")
library('ggplot2')
library('RColorBrewer')
install.packages("colortools")
library(colortools)

#working directory to get data:
setwd('~/Desktop/Desktop/1.0.Mock_Re-Seq/')
d<-read.table("Vertical_Counts.txt", sep="\t", header = T)

Condition <- read.delim("Condition_Exp1.txt", sep="\t")
Design <- read.csv("Design_Mock_Comm.txt", header=T, sep="\t")
Design_weighted <- read.csv("Design_weighted.txt", header=T, sep="\t")

#Make a color palette (custom color for each species, same as in Exp I)
palette<-c("#DBE8D3", "#D4E853", "#81778B", "#DF92A4", "#DC5BA0", "#D2A78B", "#C46FD7", "#DED997", "#E4BF52",
           "#AAD7DF", "#87A2E1", "#6DDD85", "#E048DC", "#DD5759", "#D8A2E1", "#68DEDB", "#DE8A4D", "#82EABB",
           "#7FEB50", "#7EA47F", "#853BE6", "#7174DD", "#DBC9E0", "#64B7DB", "#BAE28A")

#######
#Starting community composition:
  #Set working directory for where you'll produce graps:
  setwd('~/Desktop/Desktop/1.0.Mock_Re-Seq/Re-seq_Mock1/COI_graphs/')

            # Set-up in terms of NUMBERS:
            # Set "j" to a number of community you want to plot (1:5)
            j=5
            l <- melt(Design[j,], id.vars = "X", variable.name = "Phyla")
            
            t<-ggplot(l, aes(x = X, y=value, fill = Phyla)) + 
              geom_bar(position="fill", stat = "identity") + 
              scale_fill_manual(values = palette) +
              theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
              labs(y = "Fraction", x="") +
              ggtitle("Numbers")
            t
            ggsave(filename = paste0("SetUp_COI_", j, ".pdf"), device="pdf", width = 5.5, height = 9,  t)
            
            # Set-up in terms of BIOMASS:
            Design_weighted <- read.csv("Design_weighted.txt", header=T, sep="\t")
            j=5
            lw <- melt(Design_weighted[j,], id.vars = "X", variable.name = "Phyla")
            
            t<-ggplot(lw, aes(x = X, y=value, fill = Phyla)) + 
              geom_bar(position="fill", stat = "identity") + 
              scale_fill_manual(values = palette) +
              theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
              labs(y = "Fraction", x="") +
              ggtitle("Weight")
            t
            ggsave(filename = paste0("WEIGHT_COI_MockComm_", j, ".pdf"), device="pdf", width = 5.5, height = 9,  t)

#########
# PLotting barplots of relative read abundances:
            #set "i" to the community you want to plot (1:5); 
            i=5
            ## Try-out Community 0
            g<-ggplot(melt(d[rownames(subset(Condition, Condition$Comm==i)),], id.vars = "X", variable.name = "Phyla"),
                      aes(x = X, y = value, fill = Phyla)) + 
              geom_bar(position="fill", stat = "identity") + 
              scale_fill_manual(values = palette) +
              theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12), legend.position="bottom") +
              labs(x = "Sample", y = "Fraction") +
              scale_x_discrete(labels = subset(Condition$Label2, Condition$Comm==i)) +
              ggtitle("Mock community", i )
            g
            # saving to current working directory.
            #Size for crowded graphs
            ggsave(filename = paste0("COI_Mock_community_", i, ".pdf"), device="pdf", width = 9, height = 9,  g)
            
            
            # Plot with legend on the bottom:
            ggplot(melt(d[rownames(subset(Condition, Condition$Community==i)),], id.vars = "X", variable.name = "Phyla"),
                   aes(x = X, y = value, fill = Phyla)) + 
              geom_bar(position="fill", stat = "identity") + 
              scale_fill_manual(values = palette) +
              theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
                                 legend.position="bottom") +
              labs(x = "Sample", y = "Fraction") +
              scale_x_discrete(labels = subset(Condition$Name, Condition$Community==i)) +
              ggtitle("Mock community", i )
            ggsave(filename = paste0("Legend_bottom.pdf"), device="pdf", width = 8, height = 5,  g)


#####Plotting only homogenates  (or other lysis time) for all replicates in a given community:
i=5
## Try-out Community 
g<- ggplot(melt(d[rownames(subset(Condition, Condition$Comm==i & Condition$Dig_time=="Homog")),], id.vars = "X", variable.name = "Phyla"),
           aes(x = X, y = value, fill = Phyla)) + 
  geom_bar(position="fill", stat = "identity") + 
  scale_fill_manual(values = palette) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12), legend.position="bottom") +
  labs(x = "Sample", y = "Fraction") +
  scale_x_discrete(labels = subset(Condition$Sample, Condition$Comm==i & Condition$Dig_time=="Homog")) +
  ggtitle("Mock community", i)
# saving to current working directory.
#Size for crowded graphs
ggsave(filename = paste0("Homogenates_COI_", i, ".pdf"), device="pdf", width = 3.5, height = 9,  g)


#### Plot number of reads Initial and then COI
#
Counts<-read.delim("Counts_Tracking.txt", sep="\t")
### Plotting intial and final number of reads per sample:
# Overlaying bar plot
require(scales)
ggplot(data=Counts, aes(x=name), size=2) +
  coord_flip() + theme_bw() + labs(y = "Number of sequences", x="Sample") +
  geom_bar(aes(y=All),stat="identity",position ="identity",alpha=.3,fill='lightblue',color='lightblue4') +
  geom_bar(aes(y=COI),stat="identity",position ="identity",alpha=.3,fill='grey30',color='lightblue4') +
  geom_bar(aes(y=COI_final),stat="identity",position ="identity",alpha=.8,fill='pink',color='violetred4') +
  theme(axis.text.y = element_text(size=7)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = comma)


#### nMDS - Calculating Bray-Curtis distances and making an nMDS:
              ### An nMDS for all samples plus original community compositions (Numbers and Weights)
              setwd('~/Desktop/Desktop/1.0.Mock_Re-Seq/')
              d2<-read.table("Vertical_Counts_NUM_WEI.txt", sep="\t", header = T)
              Condition2 <- read.delim("Condition_ExpI_NUM_WEI.txt", sep="\t")
              ### making a plot for all our samples (excluding blanks and mixed samples)
              #load packages:
              install.packages("vegan")
              library("vegan")
              
              #prepare files:
              rownames(d2)<- d2$X
              d2$X<-NULL
              rownames(Condition2)<- Condition2$Sample
              
              #Calculate nMDS:
              MDS2 <- metaMDS(d2, distance = "bray")
              MDS2$stress
              MDS2
              plot(MDS2)
              par(mfrow = c(1, 1))
              plot(MDS2, type = "t")
              stressplot(MDS2)
              
              
              data.scores <- as.data.frame(scores(MDS2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
              data.scores$sample <- Condition2$Sample  # create a column of sample names, from the rownames of data.scores
              data.scores$grp1 <- Condition2$Comm_Rome
              data.scores$grp2 <- Condition2$Dig_time2
              data.scores$grp3 <- Condition2$DigTime3 #  add the grp variable created earlier 
              data.scores$Type <- Condition2$Type
              data.scores$tag <- Condition2$Rep
              head(data.scores)  #look at the data
              
              species.scores <- as.data.frame(scores(MDS2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
              species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
              head(species.scores)  #look at the data
              
              library(RColorBrewer)
              # Define the number of colors you want
              nb.cols <- 17
              mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
              # Create a ggplot with 18 colors 
              
              ggplot() + 
                geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2, colour=factor(grp1), shape=factor(Type), alpha=factor(grp3)), size=3) +
                scale_shape_manual(values=c(15,19, 8, 9)) +
                scale_alpha_discrete(labels = c("1h", "2h", "4h", "6h", "8h", "20h", "Homogenised/SetUp")) +
                geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size=2.5, alpha=0.6) +  # add the species labels
                scale_fill_manual(values = wes_palette("IsleofDogs1")) + 
                coord_equal() + labs(color="Community", shape="Treatment", alpha="Lysis Time") +
                theme_bw() 
              
              #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size=2.5, alpha=0.6) +  # add the species labels
                