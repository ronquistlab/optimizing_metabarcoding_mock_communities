#Experiment I
#DADA2 anlysis and visualization.

#loading packages:
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
install.packages("seqRFLP")
library("seqRFLP")
library("reshape2")
library('ggplot2')
library('RColorBrewer')
install.packages("colortools")
library(colortools)



#Starting with making barplots of Set-up Communities:
#Initial community composition:
setwd('~/Desktop/Desktop/3.0_Mock_Buffers/Ready_Trimmed')
colourCount<-13
col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

Design <- read.table("3.0_Design.txt", header=T, sep="\t")

#Starting community composition - ALl communities in one graph:
          l_Design <- melt(Design, id.vars = "X", variable.name = "Phyla")
          fg<-ggplot(l_Design, aes(x = X, y = value, fill = Phyla)) +
            geom_bar(position="fill", stat = "identity") +
            scale_fill_manual(values = col_vector) +
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                               axis.ticks = element_blank(), panel.grid.major=element_blank(), 
                               panel.grid.minor = element_blank()) +
            labs(y = "Number of individuals", x="Community")
          fg
          ggsave(filename = paste0("SetUp_3_0_MockComms.pdf"), device="pdf", width = 4, height = 7,  fg)
          
          d<- (subset(Design, Design$X=="Complex_L"))
          a<-melt(d, variable.name = "Phyla")
          
          gg<- ggplot(a, aes(x = X, y = value, fill = Phyla)) +
            geom_bar(position="fill", stat = "identity") +
            scale_fill_manual(values = col_vector) +
            theme_bw() + theme(axis.text.x = element_blank(), 
                               axis.ticks = element_blank(), panel.grid.major=element_blank(), 
                               panel.grid.minor = element_blank()) +
            labs(y = "Fraction of individuals", x="Community_L")
          ggsave(filename = paste0("SetUp_comm_L.pdf"), device="pdf", width = 3.5 , height = 7,  gg)
          
          
          d<- subset(Design, Design$X=="Medium_M")
          a<-melt(d, variable.name = "Phyla")
          
          gg<- ggplot(a, aes(x = X, y = value, fill = Phyla)) +
            geom_bar(position="fill", stat = "identity") +
            scale_fill_manual(values = col_vector) +
            theme_bw() + theme(axis.text.x = element_blank(), 
                               axis.ticks = element_blank(), panel.grid.major=element_blank(), 
                               panel.grid.minor = element_blank()) +
            labs(y = "Fraction of individuals", x="Community_M")
          ggsave(filename = paste0("SetUp_comm_M.pdf"), device="pdf", width = 3.5 , height = 7,  gg)
          
          
          ggplot(l_Design, aes(x = X, y = value, fill = Phyla)) +
            geom_bar(stat = "identity") + 
            scale_fill_manual(values = getPalette(colourCount)) +
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
          axis.ticks = element_blank(), panel.grid.major=element_blank(), 
          panel.grid.minor = element_blank()) +
            labs(y = "Number of individuals", x="Community") +
            ggtitle("3_0_Mock Set-up")
          
          ggplot(melt(d[subset(Condition2$num, Condition2$Community=="L"),], id.vars = "X", variable.name = "Phyla"),
                 aes(x = X, y = value, fill = Phyla)) + 
            geom_bar(position="fill", stat = "identity") + 
            scale_fill_manual(values = col_vector) +
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
            labs(x = "Sample", y = "Fraction") +
            ggtitle("Mock community L")
          
          ggsave(filename = paste0("SetUp_3_0_MockComm.pdf"), device="pdf", width = 8, height = 5,  fg)
          
          t<-ggplot(l, aes(x = X, y = value, fill = Phyla)) +
            geom_bar(stat = "identity") + 
            scale_fill_manual(values = getPalette(colourCount)) +
            theme_bw() + 
            theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                  panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
            labs(y = "Number of individuals", x="") +
            ggtitle("Mock Set-up", j-1 )
          ggsave(filename = paste0("SetUp_MockComm_", j-1, ".pdf"), device="pdf", width = 3, height = 7,  t)
          j+1 
          }


# DADA2 analysis:
#Sequences trimmed with cutadapt outside of R;
path <- "~/Desktop/Desktop/3.0_Mock_Buffers/Ready_Trimmed" 
list.files(path)

#Make list of your files:
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
fnFs
fnRs

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names


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

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[3:5])  

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filt=filtFs, fnRs, filt.rev=filtRs,
                     maxN=0, maxEE=c(2,4),
                     compress=F, verbose=TRUE)
head(out)
out
out<-as.data.frame(out)
write.csv(out, "COI_counts_perSamples.csv")

filtFs2<-filtFs[-77:-78]
filtFs3<-filtFs2[-49:-52]
filtFs2<-filtFs3[-8]

filtRs2<-filtRs[-77:-78]
filtRs3<-filtRs2[-49:-52]
filtRs2<-filtRs3[-8]

derepF1 <- derepFastq(filtFs2, verbose=TRUE)
derepR1 <- derepFastq(filtRs2, verbose=TRUE)

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
asv <-rownames(t.seqtab.nonchim)
rownames(t.seqtab.nonchim) <- paste("ASV", 1:length(asv), sep="")

### Creating FASTA file with all 482 sequences output of dada2:
ASV<-paste("ASV", 1:length(asv), sep="")
d<-data.frame(ASV, asv)
d.fasta = dataframe2fas(d, file="d.fasta")
write(d.fasta, "499_ASVs.fasta")


############ Getting ASV table with only full, errorless(?maybe/possibly/hopefully) COI seqs.
###########
count.df <- as.data.frame(t.seqtab.nonchim)
count.df$seq.len <- nchar(getSequences(seqtab.nochim))
count.df$Seq <- getSequences(seqtab.nochim)

COI.418 <- subset(count.df, count.df$seq.len==418)
d<-data.frame(rownames(COI.418), COI.418$Seq)
d.fasta = dataframe2fas(d, file="418_COI.fasta")
write.table(COI.418, file="All_418AVSs_Table.txt")

COI.418$Seq <- NULL
COI.418$seq.len <- NULL

COI.412 <- subset(count.df, count.df$seq.len==412)
f<-data.frame(rownames(COI.412), COI.412$Seq)
f.fasta = dataframe2fas(f, file="412_COI.fasta")
write.table(COI.412, file="All_412AVSs_Table.txt")

COI.412$seq.len <- NULL
COI.412$Seq <- NULL

#COI.412 <- COI.412[-3,]

subset(d, d[,1]==rownames(COI.418))
rownames(COI.418)
d[,1]

COI.415 <- subset(count.df, count.df$seq.len==415)
rowSums(COI.415)
COI.349 <- subset(count.df, count.df$seq.len==349)

a<-rowSums(count.df)
Seqs[[26]]
Seqs[[150]]
Seqs[[360]]
Seqs[[10]]
Seqs[[10]]


# Exporting ASVs sequences for 412 length as .fasta. 
COI.412 <- subset(count.df, count.df$seq.len==412)
COI.412$seq.len <- NULL
SEQs <- data.frame(asv.seqs, rownames(COI.412))
COI412.fasta = dataframe2fas(SEQs, file="COI_412.fasta")
write(COI412.fasta, "COI_412_COI_ASV.fasta")


####Track numbers of reads along pipeline:
#Track numbers of reads throught the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF1, getN), sapply(dadaR1, getN), sapply(merger1, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track
write.csv(track, "Tracking_COI_Reads_DADA2.txt")
# Plot tracking of reads as histogram:
ggplot(data=track, aes(x=nonchim, y=input, fill=state, color=input, alpha=input)) +
  geom_bar(stat="identity", position ="identity") +
  scale_colour_manual(values=c("lightblue4","red")) +
  scale_fill_manual(values=c("lightblue","pink")) +
  scale_alpha_manual(values=c(.3, .8))

#Plotting number of reads per sample:
ggplot(data=seq, aes(x=X, y=input, fill=state)) +
  theme_bw() + coord_flip() +
  geom_bar(stat="identity", position ="identity") 

ggplot(Condition, aes(x = fact, y = TotalSeq)) + 
  theme_bw() + geom_bar(stat = "identity", color="grey9", fill="darkolivegreen")  +
  coord_flip() + labs(y = "Number of sequences", x="Sample") +
  theme(axis.text.y = element_text(size=7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(data=track, aes(x=track[,1])) +
  geom_bar(aes(y=input),stat="identity",position ="identity",alpha=.3,fill='lightblue',color='lightblue4') +
  geom_bar(aes(y=nonchim),stat="identity",position ="identity",alpha=.8,fill='pink',color='red')

ggplot(data=my_data,aes(x=Block))+
  geom_bar(aes(y=Start),stat="identity",position ="identity",alpha=.3,fill='lightblue',color='lightblue4') +
  geom_bar(aes(y=End),stat="identity",position ="identity",alpha=.8,fill='pink',color='red')

#######

######
######
######

###Visualisation - Stacked barplots of relative read counts per species:
#Load data:
Counts<-read.csv("Counts_Final_horizontal.txt", header = T,  sep="\t")
rowSums(Counts[,2:23])

Vertical_Counts<-t(Counts)
write.csv(Vertical_Counts, "Counts_Final_Vertical.txt")
d<-read.table("Final_Counts_Vertical.txt", sep="\t", header = T)

Condition <- read.delim("Conditions.txt", sep="\t", header=T, row.names = 2)
Condition$num <- 1:71

#load packages
library("reshape2")
library('ggplot2')
library('RColorBrewer')
install.packages("colortools")
library(colortools)

#Create color vector for species in the mock community (same colors as in Experiment II)
colourCount<-13
col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

Condition2 <- Condition[order(Condition$Buffer, Condition$Community, Condition$Lysis_time),]

Condition[order(Condition$Buffer, Condition$Community, Condition$Lysis_time]
Condition[order(Condition$Buffer, Condition$Community, Condition$Lysis_time),]


# Plotting STacked Barplots:
#i=0
## Try-out Community L
ggplot(melt(d[subset(Condition2$num, Condition2$Community=="L"),], id.vars = "X", variable.name = "Phyla"),
       aes(x = X, y = value, fill = Phyla)) + 
  geom_bar(position="fill", stat = "identity") + 
  scale_fill_manual(values = col_vector) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  labs(x = "Sample", y = "Fraction") +
  ggtitle("Mock community L")

## Try-out Community S and only B1
ggplot(melt(d[subset(Condition2$num, Condition2$Buffer=="a" & Condition2$Community=="S"),], id.vars = "X", variable.name = "Phyla"),
       aes(x = X, y = value, fill = Phyla)) + 
  geom_bar(position="fill", stat = "identity") + 
  scale_fill_manual(values = col_vector) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  labs(x = "Sample", y = "Fraction") +
  ggtitle("Vesterinen_S")

#load additional package:
install.packages("gridExtra")
library(gridExtra)

#Set up parameters:
B<- "a"
C<- "L"
L<- "6"
p3 <- ggplot(melt(d[subset(Condition2$num, Condition2$Buffer==B & Condition2$Community==C & Condition2$Lysis_time==L),], id.vars = "X", variable.name = "Phyla"),
             aes(x = X, y = value, fill = Phyla)) + 
  geom_bar(position="fill", stat = "identity", show.legend = FALSE) + 
  scale_fill_manual(values = col_vector) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  labs(x = "Sample", y = "Fraction") +
  ggtitle(C, "Vesterinen")
grid.arrange(p1, p2, p3, p4, nrow = 1)

# saving to the current working directory.
#Size for crowded graphs
ggsave(filename = paste0("Mock_community_", i, ".pdf"), device="pdf", width = 10, height = 10,  g)
#smaller graphs (i.e. 6 communities)
ggsave(filename = paste0("Mock_community_", i, ".pdf"), device="pdf", width = 8, height = 10,  g)
#small communities (3)
ggsave(filename = paste0("Mock_community_", i, ".pdf"), device="pdf", width = 5, height = 7,  g)

# Plot Barplots with legend on the bottom:
ggplot(melt(d[rownames(subset(Condition, Condition$Community==i)),], id.vars = "X", variable.name = "Phyla"),
       aes(x = X, y = value, fill = Phyla)) + 
  geom_bar(position="fill", stat = "identity") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
                     legend.position="bottom") +
  labs(x = "Sample", y = "Fraction") +
  scale_x_discrete(labels = subset(Condition$Name, Condition$Community==i)) +
  ggtitle("Mock community", i )
ggsave(filename = paste0("Legend_bottom.pdf"), device="pdf", width = 8, height = 5,  g)


##########################3
##########################

#
#
#

#### nMDS plots
### making a plot for all our samples (excluding blanks and mixed samples)
#adding weight and number datat from starting communities

# Working directory:
setwd('~/Desktop/Desktop/3.0_Mock_Buffers/Ready_Trimmed')
#Load packages:
install.packages("vegan")
library("vegan")
install.packages("ggplot2")
library("ggplot2")

#Load data:
d<-read.table("Vertical_Percentages_plus_SETup.txt", sep="\t", header = T)
d1<-d[1:77,]
Condition <- read.delim("Conditions_plus_SETup.txt", sep="\t", header=T, row.names = 2)
Condition$num <- 1:77

#COMM "S" May2022
      # Calculate Bray-Curtis distances and make an nMDS:
      mdl_S <- metaMDS(d1[rownames(subset(Condition, Condition$Community=="S")),], distance = "bray")
      mdl_S$stress
      
      data.scores <- as.data.frame(scores(mdl_S, "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
      data.scores$sample <- rownames(subset(Condition, Condition$Community=="S"))  # create a column of sample names, from the rownames of data.scores
      data.scores$grp1 <- subset(Condition$Community, Condition$Community=="S")
      data.scores$grp2 <- subset(Condition$Buff_name, Condition$Community=="S") #  add the grp variable created earlier
      data.scores$tag <- subset(Condition$Tag, Condition$Community=="S")
      data.scores$lysis <- subset(Condition$Lysis_time, Condition$Community=="S")
      data.scores$grp3 <- factor(paste0(data.scores$grp2, ".", data.scores$grp1))
      data.scores$treat <- subset(Condition$Treat, Condition$Community=="S")
      head(data.scores)  #look at the data
      
      species.scores <- as.data.frame(scores(mdl_S, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
      species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
      head(species.scores)  #look at the data

      #Plot the NMDS - changes according to May2022 meeting - community "S"
      ggplot() +
        geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2, color=factor(grp2), shape=factor(treat), alpha=factor(lysis)), size=4)+
        coord_fixed() + ## need aspect ratio of 1!
        #geom_segment(data = species.scores,
        #            aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
        #           arrow = arrow(length = unit(0.15, "cm")), colour = "grey", alpha=0.7)+
        #geom_text(data=species.scores, aes(x=NMDS1,y=NMDS2,label=species), size=2.5) +  # add the species labels
        scale_alpha_discrete(labels=c("2h", "4h", "6h", "homogenate")) +
        scale_color_manual(values = c("peachpuff4", "black", "black","cornflowerblue"), labels=c( "B2 (Ivanova et al. 2006", "Biomass", "Numbers", "B1 (Vesterinen et al. 2016)")) + 
        scale_shape_manual(values=c(9, 15, 19, 8)) +
        coord_equal() + labs(color="Buffer", alpha="Digestion time") +
        ggtitle("Mock community S") +
        theme_bw()

#COMM "M" May 2022
        mdl_M <- metaMDS(d1[rownames(subset(Condition, Condition$Community=="M")),], distance = "bray")
        mdl_M$stress
        
        data.scores_M <- as.data.frame(scores(mdl_M, "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
        data.scores_M$sample <- rownames(subset(Condition, Condition$Community=="M"))  # create a column of sample names, from the rownames of data.scores
        data.scores_M$grp1 <- subset(Condition$Community, Condition$Community=="M")
        data.scores_M$grp2 <- subset(Condition$Buff_name, Condition$Community=="M") #  add the grp variable created earlier
        data.scores_M$tag <- subset(Condition$Tag, Condition$Community=="M")
        data.scores_M$lysis <- subset(Condition$Lysis_time, Condition$Community=="M")
        data.scores_M$grp3 <- factor(paste0(data.scores$grp2, ".", data.scores$grp1))
        data.scores_M$treat <- subset(Condition$Treat, Condition$Community=="M")
        head(data.scores_M)  #look at the data
        
        species.scores_M <- as.data.frame(scores(mdl_M, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
        species.scores_M$species <- rownames(species.scores_M)  # create a column of species, from the rownames of species.scores
        head(species.scores_M)  #look at the data
        
        #changes according to May2022 meeting - community "S"
        ggplot() +
          geom_point(data=data.scores_M, aes(x=NMDS1,y=NMDS2, color=factor(grp2), shape=factor(treat), alpha=factor(lysis)), size=4)+
          coord_fixed() + ## need aspect ratio of 1!
          #geom_segment(data = species.scores,
          #            aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
          #           arrow = arrow(length = unit(0.15, "cm")), colour = "grey", alpha=0.7)+
          #geom_text(data=species.scores, aes(x=NMDS1,y=NMDS2,label=species), size=2.5) +  # add the species labels
          scale_alpha_discrete(labels=c("2h", "4h", "6h", "homogenate")) +
          scale_color_manual(values = c("peachpuff4", "black", "black","cornflowerblue"), labels=c( "B2 (Ivanova et al. 2006", "Biomass", "Numbers", "B1 (Vesterinen et al. 2016)")) + 
          scale_shape_manual(values=c(9, 15, 19, 8)) +
          coord_equal() + labs(color="Buffer", alpha="Digestion time") +
          ggtitle("Mock community M") +
          theme_bw()

#COMM "L" May 2022
          mdl_L <- metaMDS(d1[rownames(subset(Condition, Condition$Community=="L")),], distance = "bray")
          mdl_L$stress
          
          data.scores_L <- as.data.frame(scores(mdl_L, "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
          data.scores_L$sample <- rownames(subset(Condition, Condition$Community=="L"))  # create a column of sample names, from the rownames of data.scores
          data.scores_L$grp1 <- subset(Condition$Community, Condition$Community=="L")
          data.scores_L$grp2 <- subset(Condition$Buff_name, Condition$Community=="L") #  add the grp variable created earlier
          data.scores_L$tag <- subset(Condition$Tag, Condition$Community=="L")
          data.scores_L$lysis <- subset(Condition$Lysis_time, Condition$Community=="L")
          data.scores_L$grp3 <- factor(paste0(data.scores$grp2, ".", data.scores$grp1))
          data.scores_L$treat <- subset(Condition$Treat, Condition$Community=="L")
          head(data.scores_L)  #look at the data
          
          species.scores_L <- as.data.frame(scores(mdl_L, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
          species.scores_L$species <- rownames(species.scores_L)  # create a column of species, from the rownames of species.scores
          head(species.scores_L)  #look at the data
          
          #changes according to May2022 meeting - community "S"
          ggplot() +
            geom_point(data=data.scores_L, aes(x=NMDS1,y=NMDS2, color=factor(grp2), shape=factor(treat), alpha=factor(lysis)), size=4)+
            coord_fixed() + ## need aspect ratio of 1!
            #geom_segment(data = species.scores,
            #            aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
            #           arrow = arrow(length = unit(0.15, "cm")), colour = "grey", alpha=0.7)+
            #geom_text(data=species.scores, aes(x=NMDS1,y=NMDS2,label=species), size=2.5) +  # add the species labels
            scale_alpha_discrete(labels=c("2h", "4h", "6h", "homogenate")) +
            scale_color_manual(values = c("peachpuff4", "black", "black","cornflowerblue"), labels=c( "B2 (Ivanova et al. 2006", "Biomass", "Numbers", "B1 (Vesterinen et al. 2016)")) + 
            scale_shape_manual(values=c(9, 15, 19, 8)) +
            coord_equal() + labs(color="Buffer", alpha="Digestion time") +
            ggtitle("Mock community L") +
            theme_bw()
          
          
          
          
          
##This script produces heatmaps for Mock comm Experiment I (Buffers)

#load packages:   
library(usethis)
library(devtools)
install.packages('pheatmap')
library(pheatmap)
library(RColorBrewer)
install.packages('gridBase')
library(gridBase)
          
          #Comm S only
          S<- read.csv("~/Desktop/ExpI_Counts_for_heatmap.csv", header=TRUE, sep=",", row.names=1)
          S_relAb <- (S/rowSums(S))
          S_relAb[1:3, 1:3]
          S_filtered <- S_relAb
          #data_filtered[data_filtered < 0.001]
          S_filtered[S_filtered == 0] <- NA #turn zeros to NA to apply custum color
          S_matrix <- as.matrix(t(S_filtered))
          heatmap_palette <- colorRampPalette(c("lightgrey", "yellow", "orange", "red"), interpolate = c("linear"))(n = 100)#Heatmap color palette
          labs_S <- c("2h_B1_A", "4h_B1_A", "6h_B1_A", "Hom_B1_A", "2h_B1_B", "4h_B1_B", "6h_B1_B", "Hom_B1_B", "2h_B1_C", "4h_B1_C", "6h_B1_C", "Hom_B1_C", "2h_B2_A", "4h_B2_A", "6h_B2_A", "Hom_B2_A", "2h_B2_B", "4h_B2_B", "6h_B2_B", "Hom_B2_B", "2h_B2_C", "4h_B2_C", "6h_B2_C", "Hom_B2_C")
          
          pheatmap(S_matrix,
                   border_color = "black",
                   na_col = "white",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(4, 8, 12, 16, 20), #gaps to separate data according to species in my dataset. Can be removed if not necessary. 
                   cellwidth = 16,
                   cellheight = 16,
                   color = heatmap_palette,
                   fontsize = 10,
                   labels_col = labs_S )
          
          
          #comm M only
          M<- read.csv("~/Desktop/ExpI_M_Counts_for_heatmap.csv", header=TRUE, sep=",", row.names=1)
          M_relAb <- (M/rowSums(M))
          M_relAb[1:3, 1:3]
          M_filtered <- M_relAb
          M_filtered[M_filtered == 0] <- NA #turn zeros to NA to apply custum color
          M_matrix <- as.matrix(t(M_filtered))
          heatmap_palette <- colorRampPalette(c("lightgrey", "yellow", "orange", "red"), interpolate = c("linear"))(n = 100)#Heatmap color palette
          labs_M <- c("2h_B1_A", "4h_B1_A", "6h_B1_A", "Hom_B1_A", "2h_B1_B", "4h_B1_B", "6h_B1_B", "Hom_B1_B", "2h_B1_C", "4h_B1_C", "6h_B1_C", "Hom_B1_C", "2h_B2_A", "4h_B2_A", "6h_B2_A", "Hom_B2_A", "2h_B2_B", "4h_B2_B", "6h_B2_B", "Hom_B2_B", "2h_B2_C", "4h_B2_C", "6h_B2_C", "Hom_B2_C")
          
          pheatmap(M_matrix,
                   border_color = "black",
                   na_col = "white",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(4, 8, 12, 16, 20), #gaps to separate data according to species in my dataset. Can be removed if not necessary. 
                   cellwidth = 16,
                   cellheight = 16,
                   color = heatmap_palette,
                   fontsize = 10,
                   labels_col = labs_M)
          
          #comm L only
          L<- read.csv("~/Desktop/ExpI_L_Counts_for_heatmap.csv", header=TRUE, sep=",", row.names=1)
          L_relAb <- (L/rowSums(L))
          L_relAb[1:3, 1:3]
          L_filtered <- L_relAb
          L_filtered[L_filtered == 0] <- NA #turn zeros to NA to apply custum color
          L_matrix <- as.matrix(t(L_filtered))
          heatmap_palette <- colorRampPalette(c("lightgrey", "yellow", "orange", "red"), interpolate = c("linear"))(n = 100)#Heatmap color palette
          labs_L <- c("2h_B1_A", "4h_B1_A", "6h_B1_A", "Hom_B1_A", "2h_B1_B", "4h_B1_B", "6h_B1_B", "2h_B1_C", "4h_B1_C", "6h_B1_C", "Hom_B1_C", "2h_B2_A", "4h_B2_A", "6h_B2_A", "Hom_B2_A", "2h_B2_B", "4h_B2_B", "6h_B2_B", "Hom_B2_B", "2h_B2_C", "4h_B2_C", "6h_B2_C", "Hom_B2_C")
          
          pheatmap(L_matrix,
                   border_color = "black",
                   na_col = "white",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(4, 7, 11, 15, 19), #gaps to separate data according to species in my dataset. Can be removed if not necessary. 
                   cellwidth = 16,
                   cellheight = 16,
                   color = heatmap_palette,
                   fontsize = 10,
                   labels_col = labs_L)
          
          