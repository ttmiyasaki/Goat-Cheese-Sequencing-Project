path <- "/share/magalab/Tyler/Cheese/Demux_CutAdapt3_Trimmed/Trimmed" 
list.files(path) 
#put files into path

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
#sort files into F/R

sample.names <- sapply(strsplit(basename(fnFs), "_Trim_"), `[`, 1)
#seperate sample names

png(file="qualityF1.png")
plotQualityProfile(fnFs[1:4])
dev.off()
#generate plot of Forward read quality 
png(file="qualityR1.png")
plotQualityProfile(fnRs[1:4])
dev.off()
#generate plot of Reverse read quality

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,215),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE	
keep<-out[,"reads.out"]>20
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]
sample.names <-sample.names[keep]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE	
#filter reads

errF <- learnErrors(filtFs, nbases=1e+8 ,multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e+8, multithread=TRUE)
#learn error rate

png(file="errorrate1F.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()	
png(file="errorrate1R.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()
# plot error rate, numbers correspond to quality #s 

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#dereplicate

names(derepFs) <- sample.names
names(derepRs) <- sample.names
#give dereps names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#core algorythim
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,verbose=TRUE)
#merge reads
head(mergers[[1]])
# Inspect the merger data.frame from the first sample

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#make sequence table from 

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#remove chimera

sum(seqtab.nochim)/sum(seqtab)
#check chimera prevalance

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

trackdt = as.data.frame(track)
sum(trackdt$nonchim)/sum(trackdt$input)


taxa <- assignTaxonomy(seqtab.nochim, "/share/magalab/Tyler/Cheese/Demux_CutAdapt2_Trimmed2/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.txt", sep="\t", quote=F, col.names=NA)

  # tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.txt", sep="\t", quote=F, col.names=NA)

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

samples.out <- rownames(seqtab.nochim)
box <- as.integer(substr(samples.out,2,2))
cheese <- substr(samples.out,5,5)
week <- substr(samples.out,8,8)
isCheese <- substr(samples.out,1,4)
milkID<- substr(samples.out,6,6)
samdf <- data.frame(Box=box, Cheese=cheese, Week=week, Identity=isCheese, MilkID=milkID)
samdf$When[samdf$Week==1]<-"Week 1"
samdf$When[samdf$Week==2]<-"Week 2"
samdf$When[samdf$Week==3]<-"Week 3"
samdf$When[samdf$Week==4]<-"Week 4"
samdf$When[samdf$Week==5]<-"Week 5"
samdf$When[samdf$Week==6]<-"Week 6"
samdf$When[samdf$Week==7]<-"Week 7"
samdf$When[samdf$Week==8]<-"Week 8"
samdf$When[samdf$Week=="C"]<-"Core Sample"
samdf$When[samdf$Identity=="Milk"]<-"Week 0 Milk"
samdf$Spike<- "No Spike"
samdf$Spike[samdf$Box>3] <- "Spiked"
samdf$Spike[samdf$MilkID=="S"] <- "Spiked"
samdf$SamType<-"Goat Cheese"
samdf$SamType[samdf$Identity=="Milk"] <- "Goat Milk"
samdf$SamType[samdf$Identity=="Pink"]<-"Pink Contamination"
samdf$MilkType[samdf$Box==1|samdf$Box==4]<-"Raw WT"
samdf$MilkType[samdf$Box==2|samdf$Box==5]<-"WT + LZ"
samdf$MilkType[samdf$Box==3|samdf$Box==6]<-"Raw hLZ"
rownames(samdf) <- samples.out
#head(samdf)
samdf
write.table(samdf, "sample_dataframe.txt", sep="\t", quote=F, col.names=NA)

ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(asv_tax))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
#construct phyloseq and trim mock?


png(file="AlphaDivByWeek10.png")
plot_richness(ps, x="Week", measures=c("Shannon", "Simpson"), color="When")
dev.off()

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
#ordinate data

png(file="BrayNMDS10.png")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
dev.off()

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
png(file="Top20Taxa_10F.png")
plot_bar(ps.top20, x="Week", fill="Family") + facet_wrap(~When, scales="free_x")
dev.off()



