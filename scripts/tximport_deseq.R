library(tximportData)
library(tximport)
library(DESeq2)
library(stringr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gplots)

args <- commandArgs()
print(args)
wd <- args[6]
id <- args[7]
#v_s <- args[8]

#Locate the samples.txt made according to the given template
samples <- read.table(file.path(wd, "Samples.txt"), header=TRUE)
samples$sample

#Path to salmon quant files & Naming
files <- file.path(wd, samples$sample, "quant.sf")
files
names(files) <- paste0(samples$sample)
all(file.exists(files))

#tximport the transcripts counted by salmon (txOut is specifically for transcripts, if gene level summarization is required, please go for tx2gene)
txi <- tximport(files, type= "salmon", txOut=TRUE, countsFromAbundance="scaledTPM")

#tximport is a necessary offset when quantifying transcript, therefore unnecessary during gene D.E analyses.

#Starting D.E analysis with DESeq2
sampleTable <- data.frame(condition = factor(samples$condition))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds

# Remove redundant rows
dds <- dds[rowSums(counts(dds)) > 10,]
dds <-DESeq(dds)
nrow(dds)

if (nrow(dds) > 1000){
vst <- vst(dds, blind=FALSE)
} else {
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
}
#rownames(dds)

#VST table definition
#vst <- args[8](dds, blind=FALSE)


# Plot Principal Component Analysis (PCA)
plotPCA(vst, intgroup="condition", ntop=nrow(counts(dds)))


# Plot Boxplot
a <- DESeq2::plotPCA(vst, intgroup="condition")
a + geom_label(aes(label = samples$sample),)
nudge <- position_nudge(y = 1)
a + geom_label(aes(label = samples$sample), position = nudge)
a + geom_text(aes(label = samples$sample), position = nudge, size=3 )
boxplot(assay(vst), col= c("Red", "Red", "Red", "Green", "Green", "Green"), pch=".",
        vertical=TRUE, cex.axis=0.5, main = "LncRNA_test",
        las=2, ylab="assay(vst)", xlab="Samples", ylim=c(-10,30),
        font.main= 5, font.axis=0.5, font.lab=2 )


# Plot correlation heatmap
cU <-cor( as.matrix(assay(vst)))
cols <- c("dodgerblue3", "firebrick3")[samples$sample]
heatmap.2(cU, symm=TRUE, col= colorRampPalette(c("darkblue","white"))(100),
          labCol=colnames(cU), labRow=colnames(cU),
          distfun=function(c) as.dist(1 - c),
          trace="none",
          Colv=TRUE, cexRow=0.9, cexCol=0.9, key=F,
          font=2,
          RowSideColors=cols, ColSideColors=cols)

#heatmap(test, col= colorRampPalette(c("blue","yellow","white"))(100), font=2)

# Plot dispersion plot
plotDispEsts(dds)

# parse condition representation for contrast from samples file
df <- samples[!duplicated(samples[c('condition')]), ]
df1 <- as.character(df$condition)

# define filename and condition representation for output file generation
res <- results(dds, contrast=c("condition",df1))
summary(res)
grp.mean <- sapply(levels(dds$condition),
                   function(lvl)
                     rowMeans(counts(dds,normalized=TRUE)[,dds$condition== lvl]))
norm.counts <- counts(dds, normalized=TRUE)
all <- data.frame(res, assay(vst))
csv <- paste(id, ".csv")
write.table(all, file=csv, sep=",")

