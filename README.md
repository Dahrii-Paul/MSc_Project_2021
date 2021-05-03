# MSc Project 2021

### R code from bam file extract count file
```R
#Get the working directory
getwd()
setwd("/home/user01/Paul/A_project/fastq")
#List the file in the directory
list.files()

library("Rsamtools")
#load phenotype sample table
sampleTable <- read.csv("sample_hiNPCs.CSV", row.names = 1)

list.files()
sampleTable <- read.csv("sample_hiNPCs.CSV", row.names = 1)
sampleTable
filenames <- paste0(sampleTable$SampleName, ".bam")
sampleTable
filenames
filenames <- paste0(sampleTable$SampleName, ".BAM")
file.exists(filenames)
filenames
#####################################################
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
library("GenomicFeatures")
#######################################################################
setwd("~/Suriya/Ref")
#Annotation  
gtffile <- ("GRCh38.99.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
txdb
ebg <- exonsBy(txdb, by="gene")
ebg
library("GenomicAlignments")

#Distribute the core 
library("BiocParallel")
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
setwd("/home/user01/Paul/A_project/fastq")

#Extract the count from bam files
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
se
head(assay(se))
write.csv(assay(se), "Count.csv")

```
# Extract DEGS from count file and find LncRNA
```R
getwd()
setwd("/Users/paul/Documents/work/MSc_Project/R_")
list.files()

library(DESeq2)
library(ggplot2)

#load count table
count_Data<-read.csv(file = 'Count.csv', header = TRUE, sep = ",", row.names = 1)
count_Data
#load the metaData
metaData <- read.csv('meta_zika.csv', header = TRUE, sep = ",", row.names = 1)
metaData
#%in% operator used to identify if an element belongs to a vector
all(colnames(count_Data) %in% rownames(metaData))
all(colnames(count_Data) == rownames(metaData))
#DEGs
dds <- DESeqDataSetFromMatrix(countData = count_Data,
                              colData = metaData,
                              design = ~ Condition)
head(dds)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
### factor analysis ####
dds$Condition <- factor(dds$Condition, levels = c("PE","FSS","Control"))
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds$Condition <- droplevels(dds$Condition)
### Deseq ###
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
res
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
#check the plot
plotCounts(dds, gene="ENSG00000000003", intgroup="Condition")
#PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Condition")
#############################################################################
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), 
                                    pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
###########################################################################
sum(res05$padj < 0.05, na.rm=TRUE)
df <- subset(res, res$pvalue <= 0.05 & (res$log2FoldChange > 2 | res$log2FoldChange < -2))
write.csv(df, "UpandDown.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange > 2)
write.csv(df, "Upregulated.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange < -2)
write.csv(df, "downregulated.csv")
### normalization ###
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 10)
write.csv(assay(vsd), "vst_norm.csv")
write.csv(assay(rld), "rlog_norm.csv")
## to extract count data ##
DGE.results.sorted <- read.csv("UpandDown.csv", row.names = 1)
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
head(DGE.results.sorted)
normalizeddata <- read.csv("vst_norm.csv", row.names = 1)
Visualization <- normalizeddata[DGEgenes,] %>% data.frame()
write.csv(Visualization, "UpandDown_Visualization.csv")

### Finding Gene Symbol ##
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- read.csv("UpandDown.csv")
head(df)
genes <- df$X

G_list3 <- getBM(filters= "ensembl_gene_id", 
                 attributes= c("ensembl_gene_id",
                               "hgnc_symbol",
                               "gene_biotype",
                               "ensembl_transcript_id",
                               'transcript_length',
                               'transcript_start',
                               'transcript_end'), 
                 values=genes, mart= mart)
head(G_list3)
head(df)
colnames(df)
#Rename df 1st column using index
names(df)[1]<-"ensembl_gene_id"
A <- merge(df, G_list3, by = "ensembl_gene_id")
write.csv(A, "UpandDown_withsymbol.csv")

```
# [Filter out lncRNA file]()
