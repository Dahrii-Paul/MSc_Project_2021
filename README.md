# MSc Project 2021

### R code
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
# DEGs

```R
getwd()
library(biomaRt)
listEnsembl()
ensembl = useEnsembl(biomart="ensembl")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") 
##########################################
#load CSV file
getwd()
setwd("E:\\work\\MSc_MTech\\Project")
#read all the ensemble_ID
data_ID <-read.csv(file = 'Count.csv')
head(data_ID)
class(data_ID)
dim(data_ID)
myvec<-as.vector(data_ID$ID)
class(myvec)
data_attribute<-getBM(attributes=c('ensembl_gene_id', 'gene_biotype'),
                      filters = 'ensembl_gene_id',
                      values = myvec,
                      mart = ensembl)      
write.csv(data_attribute,"data.csv")


lincrna <- getBM(attributes= c('ensembl_gene_id',
                               'gene_biotype',
                               'ensembl_transcript_id',
                               'transcript_length',
                               'chromosome_name',
                               'transcript_start',
                               'transcript_end',
                               'ensembl_exon_id',
                               'exon_chrom_start',
                               'exon_chrom_end', 
                               'strand'), filters = 'ensembl_gene_id',
                 values = myvec, mart = ensembl) 
write.csv(lincrna,"lincrna3.csv")

#extract only lncRNA from lincRNA data frame
#lincrna_gene_biotype <- lincrna$gene_biotype=="lncRNA"
######################################################
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
#load count table
countData<-read.csv(file = 'Count.csv', header = TRUE, sep = ",", row.names = 1)
#load the metaData
metaData <- read.csv('meta_zika.csv', header = TRUE, sep = ",", row.names = 1)
metaData
###
all(colnames(countData) %in% rownames(metaData))
all(colnames(countData) == rownames(metaData))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ Condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
### factor analysis ####

dds$Condition <- factor(dds$Condition, levels = c("PE","FSS","Control"))
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds$Condition <- droplevels(dds$Condition)
### Deseq ###

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
write.csv(res,"result.csv")
#########################################################################
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
df <- subset(res, res$pvalue <= 0.05 & (res$log2FoldChange > 2 | res$log2FoldChange < -2))
write.csv(df, "UpandDown.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange > 2)
write.csv(df, "Upregulated.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange < -2)
write.csv(df, "Downregulated.csv")

### Finding Gene Symbol ##

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- read.csv("UpandDown.csv")
head(df)
genes <- df$ensembl_gene_id
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"), 
                values=genes, mart= mart)
G_list
A <- merge(df, G_list, by = "ensembl_gene_id")
write.csv(A, "UpandDown_withsymbol.csv")
#########################################################################################
#Design matrix for three samples control as the base point with refernce to other 2 sample
#########################################################################################
head(data_ID)
metaData
all(colnames(countData) %in% rownames(metaData))
all(colnames(countData) == rownames(metaData))
library(edgeR)
library(limma)
library(DESeq2)
#normalizing and filtering
#limma
class(countData)
dge <- DGEList(counts=countData)
#remove rows that consistently have zero or very low counts
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=2)
class(metaData)
#####################################################################
#create shorter descriptive levels and labels
sample.bam<-c("SRR9610797.BAM","SRR9610798.BAM","SRR9610799.BAM",
              "SRR9610800.BAM","SRR9610801.BAM","SRR9610802.BAM",
              "SRR9610803.BAM","SRR9610804.BAM","SRR9610805.BAM")
bam<-metaData$Condition
design <- model.matrix(~factor(bam))
colnames(design) <-c("control","FSS","PE")
#filter remove zero
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
#DEGs
logCPM <- cpm(dge, log=TRUE, prior.count=2)
head(logCPM)
fit <- lmFit(logCPM, design)
head(fit)
fit <- eBayes(fit, trend=TRUE)
a<-topTable(fit, coef=ncol(design))
###########################################################
countData<-read.csv(file = 'Count.csv', header = TRUE, sep = ",", row.names = 1)
## Calculate the Counts Per Million measure
myCPM <- cpm(countData)
## Identify genes with at least 0.5 cpm in at least 2 samples
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countData[keep,]
## Convert to an edgeR object
dgeObj <- DGEList(counts.keep)
## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)
## Obtain corrected sample information
sampleinfo <- read.csv('meta_zika.csv', header = TRUE, sep = ",", row.names = 1)
group <- paste(sampleinfo$Condition,sampleinfo$rep,sep=".")
#Create the design matrix
#The main assumption here is that the effect of the status is the same in all type of cells

```















