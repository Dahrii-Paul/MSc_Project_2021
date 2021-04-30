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
