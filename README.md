# MSc Project 2021
# MSc_Project_Aishwarya_2021

### R
```R
getwd()
setwd("/home/user01/Paul/A_project/fastq")

library("Rsubread")
#index directory
#("/home/user01/Suriya/Book/ene/")
list.files()
#index
#setwd("/home/user01/Paul/AlignmentFiles/genome/")
#buildindex(basename = "Homo", reference = "Homo_sapiens.GRCh38.dna.toplevel.fa", memory = 1600000)

#read file
reads1 <- list.files( path = "/home/user01/Paul/A_project/fastq", pattern = "*_1.fastq.gz" )
reads2 <- list.files( path = "/home/user01/Paul/A_project/fastq", pattern = "*_2.fastq.gz" )
all.equal(length(reads1),length(reads2))
align(index = "/home/user01/Suriya/Ref/Homo", readfile1=reads1, readfile2=reads2, 
      input_format="FASTQ", output_format="BAM", nthreads=32)
########################################################################################

#align(index = "/home/user01/Suriya/Ref/Homo", "SRR9610797_1.fastq.gz", "SRR9610797_2.fastq.gz", 
#      input_format="FASTQ", output_format="BAM", nthreads=32)

#library("Rsamtools")
library("Rsamtools")
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
library("GenomicFeatures")
#######################################################################
setwd("~/Suriya/Ref")
gtffile <- ("GRCh38.99.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
txdb
ebg <- exonsBy(txdb, by="gene")
ebg
library("GenomicAlignments")
library("BiocParallel")
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
setwd("/home/user01/Paul/A_project/fastq")
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
se
head(assay(se))
write.csv(assay(se), "Count.csv")

```
