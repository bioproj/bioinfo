.libPaths(c("/home/wy/workspace/transcriptomics/transcriptomics/R_LIBS",.libPaths()))
.libPaths()
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
options()$repos  
options()$BioC_mirror 

BiocManager::version()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install("Rsamtools")
BiocManager::install("DESeq2")



setwd("transcriptomics/")
library(tximport)
library(readr)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("./testData/chr22_genes.gtf", "gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
head(tx2gene)



files = list.files(path="out/salmon",pattern = "quant.sf",recursive = T,full.names = T)
names(files) <- c("HBR_Rep1","HBR_Rep2","HBR_Rep3","UHR_Rep1","UHR_Rep2","UHR_Rep3")
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)



# library(GenomicFeatures)
library(readr)

head(txi.salmon$counts)




