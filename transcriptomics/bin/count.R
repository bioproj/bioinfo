#!/usr/local/bin/Rscript
# .libPaths(c("/home/wy/workspace/transcriptomics/transcriptomics/R_LIBS",.libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if(F){
  setwd("./out/")
  getwd()
  gtf_ <- "../testData/chr22_with_ERCC92.gtf.gz"
  quant_path_ <- "salmon"
  treatment_group_ <- "UHR_Rep1_ERCC-Mix1,UHR_Rep2_ERCC-Mix1,UHR_Rep3_ERCC-Mix1"
  control_group_ <- "HBR_Rep1_ERCC-Mix2,HBR_Rep2_ERCC-Mix2,HBR_Rep3_ERCC-Mix2"
  treatment_group_name_ <- "UHR"
  control_group_name_<- "HBR"
}
gtf_ <- args[1]
quant_path_ <- args[2]
treatment_group_ <- args[3]
control_group_ <- args[4]
treatment_group_name_ <- args[5]
control_group_name_ <- args[6]


# Rscript bin/count.R  ./testData/chr22_genes.gtf out/salmon
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)

txdb <- makeTxDbFromGFF(gtf_)
k <- keys(txdb, keytype = "GENEID")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME") |>
  dplyr::select("TXNAME","GENEID")
# head(tx2gene)
# tx2gene |>
#   filter(TXNAME %in% c("DQ459412"))


if(!is.na(control_group_) && !is.na(treatment_group_)){
  treatment_group <- str_split(treatment_group_,",")[[1]]
  control_group <- str_split(control_group_,",")[[1]]
  
  files <<- file.path( c(treatment_group,control_group), "quant.sf")
  files_exists <- all(sapply(files, function(x) file.exists(x)))
  if(!files_exists){
    stop("文件不存在:",paste0(files," "))
  }
}else{
  files <<- list.files(path=quant_path_,pattern = "quant.sf",recursive = T,full.names = T)
}
files_name <- sapply(files, function(x)str_split(x,"/")[[1]][length(str_split(x,"/")[[1]])-1])
names(files) <- files_name
# names(files) <- c("HBR_Rep1","HBR_Rep2","HBR_Rep3","UHR_Rep1","UHR_Rep2","UHR_Rep3")
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
# dim(txi.salmon$counts)
# head(txi.salmon$counts)
txi.salmon$counts |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  write_tsv( file = "count.tsv")


if(!is.na(control_group_) && !is.na(treatment_group_)){
  treatment_group <- str_split(treatment_group_,",")[[1]]
  control_group <- str_split(control_group_,",")[[1]]
  treatment_group_name <- treatment_group_name_
  if(is.na(treatment_group_name_)){
    treatment_group_name <- "treatment"
  }
  control_group_name <- control_group_name_
  if(is.na(control_group_name_)){
    control_group_name <- "control"
  }
  # sample  condition
  # HBR_Rep1    HBR
  # HBR_Rep2    HBR
  # HBR_Rep3    HBR
  # UHR_Rep1    UHR
  # UHR_Rep2    UHR
  # UHR_Rep3    UHR
  samples <- tibble(
    sample = c(treatment_group,control_group),
    condition =   c(rep(treatment_group_name,length(treatment_group)) , rep(control_group_name,length(control_group)) )
  )
  intersect_num <- intersect(colnames(txi.salmon$counts), samples$sample)
  if(length(intersect_num) !=length(colnames(txi.salmon$counts))) {
    stop("样本组长度不一致:\n输入的样本名称:",paste0(samples$sample,collapse = " "),"\n表格的样本名称:",paste0(colnames(txi.salmon$counts),collapse = " "))
  }
  if(!identical(colnames(txi.salmon$counts), samples$sample)){
    stop("样本组位置不一致:\n输入的样本名称:",paste0(samples$sample,collapse = " "),"\n表格的样本名称:",paste0(colnames(txi.salmon$counts),collapse = " "))
  }
  # txi.salmon[,samples$sample]
  
  dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  summary(res) 
  png(filename = "dispersion_plot.png", width = 7, height = 7, res = 300, units = "in")
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  rld <- rlogTransformation(dds)
  # head(assay(rld))

  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(samples$condition))])
  sampleDists <- as.matrix(dist(t(assay(rld))))
  png(filename = "heatmap.png", width = 8, height = 8, res = 300, units = "in")
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[samples$condition],
            RowSideColors=mycols[samples$condition],
            margin=c(30, 20), main="Sample Distance Matrix")
  dev.off()
  
  DESeq2::plotPCA(rld, intgroup="condition")
  ggsave(filename = "pca.png",width = 7,height = 5,units = "in")
  
  # table(res$padj<0.05)
  res <- res[order(res$padj), ]
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  # head(resdata)
  
  png(filename = "pvalue_hist.png", width = 7, height =7, res = 300, units = "in")
  hist(res$pvalue, breaks=50, col="grey")
  dev.off()
  
  png(filename = "ma.png", width = 7, height =7, res = 300, units = "in")
  DESeq2::plotMA(dds, ylim=c(-1,1))
  dev.off()
  
  # Volcano plot
  png(filename = "volcano.png", width = 7, height =7, res = 300, units = "in")
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
  with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  dev.off()
  
  if(F){
    res$symbol <- mapIds(org.Hs.eg.db,
                         keys=row.names(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    res$entrez <- mapIds(org.Hs.eg.db,
                         keys=row.names(res),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
    res$name <- mapIds(org.Hs.eg.db,
                       keys=row.names(res),
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
    # library(gageData)
    head(res)
    data(kegg.sets.hs)
    data(sigmet.idx.hs)
    kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
    head(kegg.sets.hs, 3)
    
    
    
    foldchanges <- res$log2FoldChange
    names(foldchanges) <- res$entrez
    keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
    lapply(keggres, head)
    
    
    # Get the pathways
    keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
      tbl_df() %>%
      filter(row_number()<=5) %>%
      .$id %>%
      as.character()
    keggrespathways
    
    # Get the IDs.
    keggresids <- substr(keggrespathways, start=1, stop=8)
    keggresids
    
    
    # Define plotting function for applying later
    plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
    
    # Unload dplyr since it conflicts with the next line
    detach("package:dplyr", unload=T)
    tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
    
    
  }
    
}



