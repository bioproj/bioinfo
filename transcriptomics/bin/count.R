#!/usr/local/bin/Rscript
# .libPaths(c("/ssd2/springCloudProd/pipeline/R_LIBS",.libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if(F){
  setwd("./out/")
  getwd()
  gtf_ <- "../genome/Homo_sapiens.GRCh38.113.gtf.gz"
  quant_path_ <- "salmon"
  treatment_group_ <- "SH1_1.clean,SH1_2.clean"
  control_group_ <- "TRC1.clean,TRC2.clean"
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
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

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
library(rtracklayer)
message("start------------")
gtf_data = rtracklayer::import(gtf_) 
message("gtf_data------------")
gtf_data <- gtf_data |> as.data.frame()
tx2gene <- gtf_data |>
  filter(grepl("transcript",type)) |>
  mutate(TXNAME=paste0(transcript_id,".",transcript_version),GENEID=paste0(gene_id,".",gene_version)) |>
  dplyr::select(TXNAME,GENEID)

id_map <- gtf_data |>
  filter(grepl("gene",type)) |>
  mutate(TXNAME=paste0(transcript_id,".",transcript_version),GENEID=paste0(gene_id,".",gene_version)) |>
  dplyr::select(gene_id=GENEID,gene_name,gene_length=width)

message("------------")

# txdb <- makeTxDbFromGFF(gtf_)
# k <- keys(txdb, keytype = "GENEID")
# tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = c("TXNAME","TXID","EXONNAME","CDSID") )|>
#   dplyr::select("TXNAME","GENEID")
# head(tx2gene)
# tx2gene |>
#   filter(TXNAME %in% c("DQ459412"))
# keytypes(txdb)
# columns(txdb)

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

if(!is.na(control_group_) && !is.na(treatment_group_)){
  vs_dir <- paste0(control_group_,"_VS_",treatment_group_)
  dir.create(vs_dir)
  setwd(vs_dir)
}

txi.salmon$counts |>
  as.data.frame()  |>
  rownames_to_column("gene_id") |>
  remove_rownames() |>
  inner_join(id_map,by="gene_id") |>
  relocate(gene_name,gene_length,.after="gene_id") |>
  write_tsv( file = "count.tsv")


# length(intersect(id_map$gene_id,rownames(txi.salmon$counts)))

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
  
  # as.data.frame(res) |>
  #   rownames_to_column("gene_id") |>
  #   inner_join(id_map,by="gene_id") |>
  #   relocate(gene_name,gene_length,.after="gene_id") |>
  #   write_tsv(file="deg.tsv")
  
  # a <- summary(res) |> print(file = "a.txt")
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
  names(resdata)[1] <- "gene_id"
  resdata <- resdata |>
    inner_join(id_map,by="gene_id") |>
    filter(!is.na(gene_name)) |>
    relocate(gene_name,gene_length,.after="gene_id") 
  resdata |>
    write_tsv(file="deg.tsv")
  
  
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
  
  
  # select <- order(rowMeans(counts(dds,normalized=TRUE)),
  #                 decreasing=TRUE)[1:100]
  deg_df <- resdata |>
    dplyr::filter(padj<0.01) |>
     arrange(desc(padj)) 
  
  df <- as.data.frame(colData(dds)[,c("condition")])
  ntd <- normTransform(dds)
  # vsd <- vst(dds, blind=FALSE)
  # rld <- rlog(dds, blind=FALSE)
  
  # assay(ntd)[deg_df$gene_id,] |>
  #   rownames_to_column(gene_id) 
  
  exp_normTransform <- assay(ntd) |> 
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    inner_join(dplyr::select(id_map,c("gene_name","gene_id")),by="gene_id") |>
    dplyr::select(-gene_id) %>%
     {.[!duplicated(.$gene_name),]} |>
    filter(!is.na(gene_name)) |>
    remove_rownames() |>
    column_to_rownames("gene_name")
  
  # c <- exp_normTransform[deg_df$gene_name,]
  # c <- exp_normTransform[deg_df$gene_name,]
  # 
  # deg_df$gene_name[1:100]
    
  # assay(ntd)[deg_df$gene_id[1:3],]
  # assay(vsd)[deg_df$gene_id[1:3],]
  # assay(rld)[deg_df$gene_id[1:3],]
  
  png(filename = "pheatmap_gene.png", width = 7, height =7, res = 300, units = "in")
  pheatmap(exp_normTransform[deg_df$gene_name[1:50],],scale = "row",  cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE)
  dev.off()
  
  
  
  entrez_id <- bitr(deg_df$gene_name, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)%>%
    pull("ENTREZID")%>%unique()
  message("start clusterProfiler",length(entrez_id))
  
  kegg<- clusterProfiler::enrichKEGG(gene =entrez_id,organism = 'hsa',pAdjustMethod="none")
  kegg <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
  kegg@result |>
    as.data.frame() |>
    rownames_to_column("name")  |>
    write_tsv("KEGG.tsv")
  # png(filename = "KEGG.png",width = 7,height = 6,res = 300,units = "in")
  dotplot(kegg,showCategory=20)
  # dev.off()
  ggsave(filename = "KEGG.png",width = 10,height =9,units = "in")
  
  ego<-enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = 'all',readable = T ,pAdjustMethod="BH")
  ego@result |>
    as.data.frame() |>
    rownames_to_column("name")  |>
    write_tsv("GO.tsv")
  # png(filename = "GO.png",width = 10,height = 9,res = 300,units = "in")
  colorSel="qvalue"
  dotplot(ego,showCategory = 10,split="ONTOLOGY",orderBy = "GeneRatio", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
  # dev.off()
  ggsave(filename = "GO.png",width = 10,height =9,units = "in")
  
  message("stop clusterProfiler")
  
  
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



