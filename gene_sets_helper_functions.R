initialize_conv_df <- function(){
  
  #Create an object from biomaRt that translates gene symbols between mouse and human, returns a df 
  
  require(biomaRt)
  require(tidyverse)
  
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  conv <- getLDS(attributes="entrezgene", attributesL="entrezgene", mart=human, martL =mouse)
  colnames(conv) <- c("human", "mouse")
  conv <- apply(conv, 2, as.character)
  conv <- as.data.frame(conv)
  conv
}

convert_entrez_to_human <- function(gsea_vec, conv = conv, isges=F, convFrom="mouse"){
  
  #Convert murine ENTREZIDs to human gene symbols, requires conv object
  
  require(tidyverse)
  
  if (isges==T){
    df <- data.frame(ID=names(gsea_vec), Log2=gsea_vec)
  } else {
    df <- data.frame(ID = gsea_vec)
  }
  
  j <- inner_join(df, conv, by=c("ID" = convFrom))
  j
}

to_entrezid <- function(symbolvector, ODB = org.Hs.eg.db){
  
  #Quickly convert a vector of SYMBOL to ENTREZID
  
  AnnotationDB::select(ODB, columns = "ENTREZID", keys = symbolvector, keytype = "SYMBOL")$ENTREZID
}

gene_list_from_cluster <- function(clusters, species_convert = T, top = 10){
  
  #Retrieve SYMBOLS from genes represented in gene set results
  
  df <- as.data.frame(clusters)
  if (top > 0){
    df <- df[0:top,]
  }
  
  outp_list <- lapply(df$geneID, function(x) strsplit(x, "/")[[1]])
  names(outp_list) <- df$Description
  if (species_convert){
    outp_list <- lapply(outp_list, function(x) convert_entrez_to_human(x, conv, convFrom = "human"))
    outp_list <- lapply(outp_list, function(x) as.vector(x$mouse))
  }
  outp_list <- lapply(outp_list, function(x) bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db))
  outp_list <- lapply(outp_list, function(x) as.vector(x$SYMBOL))
  
  
}

mm_symbols_to_hs_entrez <- function(symbols){
  require(clusterProfiler)
  require(org.Mm.eg.db)
  outp <- convert_entrez_to_human(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID, conv)$human
  outp
}


plot_heatmap <- function(vst_or_rlog, genes, main = ""){
  
  #plot heatmap from genes SYMBOLS, tests that genes are in assay
  
  require(pheatmap)
  
  if (!all(genes %in% rownames(plot_matrix))){warning("Not all requested genes are in the experiment")}
  plot_matrix <- assay(vst_or_rlog)
  genes <- genes[genes %in% rownames(plot_matrix)]
  pheatmap(plot_matrix[genes,], main = main)
}