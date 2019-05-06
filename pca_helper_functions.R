get_pca <- function(counts, top_n = 500, flip = F){
  
  #Get PCA from count matrix by using genes with 10 n variance, automatically centers and scales
  
  row_ix <- head(order(rowVars(counts), decreasing = T), n = top_n)
  pca_matrix <- counts[row_ix,]
  if (flip){
    pca_matrix <- t(pca_matrix)
  }
  pca <- prcomp(t(pca_matrix), center = T, scale. = T)
  pca
}

pca_plot <- function(pca, coldata, x = PC1, y = PC2){
  
  #Make PCA scatter plot, given output from get_pca. Argument coldata is sampledata table as used by DESEQ2, expected column named cond
  
  x <- enquo(x)
  y <- enquo(y)
  plot_df <- pca$x %>% as.data.frame() %>% mutate(name = rownames(.))
  coldata <- coldata %>% as.data.frame() %>% mutate(name = rownames(.))
  plot_df <- plot_df %>% inner_join(coldata, by = "name")
  plot_df %>% ggplot(aes_string(x = !! x, y = !! y, col = "cond"), size = 10) + geom_point() + labs(title = paste0("PCA plot for ", !! x, " and ", !! y)) + facet_wrap(~cond) + theme_bw()
}

plot_loadings <- function(pca, PC = PC1){
  
  #Plot loadings from a PCA object. Orders bars from lowest to highest 
  
  PC <- enquo(PC)
  pca$rotation %>% as.data.frame() %>% mutate(gene = rownames(.))  %>% ggplot(aes(x = reorder(gene, !! PC), y =!! PC)) + geom_col() + theme_bw()
}

get_loadings <- function(pca, PC = PC1){
  
  #Get the loadings from a PCA object. By default gets the rows part of the highest 20% of the entire range
  
  PC <- enquo(PC)
  pca <- pca$rotation %>% as.data.frame() %>% mutate(gene = rownames(.)) %>% dplyr::select(gene, !! PC) %>% arrange(desc(abs(!! PC)))
  pca %>% ggplot(aes(x = reorder(gene, !! PC), y = !! PC)) + geom_col()
  interval <- cut(pca$PC, 5, labels = F)
  pca[interval == 1,]
}

scree_plot <- function(pca){
  
  #Plot the scree plot of PCA object
  
  vars <- pca$sdev**2
  rel_vars <- vars / sum(vars)
  as.data.frame(rel_vars) %>% mutate(PC = paste0("PC", 1:length(rel_vars))) %>% ggplot(aes(x = reorder(PC, -rel_vars), y = rel_vars)) + geom_col() + labs(title = "Scree plot", subtitle = "Relative variance of each PC", x = "PC", y = "quotient of variance") + theme_bw()
}