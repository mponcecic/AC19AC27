################################################################
#                       FUNCTIONS
################################################################


# Summary
# ----------
#
# Contain all the functions created and run in this workflow



################################################################
#              Sample information/ Metadata
################################################################


# Script: 00_Sample_info.R


# From pdf to data frame
pdf_to_tab <- function(x, last_sample){
  
  # Description
  # 
  # Transform the table from the pdf into a data frame to select the variables of 
  # interest
  # This can not be completely automatized because the tables differ from projects
  
  # Create a character vector using \n as the row separator
  x <- map(x, ~ str_split(.x, "\\n") %>% unlist())
  x <- reduce(x, c)
  
  # Select the beginning of the table
  # Must be always the same
  if(is.na(str_which(x, "Kibrary ID")[1])){
    tab_start <- str_which(x, "GAP ID")[1]
  } else {tab_start <- str_which(x, "Library ID")[1]}
  
  # Select the row of the end
  tab_end <- str_which(x, last_sample)[1]
  
  # Select the rows of interest
  tab <- x[(tab_start):(tab_end)]
  
  # Change column separator form spaces to "|" and remove the first four rows
  tab <- str_replace_all(tab, "\\s{2,}", "|")
  tab <- tab[-(1:4)]
  
  # Create a data frame of the table using "|" as a separator 
  # Remove the first column beacuse it's empty
  tab_df <- as.data.frame(do.call(rbind, strsplit(tab, "|", fixed = TRUE)))[,-1] 
  
  return(tab_df)
}



################################################################
#               Differential Expressed Genes
################################################################


# Script: 06_DEG_vX.R


# Contrast matrix function
# The contrast matrix is a list of contrasts to perform
# The arguments of each vector should be: 
#   - First. Variable name or Experimental condition to test
#   - Second. Level for the comparison
#   - Third. Level for the comparison and the one use as baseline
create_contrast <- function(trt, lvl_ord) {
  contrast <- list()
  for (i in 2:length(lvl_ord)) {for (j in 1:(i - 1)) {contrast[length(contrast) + 1] <- list(c(trt, lvl_ord[i], lvl_ord[j]))}}
  return(contrast)
}


# Reorder cluster rows
# Control on the left and Treatment on the right size 
# Sort by the average distance
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")}
callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Congruence function
# 
# Check if the direction of an event in a tissue is the same as in 
# other tissue
# Can be more than comparisons groups
check_congruence <- function(x){
  a <- unlist(strsplit(x, split = ", "))
  if (all(a ==a[1])){
    return("Yes")
  } else {return("No")}
}


# Color list function
color_palette <- function(color_list, trt, lvl_order, palette = "Dark2"){
  
  # Description
  #
  # This function aims to automatized the generation of the color palette for 
  # the different levels in the condition contrast
  
  # RColorBrewer palette options
  # Run display.brewer.all(colorblindFriendly = TRUE) to choose the palette you 
  # want, if not the default option is Dark2
  
  # Number lists in the list
  n <- length(names(color_list))
  
  # If the number list is than 4 it means that the condition color palette was 
  # not manually added and must be added now to the color_list to visualize data
  if(n < 4){
    x <- brewer.pal(length(lvl_order), palette)
    names(x) <- lvl_order
    color_list[[trt]] <- x
    
    # To avoid problems with the color palette, we make sure the levels from the 
    # present color_list the correct levels from the condition  
  } else(if(all.equal.list(names(color_list[[trt]]), lvl_order) == TRUE){
    x <- brewer.pal(length(lvl_order), palette)
    names(x) <- lvl_order
    color_list[[trt]] <- x
  })
  
  return(color_list)
}


# Wilcoxon test
wilcoxon_test <- function(count_cpm, correction, metadata, trt, contrast){
  
  # Data frame with the p-value, adjusted p-value and the log2FC per each gene for 
  # each comparison individually after performing the Wilcoxon test
  
  # P-value
  pvalue <- sapply(1:nrow(count_cpm), function(k){
    m_comp <- cbind.data.frame(gene = as.numeric(t(count_cpm[k,])), metadata[,trt])
    p <- wilcox.test(gene ~ metadata[,trt], m_comp, paired = FALSE)$p.value
    return(p)})
  
  # Adjusted p-value 
  padj <- p.adjust(pvalue, method = correction)
  
  # Log2 Fold Change
  logFC <- log2(rowMeans(count_cpm[,metadata$Sample[which(metadata[,trt] == contrast[[i]][2])]])/rowMeans(count_cpm[,metadata$Sample[which(metadata[,trt] == contrast[[i]][3])]]))
  
  # Final result
  wilcoxon <- data.frame(logFC = logFC, pvalue = pvalue, padj = padj)
  return(wilcoxon)
}



################################################################
#                         PLOTS
################################################################


# Alternative to function ggsave()
# 
# Here we can make major modifications of the plots and adjust them 
# 
# Modifications made:
# - Background set white as default option instead of transparency
# - Units set in cm
ggsave_bg <- function(filename, plot, path, width, height) {ggsave(filename = filename, plot = plot, path = path, width = width, height = height, bg = "white", units = "cm", scale = 1, dpi = 300)}


  
################################################################
#                         PCA
################################################################

perform_pca_analysis <- function(m, trt, metadata, color_l) {
  
  # PCA plots
  #
  # To perform the PCA analysis we need a matrix with the PSI values per each 
  # sample. The samples must be place in the rows and the events in the columns. 
  # This matrix is centered and scale by using prcomp().
  #
  # A data frame with the original data and the metadata must be used to plot.
  # This data is used to increase interpretation power.
  #
  # The principal components which explains the variability of the data should
  # at least sum 65% of the variance of the data.
  
  m_t <- t(m)
  print(dim(m_t))
  m_pca <- prcomp(m_t, scale. = TRUE)
  
  # Save output
  pca1 <- m_pca$x
  pca2 <- m_pca$rotation
  
  # Scree plot
  # Percentage of variances explained by each principal component
  pca_scree <- fviz_eig(m_pca, choice = "variance", addlabels = TRUE, ggtheme = theme_classic(), main = "Screeplot")
  
  # Variable graphs
  # Contributions of the events in each principal component
  # Positive correlated variables point to the same side of the plot
  # Negative correlated variables point to opposite sides of the graph.
  pca_var <- fviz_pca_var(m_pca, col.var = "contrib", gradient.cols = c("blue", "yellow", "red")) +
    theme_classic() + 
    labs(title = "Variance Contribution")
  
  fig <- ggarrange(pca_scree, pca_var, ncol = 2, nrow = 1, widths = 10, heights = 4)
  
  # PC1 vs PC2
  pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1, 2),
                           geom.ind = "text", repel = TRUE, labelsize = 4,
                           col.ind = metadata[[trt]],
                           addEllipses = TRUE, ellipse.level = 0.95,
                           title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    theme(legend.position = "none")
  
  # PC1 vs PC3
  pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1, 3),
                           geom.ind = "text", repel = TRUE, labelsize = 4,
                           col.ind = metadata[[trt]],
                           addEllipses = TRUE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    theme(legend.position = "none")
  
  # PC1 vs PC4
  pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1, 4),
                           geom.ind = "text", repel = TRUE, labelsize = 4,
                           col.ind = metadata[[trt]],
                           addEllipses = TRUE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    theme(legend.position = "none")
  
  # Return the plots
  return(list(fig, pca_1vs2, pca_1vs3, pca_1vs4))
}


################################################################
#                     Heatmap
################################################################


heatmap_plot <- function(m, metadata, trt, color_l, callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}){
  
  # Heatmap
  #
  # Clustering of all the differentially expressed genes with different layers
  # of information such as the treatment. First, a Z-score is performed based on 
  # gene expression among the samples. Second, it's clustered by columns and rows
  # but the rows (genes) tree is not shown.
  # 
  # The heatmap function needs a
  #   - Matrix data (m): Samples as columns and genes as rows
  #   - Annotation samples (sample_col): Treatment corresponding to each sample
  
  # Column information
  # Sample columns information
  # Should only contain the samples of the tissue followed by the treatment
  sample_col <-  metadata %>% dplyr::select(all_of(trt))
  rownames(sample_col) <- metadata$Sample
  sample_col[trt] <- factor(metadata[[trt]])
  
  
  # Choosing heatmap font size
  if (dim(m)[2] <= 10){font_col = 7} else if (dim(m)[2] > 10 & dim(m)[2] < 20){font_col = 5} else {font_col = 2}
  if (dim(m)[1]>=60 & dim(m)[1]<200){font_row = 2} else if (dim(m)[1]>=200 & dim(m)[1]<250){font_row = 1} else if (dim(m)[1]<=60 & dim(m)[1]>=20){font_row = 2} else {font_row = 5} 
  
  
  # Do not show the gene name if there are more than 250 genes
  if(dim(m)[1]>=250){
   plot <- pheatmap(m,
             scale = "row",
             color = color_l$Heatmap,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_col = sample_col,
             annotation_colors = color_l,
             show_rownames = FALSE,
             fontsize_col = font_col,
             border_color = NA,
             treeheight_row = 0,
             clustering_callback = callback) 
    
  } else {
    plot <- pheatmap(m,
             scale = "row",
             color = color_l$Heatmap,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_col = sample_col,
             annotation_colors = color_l,
             show_rownames = TRUE,
             fontsize_row = font_row,
             fontsize_col = font_col,
             border_color = NA,
             treeheight_row = 0,
             clustering_callback = callback)
  }
  
  # Return the plot
  return(plot)
}


################################################################
#                     Histogram 
################################################################


hist_verif <- function(res_df = res_df, df = df){
  
  # Grid with a representation of the adjusted p-values and the log2 FC estimated 
  # for the genes and them only for the selected genes which overcame the thresholds 
  
  # All the genes
  # Histogram log 2 FC distribution
  A <- ggplot(data = res_df, aes(x = logFC)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "log2FC", y = "Counts")
  # Histogram adjusted p-value distribution
  B <- ggplot(data = res_df, aes(x = padj)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "Adjusted p-value", y = "Counts")
  
  # DEGs
  # Histogram log 2 FC distribution 
  C <- ggplot(data = df, aes(x = logFC)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "log2FC", y = "Counts")
  # Histogram p-value adjusted distribution
  D <- ggplot(data = df, aes(x = padj)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "Adjusted p-value", y = "Counts")
  
  fig <- ggarrange(A, B, C, D, ncol = 2, nrow = 2, widths = 10, heights = 10)
  
  # Return the plot
  return(fig)
}


################################################################
#                     Volcano
################################################################


volcano_plot <- function(data, color_list, lfc_cutoff = log2(1.5), fdr_cutoff = 0.05){
  
  # Volcano plots
  # 
  # Scatter-plot with the log2 fold-change and the -log10 adjusted p-value on the 
  # X and Y axes. In addition, 2 lines are drawn based on the log2 FC and false 
  # discovery rate (fdr) thresholds to highlight the significant genes. The significant 
  # genes are highlight using different colors from the color list.
  # 
  # If the thresholds are not specified we use the laboratory default thresholds 
  # such as FC = 1,5 and fdr = 0.05.
  # 
  # Input data: data = cont_res <- res_df
  
  # Volcano plot with the DEGs in blue and the non-significant in grey
  A <- ggplot(data = data, aes(x = logFC, y = -log10(padj), col = DEG))+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
    geom_point(size = 2, alpha = 0.6)+
    scale_color_manual(values = c("grey", "blue"), labels = c("Not significant", "Significative"))+
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)")+
    theme(text = element_text(size = 8), plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)), legend.position = "none")
  
  # Volcano plot with the DEGs in separate in tw color for the up and downregulated 
  # genes, the non-significant in grey
  B <- ggplot(data = data, aes(x = logFC, y = -log10(padj), col = Direction))+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
    geom_point(size = 2, alpha = 0.6)+
    scale_color_manual(values = as.vector(color_list$Direction))+
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "")+
    theme(text = element_text(size = 8), plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)), legend.position = "none")
  
  # Return the plots
  return(list(A = A, B = B))
}



################################################################
#                     Waterfall
################################################################


waterfall_plot <- function(data, color_list) {
  
  # Waterfall plot 
  # 
  # Barplot which represent the change in the log2 FC in the differentially 
  # expressed genes. In this figure, you can differentiate the pattern between 
  # the up and downregulated genes. 
  # 
  # Input data: data = df
  
  # Remove genes with duplicated gene ID
  wf <- df[which(!(duplicated(df$GeneID))),]
  # Sort the data based on the log2 FC
  wf$GeneID <- factor(wf$GeneID, levels = wf$GeneID[order(wf$logFC, decreasing = FALSE)])
  
  ggplot(wf, aes(x = GeneID, y = logFC, fill = Direction)) +
    geom_bar(stat = "identity")+
    scale_fill_manual(values = as.vector(color_l$Direction)[-2])+
    xlab("Differentially expressed genes")+
    ylab("Log2 Fold Change")+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
  
}

waterfall_top <- function(data, color_list, top_genes = 30) {
  
  # Waterfall plot top genes
  # 
  # Barplot which represent the change in the log2 FC in the differentially 
  # expressed genes in the top genes selected. This function looks for the top 
  # genes selected in the up and downregulated genes.
  # 
  # Default option is 30 top genes based on the log2 FC value
  # 
  # Input data: data = df
  
  # Top genes selected
  top <- top_genes/2 
  
  # Remove genes with duplicated gene ID
  wf <- df[which(!(duplicated(df$GeneID))),]
  # Sort the data based on the log2 FC
  wf$GeneID <- factor(wf$GeneID, levels = wf$GeneID[order(wf$logFC, decreasing = FALSE)])
  
  # Select the top 30 genes 
  wf_sel <- rbind(wf[1:top,], wf[(dim(wf)[1]-top):dim(wf)[1],]) 
  
  ggplot(top_up_down, aes(x = logFC, y = GeneID, fill = Direction)) +
    geom_bar( stat = "identity") +
    xlab("Log2 Fold Change") +
    ylab("Gene name") +
    scale_fill_manual(values = as.vector(color_l$Direction)[-2])+
    theme(legend.position = "none")
}






ggsave_bg(filename = paste("00_Validation_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), width = 10, height = 10, plot = fig, path = dir_fig)

pdf(paste(dir_fig, "/Heatmap_Zscore_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""), height = 4, width = 4, bg = "white", compress = TRUE)
dev.off()

pdf(paste(dir_fig, "/Volcano_DEGs_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 4, width = 4, bg = "white", compress = TRUE)
dev.off()

pdf(paste(dir_fig, "/Volcano_DEGs_Direction_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 4, width = 4, bg = "white", compress = TRUE)
dev.off()

pdf(paste(dir_fig, "/Waterfall_DEGs_",analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 4, width = 4, bg = "white", compress = TRUE)
dev.off()

pdf(paste(dir_fig, "/Waterfall_DEGs_top_",analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 4, width = 4, bg = "white", compress = TRUE)
dev.off()

ggsave_bg(filename = paste("PCA_params_", analysis, "_", threshold, "_", name, "_", project, "_scree.pdf", sep = ""), plot = results[[1]], path = dir_fig, height = 5, width = 6)
ggsave_bg(filename = paste(deparse(substitute(pca_1vs2)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""), plot = results[[2]], path = dir_fig, height = 5, width = 6)
ggsave_bg(filename = paste(deparse(substitute(pca_1vs3)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""), plot = results[[3]], path = dir_fig, height = 5, width = 6)
ggsave_bg(filename = paste(deparse(substitute(pca_1vs4)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""), plot = results[[4]], path = dir_fig, height = 5, width = 6)