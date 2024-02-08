################################################################################
#                             FUNCTIONS
################################################################################

# Summary
# ----------
#
# Contain all the functions created and run in the differentially expressed 
# genes analysis scripts (06_DEG_vX.R)

# List of functions
# 
# - create_contrast
# - callback
# - check_congruence
# - color_palette
# - wilcoxon_test
# - filter_genecounts
# 
# Plots function list
#   - ggsave_bg
#   - pca_plot
#   - heatmap_plot
#   - hist_verif 
#   - volcano_plot
#   - waterfall_plot
#   - waterfall_top




################################################################################
#                             CONTRAST
################################################################################


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



################################################################################
#                             DESIGN FORMULA
################################################################################


# This function generates the design formula for the different analysis 
# 
# Steps
# 1. Check there are covariates
# 2. Check the values of the covariate are not equal to avoid colinearity in 
#   the model. If their values are equal the variable is not included in the 
#   model.
# 3. Create the design formula

design_condition <- function(analysis, trt, var_exp, metadata){
  
  var_design <- NULL
  # Covariates  
  if(is.null(var_exp) == FALSE){
    # Check if the covariates present equal values or not 
    for (k in 1:length(var_exp)){if(length(unique(metadata[, var_exp[k]]))>1){var_design <- c(var_design, var_exp[k])}}
    # Design model for DESeq2
    if(analysis == "DESeq2"){
      design_cond = ifelse(is.null(var_design) == FALSE, paste("~", paste(var_design, "_zscore", " +", sep = "", collapse = " "), trt, sep = " "), paste("~", trt, sep = " "))
    # Design model for limma-voom and EdgeR
    }else{design_cond <- ifelse(!is.null(var_design), paste("~ 0 +", paste0(var_design, "+", collapse = " "), trt), paste("~ 0 +", trt))}
    
  # No Covariate
  } else {design_cond = ifelse(analysis == "DESeq2", paste("~", trt), paste("~ 0 +", trt))}
  
  return(design_cond)
}



################################################################################
#                               CALLBACK
################################################################################

# Reorder cluster rows
# Control on the left and Treatment on the right size 
# Sort by the average distance
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")}
callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



################################################################################
#                               CONGRUENCE
################################################################################


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



################################################################################
#                               COLOR PALETTE
################################################################################


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
    suppressWarnings(x <- brewer.pal(length(lvl_order), palette))
    names(x) <- lvl_order
    color_list[[trt]] <- x
    
    # To avoid problems with the color palette, we make sure the levels from the 
    # present color_list the correct levels from the condition  
  } else(
    if(all.equal.list(names(color_list[[trt]]), lvl_order) == FALSE){
    suppressWarnings(x <- brewer.pal(length(lvl_order), palette))
    names(x) <- lvl_order
    color_list[[trt]] <- x
  })
  
  # Remove the third element of the condition colors when the condition levels 
  # are only two and the list present three due to brewer.pal which give a minimun 
  # of three values
  if(length(lvl_order)<3 & length(lvl_order)!= length(color_list[[trt]])){color_list[[trt]] <- color_list[[trt]][-3]}
  
  return(color_list)
}



################################################################################
#                             WILCOXON TEST
################################################################################


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
  log2FC <- log2(rowMeans(count_cpm[,metadata$Sample[which(metadata[,trt] == contrast[[i]][2])]])/rowMeans(count_cpm[,metadata$Sample[which(metadata[,trt] == contrast[[i]][3])]]))
  
  # Final result
  wilcoxon <- data.frame(log2FC = log2FC, pvalue = pvalue, padj = padj)
  return(wilcoxon)
}



################################################################################
#                       FILTERING GENE COUNTS
################################################################################


# Filter gene count matrix
filter_genecounts <- function(counts_comp, metadata, trt, min_count = 10, min_prop = 0.7, n_large = 30, min_total = 15){
  
  # The filtering process differed depending on the group sample size. The group 
  # sample size is the number of samples corresponding to one level of the 
  # condition. When the group sample size is over the large sample size, the 
  # filtering step will be more restrictive. The default large sample size 
  # (n_large) is set in 30.
  # For groups with large samples size, 
  
  
  # This function is based on the function filterByExpr used in EdgeR and 
  # limma-voom. Link: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
  
  
  # minimun count (min_count)
  # minimun total count (min_total)
  
  # DEGList
  comp <- DGEList(counts = counts_comp, group = metadata[,trt])
  lib_size <- comp$samples$lib.size
  cpm_counts <- cpm(comp, lib.size = lib_size)
  
  # Group sample size over the large sample size
  # 0 when no group is over the large sample size 
  # >0 when a group is over the large sample size
  sample_size <- as.vector(tabulate(metadata[[trt]]))
  h <- sum(sample_size >= n_large)
  
  # Minimum sample size
  min_size <- sample_size[which.min(sample_size)]
  
  # Cutoff CPM
  cpm_cutoff <- min_count/mean(lib_size)*1e6
  
  # Minimum sample size in large samples
  if(h > 0){min_size <- n_large + (min_size - n_large)*min_prop}
  
  # Filtering data
  keep_counts <- rowSums(cpm_counts >= cpm_cutoff) >= (min_size - 1e-14)
  keep_total <- (rowSums(counts_comp) >= min_total -1e-14)

  
  # Output
  keep_total & keep_counts
}



# Filter gene count matrix
# filter_genecounts <- function(counts_comp, metadata, trt, min_count = 10, min_prop = 0.7, n_large = 30, min_total = 15){
#   
#   # The filtering process differed depending on the group sample size. The group 
#   # sample size is the number of samples corresponding to one level of the 
#   # condition. When the group sample size is over the large sample size, the 
#   # filtering step will be more restrictive. The default large sample size 
#   # (n_large) is set in 30.
#   # For groups with large samples size, 
#   
#   
#   # This function is based on the function filterByExpr used in EdgeR and 
#   # limma-voom. Link: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
#   
#   
#   # minimun count (min_count)
#   # minimun total count (min_total)
#   
#   # DEGList
#   comp <- DGEList(counts = counts_comp, group = metadata[,trt])
#   lib_size <- comp$samples$lib.size
#   cpm_counts <- cpm(comp, lib.size = lib_size)
#   
#   # Group sample size over the large sample size
#   # 0 when no group is over the large sample size 
#   # > 0 when a group is over the large sample size
#   sample_size <- as.vector(tabulate(metadata[[trt]]))
#   h <- sum(sample_size >= n_large)
#   
#   # Minimum sample size
#   min_size <- sample_size[which.min(sample_size)]
#   
#   # Cutoff CPM
#   cpm_cutoff <- min_count/mean(lib_size)*1e6
#   
#   # Modifications when sample group >  large samples size
#   if(h > 0){
#     # This step is done by Ivana 
#     keep_total <- (rowSums(counts_comp > 5) >= (min_size*min_prop))
#     
#     # Minimum sample size in large samples
#     min_size <- n_large + (min_size - n_large)*min_prop
#   } else {keep_total <- (rowSums(counts_comp) >= min_total -1e-14)}
#     
#   # Filtering data
#   keep_counts <- rowSums(cpm_counts >= cpm_cutoff) >= (min_size - 1e-14)
# 
#   
#   # Output
#   keep_counts & keep_total
# }



################################################################################
#                               PLOTS
################################################################################

################################################################
#                         MA plot
################################################################

MA_plot <- function(res, analysis, fdr_cutoff){
  
  # MA plot
  # 
  # The MA plot visualizes the differences between measurements taken in two 
  # samples groups, by transforming the data onto M (log ratio) and A (mean 
  # average) scales. The significant genes, based on the adjusted p-value, are 
  # plotted with a different color.
  # 
  # The aim of this function is to generate a MA plot for three of the method 
  # used in a single function and with the same aesthetic. 
  # 
  # The MA plot is slightly different between DESeq2 and EdgeR/limma-voom.
  # In DESeq2 data, we plot the log10 mean of normalized counts against the 
  # log2 fold change. However, the  data must be scaled to log10.
  # In EdgeR/limma-voom data, we plot the log10 counts per million against the 
  # log2 fold change.
  # 
  # In both graphs, a horizontal line with the intercept in 0 is added.
  # 
  # Point symbols can change from a dot (16) to a triangle when they are out of 
  # the y-axis limits, when log2FC is positive the triangle point up (2) and 
  # negative, the triangle point down (6)
  # 
  # This function is based on `MAplot` from DESeq2 which was dapted from 
  # `geneplotter`. 
  # Link to source code https://github.com/thelovelab/DESeq2/blob/devel/R/plots.R
  
  # Column with the adjusted p-value
  test.col <- "padj"
  
  # Y axis column
  lfc.col <- "log2FC"
  ylab <- "log 2 fold change"
  
  # Select X axis column and label for the different methods and DEGs color
  if (analysis == "DESeq2" | analysis == "DESeq2_NoFilter"){
    # X axis column
    counts.col <- "MeanExp"
    xlab <- "Mean of normalized counts"
    # Color DEGs
    colSig <- "blue"
  } else {   
    # EdgeR and limma-voom
    # X axis column
    counts.col <- "log2CPM"
    xlab <- "log 2 Counts per million"
    # Color DEGs
    colSig <- ifelse(analysis == "EdgeR", "red4", "darkgreen")
  } 
  
  # Line color
  colLine <- "grey40"
  # Color non-significant genes
  colNonSig <- "gray60" 
  
  
  # Plot data frame
  # 
  # From the input object, res, create a new data frame used to plot the data.
  # The data frame with 3 columns named mean, lfc and sig, and the rows corresponds 
  # to genes. In addition, all the genes with a log2FC value of 0 are not considered.
  # 
  #   - mean: Can be mean of normalized counts in DESeq2 or to log2 counts 
  #     per million in EdgeR and limma-voom. It is the plot x-axis.
  #   - lfc: Log2 fold change
  #   - sig: Boolean variable, when TRUE is a significant gene (p-adjusted > fdr) 
  #     and FALSE, non-significant gene
  df <- data.frame(mean = res[[counts.col]],
                   lfc = res[[lfc.col]],
                   sig = ifelse(res[[test.col]] <= fdr_cutoff, TRUE, FALSE))
  df <- subset(df, mean != 0)
  
  
  # LogFC variable 
  py <-df$lfc
  
  
  # Plot
  # Set y axis limits
  ylim <- c(-1,1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 1.1
  
  cex = 0.45
  
  if(analysis == "DESeq2"){
    y_axis <- pmax(ylim[1], pmin(ylim[2], py))
    
    # Data is transform into log10 scale to improve the interpretation in DESeq2
    plot(df$mean, y_axis, log = "x", pch = ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
         cex = cex, col = ifelse(df$sig, colSig, colNonSig), ylim = ylim,
         xlab = xlab, ylab = ylab); abline(h = 0, lwd = 4, col = colLine)
    
  } else { 
    plot(df$mean, df$lfc, pch = ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
         cex = cex, col = ifelse(df$sig, colSig, colNonSig), ylim = ylim,
         xlab = xlab, ylab = ylab); abline(h = 0, lwd = 4, col = colLine)
  } 
}



################################################################
#                         PCA
################################################################

pca_plot <- function(m, trt, metadata, color_l) {
  
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
  
  var_perc <- get_eigenvalue(m_pca)[,2]
  
  # Scree plot
  # Percentage of variances explained by each principal component
  pca_scree <- fviz_eig(m_pca, choice = "variance", addlabels = TRUE, ggtheme = theme_classic(), main = "Screeplot")
  
  # PC1 vs PC2
  pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1, 2),
                           geom.ind = c("point", "text"), 
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]],  fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC2 (", round(var_perc[2],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC3
  pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1, 3),
                           geom.ind = c("point", "text"),  
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC3 (", round(var_perc[3],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC4
  pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1, 4),
                           geom.ind = c("point", "text"),  
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC4 (", round(var_perc[4],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  
  # How to get the contribution (%) of each parameter to a PC?
  # 
  # First: Select the "rotation" and absolute value
  # --------------------------------------------------
  # The rotation are the principal components (the eigenvectors of the covariance 
  # matrix), in the original coordinate system. Typically a square matrix (unless 
  # you truncate it by introducing tolerance) with the same number of dimensions 
  # your original data had.
  pca_abs <- abs(m_pca$rotation)
  
  # Second: Scale data between 0 and 1
  # -----------------------------------
  pca_vals <- sweep(pca_abs, 2, colSums(pca_abs), "/")
  
  
  # Return the plots
  return(list(pca_scree, pca_1vs2, pca_1vs3, pca_1vs4, m_pca$rotation, pca_vals))
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
             clustering_method = "ward.D2", 
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
             clustering_method = "ward.D2", 
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
  A <- ggplot(data = res_df, aes(x = log2FC)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "log2FC", y = "Counts")+
    theme(text = element_text(size = 6))
  # Histogram adjusted p-value distribution
  B <- ggplot(data = res_df, aes(x = padj)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "Adjusted p-value", y = "Counts")+
    theme(text = element_text(size = 6))
  
  # DEGs
  # Histogram log 2 FC distribution 
  C <- ggplot(data = df, aes(x = log2FC)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "log2FC", y = "Counts")+
    theme(text = element_text(size = 6))
  # Histogram p-value adjusted distribution
  D <- ggplot(data = df, aes(x = padj)) +
    geom_histogram( fill = "#6696CC", color = "black")+
    labs(x = "Adjusted p-value", y = "Counts")+
    theme(text = element_text(size = 6))
  
  fig <- ggarrange(A, B, C, D, ncol = 2, nrow = 2)
  
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
  
  
  # Variable for the degs Symbol names used in volcano plots
  data$deglabel <- NA
  data$deglabel[which(data$DEG == "YES")] <- data$Symbol[which(data$DEG == "YES")]
  
  
  # Volcano plot with the DEGs in blue and the non-significant in grey
  A <- ggplot(data = data, aes(x = log2FC, y = -log10(padj), col = DEG))+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
    geom_point(size = 2, alpha = 0.6)+
    scale_color_manual(values = c("grey", "blue"), labels = c("Not significant", "Significative"))+
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)")+
    theme(text = element_text(size = 6), axis.title = element_text(size = rel(1.25)), legend.position = "none")
  
  # Volcano plot with the DEGs in separate in two color for the up and downregulated 
  # genes, the non-significant in grey
  B <- ggplot(data = data, aes(x = log2FC, y = -log10(padj), col = Direction))+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
    geom_point(size = 2, alpha = 0.6)+
    scale_color_manual(values = as.vector(color_list$Direction))+
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "")+
    theme(text = element_text(size = 6), axis.title = element_text(size = rel(1.25)), legend.position = "none")
  
  # Volcano plot with the DEGs in separate in two color for the up and downregulated 
  # genes, the non-significant in grey, with gene names for the degs
  C <- ggplot(data = data, aes(x = log2FC, y = -log10(padj), col = Direction, label = deglabel))+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
    geom_point(size = 2, alpha = 0.6)+
    geom_text_repel(hjust = "outward", max.overlaps = 10, size = 2, na.rm = TRUE, min.segment.length = 0.6)+
    scale_color_manual(values = as.vector(color_list$Direction))+
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "")+
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title = element_text(size = rel(1.25)), legend.position = "none")
  
 
  
  # Return the plots
  return(list(A, B, C))
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
  wf <- df[which(!(duplicated(df$Ensembl))),]
  # Sort the data based on the log2 FC
  wf$Ensembl <- factor(wf$Ensembl, levels = wf$Ensembl[order(wf$log2FC, decreasing = FALSE)])
  
  ggplot(wf, aes(x = Ensembl, y = log2FC, fill = Direction)) +
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
  wf <- df[which(!(duplicated(df$Ensembl))),]
  # Sort the data based on the log2 FC
  wf$Ensembl <- factor(wf$Ensembl, levels = wf$Ensembl[order(wf$log2FC, decreasing = FALSE)])
  
  # Select the top 30 genes 
  wf_sel <- rbind(wf[1:top,], wf[(dim(wf)[1]-top):dim(wf)[1],]) 
  
  ggplot(wf_sel, aes(x = log2FC, y = Ensembl, fill = Direction)) +
    geom_bar( stat = "identity") +
    xlab("Log2 Fold Change") +
    ylab("Gene name") +
    scale_fill_manual(values = as.vector(color_l$Direction)[-2])+
    theme(legend.position = "none")
}

