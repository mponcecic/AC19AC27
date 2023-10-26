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


# Wilcoxon test
# Data frame with the p-value, adjusted p-value and the log2FC per each gene for 
# each comparison individually after performing the Wilcoxon test
wilcoxon_test <- function(count_cpm, correction, metadata, trt, contrast){
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
