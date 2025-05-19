################################################################################
#                   TRAJECTORIES SUBLUSTERING PLOTS 
################################################################################

# Date: Mon May 19 15:49:03 2025
# Author: Maria Ponce

# Summary
# ---------------------

# Subdivided the different clusters into subgroups using hierchical clustering. 
# Indeed, we plot the data with dendrograms and heatmaps.

# Input folder: W:/mponce/project/06_TRAJECTORIES/analysis_ID/Analysis/Algorithm/Results
# Output folder: W:/mponce/project/06_TRAJECTORIES/analysis_ID/Analysis/Algorithm/padj_0.05_log2FC_0.58_c_5/Subclustering_k



################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC19AC27"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "W:/mponce/"

# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID_degs <- "20240801092404"
analysis_ID_traj <- "20250319124539"

# Select analysis/method performed
analysis <- "DESeq2"

# Method used for trajectories
clustering <- "K-means"

# Optimal number of clusters
c <- 5

# Fuzziness coefficient
m <- NULL 

# Number of subclustering !!!!!
# Default used in trajectories analysis is 3
k <- 3
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# General output directory
dir_out <- paste(path, project, sep = "")

# Load libraries and functions
source(paste(dir_out, "/utils/libraries_trajectories.R", sep = ""))


# Load log file 
logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID_degs, ".log", sep = ""), header = TRUE)


# COndition
trt <- logfile$condition

# Condition levels
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))

# Significance level
fdr_cutoff <- logfile$fdr_cutoff

# Log2 fold change threshold
lfc_cutoff <- round(logfile$lfc_cutoff, 2)

# Organism
organism <- logfile$Organism

# Threshold label
threshold <- paste("padj_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep = "")

# Variance stabilization
md <- logfile$Variance

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")))
names(color_list) <- c(trt, "Heatmap")
names(color_list[[trt]]) <- lvl_ord

# File name
if(clustering == "Mfuzz"){
  load_file <- paste(clustering, "_", analysis, "_", project, ";Norm_", threshold, "_m_", m,"_c_", c, "_", analysis_ID_traj, sep = "")
} else {load_file <- paste(clustering, "_", analysis, "_", project, ";Norm_", threshold, "_c_", c, "_", analysis_ID_traj, sep = "")}



################################################################################
#                        SET WORKING DIRECTORY
################################################################################


# Trajectory directory
dir_output <- paste(dir_out, "/06_TRAJECTORIES/", analysis_ID_traj, "/", analysis, "/", clustering, sep = "")

# Input directory 
dir_input <- paste(dir_output, "/Results", sep = "")

# Figure directory
dir_fig <- paste(dir_output, "/", threshold, "_c_", c, sep = "")
dir.create(file.path(dir_fig, paste("Subclustering_", k, sep = "")), showWarnings = FALSE)
dir_fig <- paste(dir_fig, "/Subclustering_", k, sep = "")
setwd(dir_fig)

# Input metadata directory
dir_metadata <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID_degs, "/", "Results", sep = "")



################################################################################
#                           METADATA DATA
################################################################################


# Load metadata
metadata <- read.table(paste(dir_metadata, "/Metadata_", project, "_", analysis_ID_degs, ".txt",  sep = ""))
head(metadata)
dim(metadata)

# Process metadata
metadata[[trt]] <- as.factor(as.character(metadata[[trt]]))

# Column name for heatmap
col_met <- data.frame(Condition = metadata[[trt]])
rownames(col_met) <- metadata$Sample




################################################################################
#                              CLUSTERS
################################################################################

# Result file
exc <- createWorkbook()


for (i in 1:c) {
  
  sheet <- i
  data <- read.xlsx(paste(dir_input, "/", load_file, ".xlsx", sep = ""), sheet = sheet)
  
  # 1. Select normalized counts
  # ----------------------------------------------------
  m <- data[, paste("Norm_", metadata$Sample, sep = "")]
  colnames(m) <- gsub("Norm_", "", colnames(m))
  rownames(m) <- data$Name

  # 2. Transform data (log2)
  # ----------------------------------------------------
  m <- log2(m+1)
  
  # 3. Plot heatmap
  # ----------------------------------------------------
  pheat <- pheatmap(mat = m,
                    scale = "row",
                    color = color_list$Heatmap,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    annotation_col = col_met,
                    annotation_colors = color_list,
                    annotation_names_row = FALSE,
                    annotation_names_col = FALSE,
                    annotation_legend = FALSE,
                    fontsize_row = 1,
                    show_rownames = FALSE,
                    border_color = NA,
                    clustering_method = "ward.D2", 
                    clustering_callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}
                    
  )
  
  
  # 4. Plot dendrogram
  # ----------------------------------------------------
  
  # Hierarchical clustering
  hc <- pheat$tree_row
  
  # Subclustering groups
  clust <- dendextend:::cutree(hc, k = k, order_clusters_as_data = FALSE)
  branch_colors <- colorspace::rainbow_hcl(n = k)
  
  # Plot subclustering
  pdf(paste(dir_fig, "/Dendrogram_n", nrow(m), "_cl_",sheet, "_", project, "_", analysis_ID_traj, ".pdf", sep = ""), height = 5, width = 7, bg = "white")
  plot(as.dendrogram(hc) %>% dendextend::color_branches(k = k), lwd = 2, cex = 0.01, leaflab = "none", main = paste("Cluster", sheet, sep = " "))
  rect.hclust(hc,k = k, border = branch_colors)
  legend_data <- data.frame(Cluster = 1:k, Color = branch_colors)
  legend("topleft", legend = legend_data$Cluster, fill = legend_data$Color, title = "Clusters", xpd = TRUE)
  dev.off()

  
  
  
  # Heatmap with dendrogam colored
  # 
  # Link: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  # 
  # BiocManager::install("ComplexHeatmap")
  row_dend = as.dendrogram(hc)
  row_dend = dendextend::color_branches(row_dend, k = k) # `color_branches()` returns a dendrogram object
  
  # Scale data by row
  m_scaled <- t(scale(t(m)))
  
  # # Column cluster order
  # col_dend <- pheat$tree_col  # Column clustering
  # 
  
  # Sample colors
  col_met <- data.frame(Condition = factor(metadata[[trt]]))
  rownames(col_met) <- metadata$Sample
  
  # 
  annotation_col <- ComplexHeatmap::HeatmapAnnotation(
    Condition = metadata[[trt]], 
    col = list(Condition = color_list[[trt]]), 
    show_legend = FALSE, 
    show_annotation_name = FALSE, 
    simple_anno_size = unit(2, "mm"), 
    annotation_name_gp = gpar(fontsize = 2))
  # 
  # library(ComplexHeatmap)
  # library(circlize)
  # 
  
  fixed_order <- colnames(m_scaled)
  
  pdf(paste(dir_fig, "/Heatmap_n", nrow(m), "_cl_", sheet, "_", project, "_", analysis_ID_traj, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(ComplexHeatmap::Heatmap(m_scaled, 
          name = "mat",
          cluster_rows = row_dend, 
          column_order = fixed_order,
          show_row_names = FALSE, 
          top_annotation = annotation_col,
          col = color_list$Heatmap,
          column_names_gp = gpar(fontsize = 6),
          row_dend_gp = gpar(lwd = 0.5),  # Change row dendrogram line width
          column_dend_gp = gpar(lwd = 0.5), # Change column dendrogram line width
          heatmap_legend_param = list(title_gp = gpar(fontsize = 5),  # Title size
                                      labels_gp = gpar(fontsize = 4), # Label size
                                      legend_width = unit(5, "mm"))
  ))
  dev.off()

    
  # 5. Process results data
  # ----------------------------------------------------
  
  # Sort genes to add the column
  it <- match(rownames(m), names(clust))
  clust_ord <- as.vector(clust[it])
  
  # Add column to the output
  m$Name <- rownames(m)
  m$Subgroup <- clust_ord
  Cluster <- m[, c(ncol(m)-1,ncol(m))]
  head(Cluster)

  # Prepare results data
  res <- merge(Cluster, data, by = "Name")
  res <- res %>% select(Name, Symbol, Ensembl, Biotype, Cluster, Subgroup, Dist_Centroid, everything())
  
  # Order rows based on heatmap order
  res <- res[pheat$tree_row$order,]
  
  addWorksheet(exc, sheet)
  writeData(exc, sheet, res)
}

# Save data
# ----------------------------------------------------
saveWorkbook(exc, file = paste(dir_fig, "/Subclustering_", k, ";", load_file, ".xlsx", sep = ""), overwrite = TRUE)

  

