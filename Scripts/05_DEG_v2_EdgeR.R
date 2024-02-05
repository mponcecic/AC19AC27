################################################################################
#                               RUN EdgeR
################################################################################


# Summary 
#------------

# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory. 
# 
# This analysis is divided in different scripts. First, perform the quality 
# control of the data (05_DEG_v1_qc). Second, perform the methods to estimate the 
# DEGs which are DESeq2 (05_DEG_v2_DESeq2), EdgeR (05_DEG_v2_EdgeR), limma-voom 
# (05_DEG_v2_limma-voom) and Wilcoxon rank-sum test (05_DEG_v2_wilcoxon). Third,
# compare the results among the different methods used using 05_DEG_v3_Comparison 
# script.
# 
# In this script, EdgeR method is perform based on the gene count matrix resulting 
# from the previous script.


# Recommendations and warnings
# ---------------------------------
# 
# If you want all the comparison to be included in the same input for applying 
# DESeq2, EdgeR, limma-voom or Wilcoxon test, this is not possible in this script. 
# We generate a unique matrix per each contrast which only includes the samples 
# of that contrast
# 
# In the design formula, interactions are not considered.
# 
# You will find a chunk of code where you can adjust all the parameters, 
# filtering, thresholds and colors needed to perform in this analysis. This 
# variables are deeply detailed in the code. Following the code, it's time to 
# load the data sets, be careful everything is matched.
# 
# Code chunk contain between the following separators can be modify and must be 
# modify to perform the correct analysis. Two chunk will be found one referring
# to the parameter selection and the order for loading the data.


# Authors: María Ponce




################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/user/"
path <- "W:/user/"

# Date of the log file 5_DEG_qc_XXXX.log
logdate <- "20231204"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/log/5_DEG_qc_", logdate, ".log", sep = ""), header = TRUE)

# Output directory
dir_out <- paste(path, project, "/05_DEG_ANALYSIS", sep = "")

# Input directory. Raw gene counts  
dir_infiles <- paste(dir_out,  "/Results/", sep = "")

# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- logfile$condition

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))

# Variance sources to include in the model
# Can be set NULL
# Mouse: RIN
# Human: RIN, dv200, Age, ... 
#
# Options
# var_exp <- c("Age", "dv200")
# var_exp <- NULL
var_exp <-  unlist(strsplit(logfile$Varexp, split = ","))

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}


# Filtering parameters 

# Minimum number of counts per gene
min_total <- logfile$min_total
# Minimum number of counts per sample
min_count <- logfile$min_count
# Large sample size 
n_large <- logfile$n_large
# Proportion
min_prop <- logfile$min_prop


### Threshold criteria 

## Significance level
# Default. pvalue < 0.05 
fdr_cutoff <- logfile$fdr_cutoff

## Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$lfc_cutoff

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- logfile$correction

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")),
                   Direction = unlist(str_split(logfile$colordir, pattern = ",")),
                   Shared = unlist(str_split(logfile$colorsh, pattern = ",")))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")
names(color_list[[trt]]) <- lvl_ord
names(color_list[["Direction"]]) <- c("Downregulated", "Not significant", "Upregulated")

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Specie 
# specie = "Human"
specie <- logfile$Organism 

# Method used to study DEGs
analysis <- "EdgeR"


################################################################################
#                               LOAD DATA
################################################################################


## Load metadata file 
sample_info <- read.table(paste(dir_infiles, "/Metadata_", project, ".txt", sep = ""))


## Load gene counts
raw_counts <- read.table(file = paste(dir_infiles, "GeneCount_filter_mincount_", min_count, "_mintotal_", min_total, "_", project, ".txt", sep = ""))


## Load transformed data 
m_vst <- read.table(file = paste(dir_infiles, "GeneCount_filter_CPM_", project, ".txt", sep = ""))


## Load annotation data
# Genome annotation path
anot_path <- sub(pattern = "/.*", replacement = "", path)
if (anot_path == ""){anot_path <- "/vols/GPArkaitz_bigdata"}
# Reference genome
ref_genome <- read.table(paste(anot_path, "/DATA_shared/Genomes_Rocky/DEG_Annotation/Annotated_Genes_20231121.txt", sep = ""), header = TRUE)
genome <- ref_genome[which(ref_genome$Specie == specie), -4]
genome$Name <- paste(genome$Symbol, genome$Ensembl, sep = "_")
print(dim(genome))
# Annotate all the genes in the matrix with the Symbol annotation
# CAUTION: Symbol id present more than one Ensembl identifier. 
# I propose to perform the analysis using the Ensembl identifier, if not the 
# mean value of the genes with the same Symbol id should be performed.
annot <- genome[which(genome$Ensembl %in% rownames(raw_counts)),]
gene_names <- annot[match(rownames(raw_counts), annot$Ensembl), ]
print(dim(gene_names))



################################################################################
#                             DATA PROCESSING
################################################################################

# Change name
raw_genes <- raw_counts

# Generate column names
col_raw <- c(colnames(raw_genes), paste("CPM", colnames(m_vst), sep = "_"))

# Bind raw, vst/rlog and normalized counts
raw_genes <- cbind(raw_genes, m_vst)  
colnames(raw_genes) <- col_raw

# Ensembl id column to merge with results
raw_genes$Ensembl <- rownames(raw_genes)



################################################################################
#                             CREATE DATAFRAMES
################################################################################

# Summary of the different comparisons results for this analysis found in the 
# 05_DEG_ANALYSIS/Result folder
sum_res <- data.frame()

# Create Workbook for the researchers
# Each sheet corresponds to one of the comparison results
# Columns: Name, Symbol, Ensembl, DEG, Direction, log2FC, padj, variables, 
# normalized counts, raw count
exc <- createWorkbook()

# Final results in *.txt
# Data frame with all the comparison results, same as exc
final_data <- data.frame()


################################################################################
#                               COMPARISONS
################################################################################


for (i in 1:length(contrast)){
  
  # Contrast 
  name <- paste(contrast[[i]][2], "vs", contrast[[i]][3], sep = "")
  print(name)
  
  # Contrast levels 
  control <- contrast[[i]][3]
  experimental <- contrast[[i]][2]
  
  # Comparison levels
  comp_lvl <- c(control, experimental)
  
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  
  # Create a folder for the analysis. This folder contains several folders 
  # classified in Results and Figures. 
  
  # Load output directory
  dir.create(file.path(dir_out , name), showWarnings = FALSE)
  dir_outfolder <- paste(dir_out, "/", name, sep='')
  setwd(dir_outfolder)
  
  # Save files with comparison results separate
  dir.create(file.path(dir_outfolder , "Results"), showWarnings = FALSE)
  dir_files <- paste(dir_outfolder, "/Results", sep='')
  
  # Figures folder
  dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
  
  
  
  ##############################################################################
  #                             Data processing
  ##############################################################################
  
  
  # Metadata
  metadata <- sample_info[which(sample_info[[trt]] %in% comp_lvl),]
  metadata[,trt] <- factor(metadata[,trt], levels = comp_lvl)
  
  # Gene count matrix per comparison
  gene_counts <- raw_counts[, metadata$Sample]
  gene_counts <- gene_counts[which(rowSums(gene_counts) != 0), ] 
  
  # Transformed data per comparison
  res_log2 <- m_vst[which(rownames(m_vst) %in% rownames(gene_counts)), metadata$Sample]
  
  # Select the contrast levels
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% c(experimental, control))]
  
  # Design formula 
  design_cond <- design_condition(analysis, trt, var_exp, metadata)
  print(design_cond)
  
  # Build a model matrix 
  m_model = model.matrix(as.formula(design_cond), metadata)
  
  # Annotation for all the results 
  ref <- paste(analysis, "_", name, "_", project, sep = "")
  print(ref)
  
  
  ##############################################################################
  #                                 EdgeR
  ##############################################################################
  
  
  # Step 1: Create DGEList object
  # ----------------------------------------------------------------------------
  #
  # DGEList contains all the information necessary to run EdgeR and limma-voom 
  # - Count matrix
  # - Metadata with the following columns: Sample, Time and RIN. Data associated 
  #       to the samples valuable to perform the contrast
  
  # The options of the SummarizedExperiment can be applied in this object.
  deg <- DGEList(counts = gene_counts, group = metadata[,trt])
  if(!is.null(var_exp)){deg$samples[var_exp] <- metadata[,var_exp]}
  
  
  # Step 2: Normalization
  # ----------------------------------------------------------------------------
  # 
  # This step is common between EdgeR and limma-voom 
  # 
  # Minimized the log2FC between samples for most genes using weighted Trimmed Mean of M-values 
  # between each pair of samples (TMM method) to normalized the gene counts based 
  # on the library size.
  # Corrects for sequencing detph, RNA composition and gene length
  # Can be used to perform DE analysis and within sampe analysis
  # 
  # All the parameters used are the default options
  deg <- normLibSizes(deg, method = "TMM", refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10, p = 0.75)
  
  
  # Step 3: Run EdgeR 
  # ----------------------------------------------------------------------------
  #
  # In this NB model, we use a GLM (Generalized linear model) which is used when 
  # multiple factor variables are used in the analysis. In this case, we have the 
  # experimental variables and the covariate correction.
  # 
  # First, estimates a common NB regression based on the the matrix model
  # Second, estimates the abundance-dispersion trend by Cox-Reid
  # Third, computes empirical Bayes estimate of the NB dispersion parameter for 
  # each gene (= tag), with expression levels specified by a log-linear model.
  #
  # This analysis can be performed using: 
  #       estimateDisp(deg, design = m_model, tagwise = TRUE, trend.method = "locfit", robust = TRUE)
  # The results are very similar to the classic method used here but are not equal
  
  degde <- estimateGLMCommonDisp(deg, design = m_model, method = "CoxReid", verbose = FALSE)
  degde <- estimateGLMTrendedDisp(degde, design = m_model, method = "auto")
  degde <- estimateGLMTagwiseDisp(degde, design = m_model)
  
  # Biological coefficient variation per gene 
  pdf(file = paste(dir_fig, "/00_BCV_", ref, ".pdf", sep = ""), height = 5 , width = 5)
  plotBCV(degde)
  dev.off()
  
  
  # Step 4: Perform the contrast 
  # ----------------------------------------------------------------------------
  # Perform contrast
  fit <- exactTest(degde, pair = c(control, experimental), dispersion = "auto", rejection.region = "doubletail")
  
  # Estimate adjusted p-value
  res <- as.data.frame(topTags(fit, adjust.method = correction, n = dim(gene_counts)[1]))
  
  # Change columns names to plot data 
  res <- res[match(rownames(gene_counts), rownames(res)),]
  colnames(res) <- c("log2FC", "log2CPM", "pvalue", "padj")
  
  
  # MA plot 
  pdf(paste(dir_fig, "/00_MA_plot_", ref,".pdf", sep = ""), height = 4, width = 5)
  MA_plot(res, analysis, fdr_cutoff)
  dev.off()
  
 
  
  ##############################################################################
  #                            Data Processing
  ##############################################################################
  
  
  ## Results as a data frame
  res_df <- as.data.frame(res) 
  res_df$Ensembl <- rownames(res_df)
  print(dim(res_df))
  
  ## Threshold label
  threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")
  
  ## Significative genes
  # Events with p-val NA are saved too
  # Based on the alpha and log2 FC thresholds  
  res_df$DEG[res_df$padj <= fdr_cutoff & abs(res_df$log2FC) > lfc_cutoff] <- "YES"
  res_df$DEG[res_df$padj > fdr_cutoff & abs(res_df$log2FC) <= lfc_cutoff] <- "NO"
  res_df$DEG[is.na(res_df$DEG)] <- "NO"
  
  ## Gene direction
  # Contains if DEGs are up (Upregulated) or downregulated (Downregulated)
  # Based on the alpha and log2 FC 
  res_df$Direction[res_df$padj <= fdr_cutoff & res_df$log2FC > lfc_cutoff] <- "Upregulated"
  res_df$Direction[res_df$padj <= fdr_cutoff & res_df$log2FC < -lfc_cutoff] <- "Downregulated"
  res_df$Direction[res_df$padj > fdr_cutoff & res_df$log2FC <= lfc_cutoff] <- "Not significant"
  res_df$Direction[is.na(res_df$Direction)] <- "Not significant"
  
  
  ## Annotated gene names in Symbol
  res_df <- merge(res_df, gene_names, by = "Ensembl") 
  print(head(res_df))
  print(dim(res_df))
  
  
  ## MERGE WITH GENE COUNTS
  # Row names to a variable
  genes <- gene_counts[, metadata$Sample]
  genes$Ensembl <-  rownames(genes)
  # Merge gene_counts and comparison results
  result <- merge(x = res_df, y = genes, by = "Ensembl")
  
  
  ## Differential expressed genes
  # Select differentially expressed genes
  df <- result[which(result$DEG == "YES"),]
  
  
  ## Transform matrix 
  # Select the differentially expressed genes that overcame the test
  # Used to plot the data 
  m <- res_log2[which(rownames(res_log2) %in% df$Ensembl),]
  
  if(identical(rownames(m), df$Ensembl) == FALSE){m <- m[match(rownames(m), df$Ensembl),]}
  
  
  ##############################################################################
  #                                 Plot
  ##############################################################################
  
  
  ## HISTOGRAMS
  # Representation of the adjusted p-value and log2 fold-change for all and 
  # significant genes
  
  plot_hist <- hist_verif(res_df, df)
  ggsave(filename = paste("01_Histogram_verif_", ref, ".pdf", sep = ""), plot = plot_hist, path = dir_fig, height = 4, width = 4, bg = "white")
  
  
  ##  CORRELATION 
  #
  # Execute a Pearson correlation which accept possible NAs.
  # The acceptance of the NA should be valuable in the future in case, we accept
  # comparisons when a sample is missing a PSI value.
  #
  # The results of the correlation matrix must be aligned with the results
  # in the heatmap and PCA.
  m0 <- m[which(rowSums(m) != 0), ] 
  pem <- cor(m0, method = "pearson", use = "na.or.complete")
  plot_h <- pheatmap(pem, color = colorRampPalette(brewer.pal(9, "Blues"))(255),
                     cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 6, fontsize_col = 6, border_color = NA, treeheight_row = 0, treeheight_col = 0)
  
  pdf(paste(dir_fig, "/Cor_pearson_", ref, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_h)
  dev.off()
  
  
  ## PCA PLOTS
  plot_pcas <- pca_plot(m0, trt, metadata, color_l)
  
  ggsave(filename = paste("PCA_params_", ref, ".pdf", sep = ""), plot = plot_pcas[[1]], path = dir_fig, height = 4, width = 4, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_", ref, ".pdf", sep = ""), plot = plot_pcas[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_", ref, ".pdf", sep = ""), plot = plot_pcas[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_", ref, ".pdf", sep = ""), plot = plot_pcas[[4]], path = dir_fig, height = 5, width = 6, bg = "white")
  
  ## HEATMAP
  plot_heatmap <- heatmap_plot(m0, metadata, trt, color_l)
  
  pdf(paste(dir_fig, "/Heatmap_zscore_", ref, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_heatmap)
  dev.off()
  
  
  ## VOLCANO
  
  volcano <- volcano_plot(res_df, color_list = color_l, lfc_cutoff, fdr_cutoff)
  
  ggsave(filename = paste("Volcano_", ref, ".pdf", sep = ""), plot = volcano[[1]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste("Volcano_color_", ref, ".pdf", sep = ""), plot = volcano[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste("Volcano_lablels_", ref, ".pdf", sep = ""), plot = volcano[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
  
  
  ## WATERFALL
  
  waterfall_p <- waterfall_plot(df, color_l)
  waterfall_plot_top <- waterfall_top(df, color_l)
  
  ggsave(filename = paste("Waterfall_", ref, ".pdf", sep = ""), plot = waterfall_p, path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste("Waterfall_top_genes_", ref, ".pdf", sep = ""), plot = waterfall_plot_top, path = dir_fig, height = 5, width = 6, bg = "white")
  
  
  
  ##############################################################################
  #                               Save data 
  ##############################################################################
  
  
  # Summary table
  # Summary table
  sum_line <- c(name, analysis, design_cond, "CPM", dim(gene_counts)[1], dim(res_df$DEG), sum(res_df$DEG == "YES"), sum(res_df$Direction == "Upregulated"), sum(res_df$Direction == "Downregulated"))
  sum_res <- rbind(sum_res, sum_line)
  
  
  # All results
  colnames(res_log2) <- paste("CPM", colnames(res_log2), sep = "_")
  data <- cbind(result, res_log2)
  data <- data %>% select(Name, Symbol, Ensembl, DEG, Direction, log2FC, padj, log2CPM, pvalue, everything())
  write.table(data, paste(dir_files, "/", ref, ";All_CPM_", threshold, ".txt", sep = ""), row.names = FALSE)
  
  # Save data in the workbook
  addWorksheet(exc, name)
  writeData(exc, data, sheet = name)
  
  # All comparisons results
  result2 <- merge(x = res_df, y = raw_genes, by = "Ensembl")
  result2$Comparison <- name
  result2 <- result2 %>% select(Comparison, Name, Symbol, Ensembl, Biotype, DEG, Direction, log2FC, padj, log2CPM, pvalue, everything())
  final_data <- rbind(final_data, result2)
  
}


# Save Workbook
saveWorkbook(exc, file =  paste(dir_infiles, "/", analysis, "_", project, ";All_CPM_", threshold, ".xlsx", sep = ""), overwrite = TRUE)

# Save all comparisons 
colnames(final_data) <- c("Comparison", "Name", "Symbol", "Ensembl", "Biotype", "DEG", "Direction", "log2FC", "padj", "log2CPM", "pvalue", col_raw)
write.table(final_data, file = paste(dir_infiles, analysis, "_", project, ";All_CPM", threshold, ".txt", sep = ""), sep = " ", row.names = FALSE, col.names = TRUE)

# Save selected genes in all comparison
final_data <- final_data[which(final_data$DEG == "YES"), ]
write.table(final_data, file = paste(dir_infiles, analysis, "_", project, ";Selected_CPM", threshold, ".txt", sep = ""), sep = " ", row.names = FALSE, col.names = TRUE)

# Save Summary table 
colnames(sum_res) <- c("Comparison", "Analysis", "Design", "Transformation", "Genes", "DEGs", "Upregulated", "Downregulated")
write.csv(sum_res, paste(dir_infiles, "/Summary_tab_", analysis, "_", project, "_", threshold, ".csv", sep = ""), row.names = FALSE)



################################################################################
#                                 LOG FILE 
################################################################################

# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- dir_out
log_data$Analysis <- analysis
log_data$project_name <- project
log_data$Organism <- specie
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_ord, collapse =",")
log_data$Outliers <- paste(outliers, collapse = ",") 
log_data$Varexp <- paste(var_exp, collapse = ",") 
log_data$min_count <- logfile$min_count
log_data$min_prop <- logfile$min_prop
log_data$n_large <- logfile$n_large
log_data$min_total <- log_data$min_total
log_data$fdr_cutoff <- fdr_cutoff
log_data$lfc_cutoff <- lfc_cutoff
log_data$correction <- correction
log_data$contrast <- paste(unlist(contrast), collapse = ",")
log_data$colortrt <- paste(color_list[[1]], collapse = ",")
log_data$colorheat <- paste(color_list[[2]], collapse = ",")
log_data$colordir<-  paste(color_list[[3]], collapse = ",")
log_data$colorsh <- paste(color_list[[4]], collapse = ",")

write.table(as.data.frame(log_data), paste(path, project, "/log/5_DEG_v2_", analysis, "_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")

