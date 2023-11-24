################################################################################
#                               RUN DESEeq2
################################################################################

# Summary
# ---------- 
# 
# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory.
# 
# After performing the analysis 


################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Date of the log file 5_DEG_qc_XXXX.log
logdate <- "202311"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/log/5_DEG_qc_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/04_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
dir_out <- paste(path, project, "/05_DEG_ANALYSIS", sep = "")

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
var_exp <- c("RIN")

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}


# Filtering lowly expressed genes
# Default: TRUE
filter_cutoff <- logfile$filter_cutoff


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
analysis <- "DESeq2"



################################################################################
#                               LOAD FILES
################################################################################


# Load annotation data 
gen_annot <- read.table(paste(dir_out,"/Annotation_", project,".txt", sep = ""), header = TRUE)

  
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
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  
  # Create a folder for the analysis. This folder contains several folders 
  # classified in Results and Figures. 
  
  # Load output directory
  dir_outfolder <- paste(dir_out, "/", name, sep='')
  setwd(dir_outfolder)
  
  # Files folder
  dir_output <- paste(dir_outfolder,"/Results", sep='')
  # Quality control figures folder
  dir_fig_qc <- paste(dir_outfolder, "/Figures", sep='')
  
  # Analysis figures folder
  dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
  
  
  ##############################################################################
  #                               Load Data
  ##############################################################################
  
  
  # Load gene count filtered
  gene_counts <- read.table(paste(dir_output, "/GeneCount_", name, "_", project, ".txt", sep = ""))
  
  # Load sample information per comparison
  metadata <- read.table(paste(dir_output, "/Metadata_", name, "_", project, ".txt", sep = ""))
  metadata[,trt] <- factor(metadata[,trt], levels = c(control, experimental))
  
  # Select the contrast levels
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% c(experimental, control))]
  
  # DESeq2 design formula
  # Used to estimate the variance stabilization (VST) method proposed in DESeq2
  # 
  # Steps
  # 1. Check there are covariates
  # 2. Check the values of the covariate are not equal to avoid colinearity in 
  #   the model. If their values are equal the variable is not included in the 
  #   model.
  # 3. Create the design formula
  design_cond <- design_condition(analysis, trt, var_exp, metadata)
  print(design_cond)
  
  # Annotation for all the results 
  ref <- paste(analysis, "_", name, "_", project, sep = "")
  
  
  
  ##############################################################################
  #                                 DESeq2
  ##############################################################################
  
  
  # Step 1: Create DESeq object
  # ----------------------------------------------------------------------------
  #
  # DESeqDataSetFromMatrix contains all the information necessary to run DESeq. 
  # - Count matrix
  # - Metadata with the following columns: Sample, Time and RIN. Data associated 
  #       to the samples valuable to perform the contrast
  # - Design formula with the experimental conditions (Time) and the covariates (RIN)
  
  # The options of the SummarizedExperiment can be applied in this object.
  dds <- DESeqDataSetFromMatrix(countData = gene_counts,                     # Input: Count matrix 
                                colData = metadata,                          # Input: Metadata
                                design = eval(parse(text = design_cond)),    # Design matrix
                                tidy = FALSE,                                # Default; Option: TRUE = First column of count data is the row names for the count matrix
                                ignoreRank = FALSE                           # Default: Reserved for DEXSeq developers
  )
  
  
  # Step 2: Run DESeq
  # ----------------------------------------------------------------------------
  # 
  # Filtering by default
  #   - Genes with zero counts are remove
  #   - Gene with an extreme count outlier
  #   - Genes with low mean normalized counts
  
  dds <- DESeq(dds,                                                  # DESeqDataSet
               test = "Wald",                                        # default; Wald Test
               fitType = "parametric",                               # default; Dispersion estimates
               sfType = "ratio",                                     # default; Size factors estimates
               betaPrior = FALSE,                                    # default; betaPrior = FALSE
               minReplicatesForReplace = 7,                          # default; Outliers replaced
               useT = FALSE,                                         # default; Std Normal distr for Wald statistics
               # minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,   # default;  lower bound on the estimate count while fitting the GLM
               quiet = FALSE,                                        # default; Print a message per step performed
               
               ## Not used
               # full = design(object),                                # default; Not used
               # reduced,                                              # defaultt; Not used
               # Parallelization
               # parallel = FALSE,                                     # default; Not used
               # BPPARAM = bpparam()                                   # default; Not used
  )
  
  
  # Fit of the dispersion estimates
  pdf(file = paste(dir_fig, "/00_Dispersion_DESeq2_", project, ".pdf", sep = ""), width = 6, height = 5)
  plotDispEsts(dds, genecol = "black", fitcol = "red", finalcol = "dodgerblue", legend = TRUE, log = "xy", cex = 0.45)
  dev.off()
  
  
  # Step 3: Contrast genes with Wald test
  # ----------------------------------------------------------------------------
  
  # Log2 fold change result table for an specific comparison
  resl <- results(object = dds,                      # DESeqDataSet  
                  contrast = contrast[[i]],          # Constrast 
                  test = "Wald",                     # default; Wald Test
                  minmu = 0.5,                       # default; lower bound on the estimate count while fitting the GLM
                  
                  ## Independent filtering  
                  alpha = 0.1,                       # default = 0.1
                  pAdjustMethod = correction,        # default = "BH"
                  lfcThreshold = 0,                  # default
                  independentFiltering = FALSE,      # Independent filtering 
                  format = "DataFrame"               # default, Result format
                  # theta,                             # Quantiles at which to assess the number of rejections
                  # filter,                            # default, mean of normalized counts
                  # filterFun,                         # optional; for performing independent filtering and p-value adjustment
                  # 
                  # altHypothesis = "greaterAbs",     #c("greaterAbs", "lessAbs", "greater", "less"),
                  # listValues = c(1, -1),
                  # cooksCutoff,                      # Threshold Cook's distance
                  # 
                  
                  # 
                  # ## Not used
                  # addMLE = FALSE,                   # default; Not used
                  # tidy = FALSE,                     # default; Not used
                  # saveCols = NULL,                  # default; Not used
                  # name,                             # default; Not used
                  # 
                  # ## Parallelization
                  # parallel = FALSE,                 # default; Not used
                  # BPPARAM = bpparam(),              # default; Not used
  )
  
  
  # Results as data frame 
  res <- as.data.frame(resl)
  res <- res[match(rownames(gene_counts), rownames(res)),]
  
  
  # Change columns names to plot data 
  colnames(res) <- c("MeanExp","logFC", "lfcSE", "stat", "pvalue", "padj")

  # MA plot 
  pdf(paste(dir_fig, "/00_MA_plot_", ref,".pdf", sep = ""), height = 4, width = 5)
  MA_plot(res, analysis, fdr_cutoff)
  dev.off()
  
  
  
  ##############################################################################
  #                            Data Processing
  ##############################################################################

  
  ## Results as a data frame
  res_df <- res
  res_df$Ensembl <- rownames(res_df)
  print(dim(res_df))
  
  ## Threshold label
  threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")
  
  ## Significative genes
  # Events with p-val NA are saved too
  # Based on the alpha and log2 FC thresholds  
  res_df$DEG[res_df$padj <= fdr_cutoff & abs(res_df$logFC) > lfc_cutoff] <- "YES"
  res_df$DEG[res_df$padj > fdr_cutoff & abs(res_df$logFC) <= lfc_cutoff] <- "NO"
  res_df$DEG[is.na(res_df$DEG)] <- "NO"
  
  ## Gene direction
  # Contains if DEGs are up (Upregulated) or downregulated (Downregulated)
  # Based on the alpha and log2 FC 
  res_df$Direction[res_df$padj <= fdr_cutoff & res_df$logFC > lfc_cutoff] <- "Upregulated"
  res_df$Direction[res_df$padj <= fdr_cutoff & res_df$logFC < -lfc_cutoff] <- "Downregulated"
  res_df$Direction[res_df$padj > fdr_cutoff & res_df$logFC <= lfc_cutoff] <- "Not significant"
  res_df$Direction[is.na(res_df$Direction)] <- "Not significant"
  
  
  ## Annotated gene names in Symbol
  res_df <- merge(gen_annot, res_df, by = "Ensembl") 
  print(head(res_df))
  print(dim(res_df))
  
  ## MERGE WITH GENE COUNTS
  # Row names to a variable
  genes <- gene_counts
  genes$Ensembl <-  rownames(genes)
  # Merge gene_counts and comparison results
  result <- merge(x = res_df, y = genes, by = "Ensembl")
  
  ## Differential expressed genes
  # Select differentially expressed genes
  df <- result[which(result$DEG == "YES"),]
  dim(df)
  
  
  ## Variance stabilizing transformation
  # 
  # Variance stabilization methods in log2 scale to interpret the data
  # 
  # Choose VST for samples size group smaller than 30. Why? 
  # 
  # if you have many samples (e.g. 100s), the rlog function might take too long, and so the vst function will be 
  # a faster choice. The rlog and VST have similar properties, but the rlog requires fitting a shrinkage term for 
  # each sample and each gene which takes time. See the DESeq2 paper for more discussion on the differences 
  # (Love, Huber, and Anders 2014)
  
  # Why blind = FALSE 
  # 
  
  # Output
  # - Columns are samples
  # - Rows are genes
  # 
  #                    N_3_E1    N_2_E4    N_2_E3   N_2_E2   H4_2_E4   H4_2_E3   H4_2_E2   H4_2_E1
  # ENSG00000000003 11.075692 10.993586 11.209109 11.03978 11.136152 11.147595 11.118555 11.125602
  # ENSG00000000419 11.959761 11.938768 12.033994 11.84389 11.846617 11.951190 11.951552 11.888236
  # ENSG00000000457 10.476816 10.557566 10.470185 10.47396 10.568444 10.488905 10.569326 10.394093
  
  # Estimate the biggest group sample size
  group_n <- max(as.vector(tabulate(metadata[[trt]])))
  
  # Data transformation
  if(group_n < 30){
    res_log2 <- as.data.frame(assay(vst(dds, blind = FALSE)))
    md <- "VST"
  }else{
    res_log2 <- as.data.frame(assay(rlog(dds, blind = FALSE)))
    md <- "RLOG"}
  print(md) 
  
  ## Transform matrix 
  # Select the differentially expressed genes that overcame the test
  # Used to plot the data 
  m <- res_log2[which(rownames(res_log2) %in% df$Ensembl), ]
  
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
  pem <- cor(m, method = "pearson", use = "na.or.complete")
  plot_h <- pheatmap(pem, color = colorRampPalette(brewer.pal(9, "Blues"))(255),
                     cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 6, fontsize_col = 6, border_color = NA, treeheight_row = 0, treeheight_col = 0)
  
  pdf(paste(dir_fig, "/Cor_pearson_", ref, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_h)
  dev.off()
  
  
  ## PCA PLOTS
  plot_pcas <- pca_plot(m, trt, metadata, color_l)
  
  ggsave(filename = paste("PCA_params_", ref, ".pdf", sep = ""), plot = plot_pcas[[1]], path = dir_fig, height = 4, width = 4, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs2)), ref, ".pdf", sep = ""), plot = plot_pcas[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs3)), ref, ".pdf", sep = ""), plot = plot_pcas[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs4)), ref, ".pdf", sep = ""), plot = plot_pcas[[4]], path = dir_fig, height = 5, width = 6, bg = "white")
  
  ## HEATMAP
  plot_heatmap <- heatmap_plot(m, metadata, trt, color_l)
  
  pdf(paste(dir_fig, "/Heatmap_zscore_", ref, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_heatmap)
  dev.off()
  
  
  ## VOLCANO
  
  volcano <- volcano_plot(res_df, color_list = color_l, lfc_cutoff, fdr_cutoff)
  
  ggsave(filename = paste("Volcano_", ref, ".pdf", sep = ""), plot = volcano[[1]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste("Volcano_color_", ref, ".pdf", sep = ""), plot = volcano[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  
  
  ## WATERFALL
  
  waterfall_pl <- waterfall_plot(df, color_l)
  waterfall_plot_top <- waterfall_top(df, color_l)
  
  ggsave(filename = paste("Waterfall_", ref, ".pdf", sep = ""), plot = waterfall_pl, path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste("Waterfall_top_genes_", ref, ".pdf", sep = ""), plot = waterfall_plot_top, path = dir_fig, height = 5, width = 6, bg = "white")
  
  
  ##############################################################################
  #                               Save data 
  ##############################################################################
  
  
  # Summary table
  sum_res <- c()
  sum_res$Contrast <- name
  sum_res$Method <- analysis
  sum_res$Design <- design_cond
  sum_res$Transformation <- md
  sum_res$Genes <- dim(res_df$DEG)
  sum_res$DEG <- sum(res_df$DEG == "YES")
  sum_res$Up <- sum(res_df$Direction == "Upregulated")
  sum_res$Down <- sum(res_df$Direction == "Downregulated")
  write.csv(sum_res, paste(dir_output, "/Summary_tab_", analysis, ";", ref, "_", threshold,".csv", sep = ""))
  
  # All results
  colnames(res_log2) <- paste(md, colnames(res_log2), sep = "_")
  data <- cbind(result, res_log2)
  data <- data %>% select(Name, Symbol, Ensembl, DEG, Direction, logFC, padj, MeanExp, lfcSE, stat, pvalue, everything())
  
  write.table(data, paste(dir_output, "/", ref, ";All_", md, "blindFALSE_", threshold,".txt", sep = ""), row.names = FALSE)
  write.xlsx(data, paste(dir_output, "/", ref, ";All_", md, "blindFALSE_", threshold,".xlsx", sep = ""), overwrite = TRUE)
  
  # Differential expressed genes
  colnames(m) <- paste(md, colnames(m), sep = "_")
  sel <- cbind(df, m)
  sel <- sel %>% select(Name, Symbol, Ensembl, DEG, Direction, logFC, padj, MeanExp, lfcSE, stat, pvalue, everything())
  write.table(data, paste(dir_output, "/", ref, ";DEGs_", md, "blindFALSE", threshold,".txt", sep = ""), row.names = FALSE)

    }

 

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

