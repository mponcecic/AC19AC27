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
project <- "PRUEBA"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Date of the log file
logdate <- "20231106"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/1_DEG_qc_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/03_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
# dir_out <- paste(path, project, sep = "")   # Default option
dir_out <- paste(path, project, sep = "")

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
#                               COMPARISONS
################################################################################


for (i in 1:length(contrast)){
  
  # Contrast 
  name <- paste(contrast[[i]][2], "vs", contrast[[i]][3], sep = "")
  print(name)
  
  # Comparison levels
  comp_lvl <- c(contrast[[i]][2],  contrast[[i]][3])
  
  # Select the contrast levels
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% comp_lvl)]
  
  # Design formula 
  design_cond = ifelse(is.null(var_exp) == TRUE, paste("~", trt, sep = " "), paste("~", paste(var_exp, "_zscore", " +", sep = "", collapse = " "), trt, sep = " "))
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  
  # Create a folder for the analysis. This folder contains several folders 
  # classified in Results and Figures. 
  
  # Load output directory
  dir_outfolder <- paste(dir_out,"/", name, sep='')
  setwd(dir_outfolder)
  
  # Files folder
  dir_output <- paste(dir_outfolder,"/Results", sep='')
  # Quality control figures folder
  dir_fig_qc <- paste(dir_outfolder, "/Figures", sep='')
  
  # Analysis figures folder
  dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
  
  
  ##############################################################################
  #                               LOAD DATA 
  ##############################################################################
  
  
  # Load gene count filtered
  gene_counts <- read.table(paste(dir_outfolder, "/", "GeneCount_", name, "_filtered_", project, ".txt", sep = ""))
  
  # Load sample information per comparison
  metadata <- read.table(paste(dir_outfolder, "/", "Metadata_", name, "_", project, ".txt", sep = ""))
  metadata <- metadata[which(metadata[[trt]] %in% comp_lvl),]
  metadata[,trt] <- factor(metadata[,trt])
  
  
  
  ##############################################################################
  #                           PERFORM ANALYSIS
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
  
  
  # Step 3: Contrast genes with Wald test
  # ----------------------------------------------------------------------------
  # 
  #
  
  # Log2 fold change result table for an specific comparison
  resl <- results(object = dds,                      # DESeqDataSet  
                  contrast = contrast[[i]],          # Constrast 
                  test = "Wald",                     # default; Wald Test
                  minmu = 0.5,                       # default; lower bound on the estimate count while fitting the GLM
                  
                  ## Independent filtering  
                  alpha = 0.1,                       # default = 0.1
                  pAdjustMethod = correction,        # default = "BH"
                  lfcThreshold = 0,                  # default
                  independentFiltering = TRUE,       # Independent filtering 
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
  
  
  # MA plot 
  pdf(paste(dir_fig, "/00_MA_plot_", analysis, "_", project, "_", name ,".pdf", sep = ""), height = 4, width = 5)
  DESeq2::plotMA(resl)
  dev.off()
  
  res <- as.data.frame(resl)
  res <- res[match(rownames(gene_counts), rownames(res)),]
  colnames(res) <- paste(analysis, colnames(res), sep = "_")
  res_all <- res
  
  
  # Change columns names to plot data 
  colnames(res) <- c("MeanExp","logFC", "lfcSE", "stat", "pvalue", "padj")
  
  
  
  ##############################################################################
  #                            DATA PROCESSING
  ##############################################################################

  
  ## Results as a data frame
  res_df <- as.data.frame(res) 
  res_df$GeneID <- rownames(res_df)
  print(dim(res_df))
  
  ## Threshold label
  threshold <- paste("alpha_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep ="")
  
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
  res_df <- merge(res_df, gene_names, by = "GeneID") 
  print(head(res_df))
  print(dim(res_df))
  
  ## MERGE WITH GENE COUNTS
  # Row names to a variable
  genes <- gene_counts[, metadata$Sample]
  genes$GeneID <-  rownames(genes)
  # Merge gene_counts and comparison results
  result <- merge(x = res_df, y = genes, by = "GeneID")
  
  
  ## Differential expressed genes
  # Select differentially expressed genes
  df <- result[which(result$DEG == "YES"),]
  # DEDs with the normalized counts values
  df_norm <- merge(df, dds_norm, by = "GeneID")
  
  # Matrix
  # Remove GeneID column
  #----------------------------------------------------------------------------------------------------------------------------------
  m <- dds_vst[df$GeneID, which(gsub(pattern = "VST_", replacement = "", x = colnames(dds_vst)) %in% metadata$Sample)]
  colnames(m) <- metadata$Sample
  #----------------------------------------------------------------------------------------------------------------------------------
  
  ##############################################################################
  #                               PLOT DATA 
  ##############################################################################

  ##############################################################################
  #                               SAVE DATA 
  ##############################################################################
  
  }













################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project_name <- project
log_data$Organism <- specie
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_ord, collapse =",")
log_data$Outliers <- paste(outliers, collapse = ",") 
log_data$Varexp <- paste(var_exp, collapse = ",") 
log_data$filter_cutoff <- filter_cutoff 
log_data$fdr_cutoff <- fdr_cutoff
log_data$lfc_cutoff <- lfc_cutoff
log_data$correction <- correction
log_data$contrast <- paste(unlist(contrast), collapse = ",")
log_data$colortrt <- paste(color_list[[1]], collapse = ",")
log_data$colorheat <- paste(color_list[[2]], collapse = ",")
log_data$colordir<-  paste(color_list[[3]], collapse = ",")
log_data$colorsh <- paste(color_list[[4]], collapse = ",")

write.table(as.data.frame(log_data), paste(path, project, "/1_DEG_v2_DESeq_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")

