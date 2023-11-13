################################################################################
#                               RUN EdgeR
################################################################################

# Summary
# ---------- 
# 
# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory.



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
logdate <- "20231110"
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
analysis <- "EgdeR"



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
  dir_outfolder <- paste(dir_out, "/04_DEG_ANALYSIS/", name, sep='')
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
  metadata[,trt] <- factor(metadata[,trt])
  
  
  # Comparison levels
  comp_lvl <- c(contrast[[i]][2],  contrast[[i]][3])
  
  # Select the contrast levels
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% comp_lvl)]
  
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
  deg$samples[var_exp] <- metadata[,var_exp]
  
  
  # Step 2: Normalization
  # ----------------------------------------------------------------------------
  # 
  # This step is common between EdgeR and limma-voom 
  # 
  # Minimized the logFC between samples for most genes using weighted Trimmed Mean of M-values 
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
  
  
  # Step 4: Contrast gene
  # ----------------------------------------------------------------------------
  
  # Perform contrast
  fit <- exactTest(degde, pair = c(control, experimental), dispersion = "auto", rejection.region = "doubletail")
  
  # MA plot 
  # In the x axis, average logCPM(counts per million) but in DESeq2 is mean of normalized counts
  pdf(file = paste(dir_fig, "/00_MA_plot_", analysis, "_", project,".pdf", sep = ""), height = 5 , width = 5)  
  plotSmear(fit, de.tags = labels_sig, ylab = "log fold change"); abline(h = 0, col  ="black", lty = 1, lwd = 2)
  dev.off()
  
  
  
  
  
  
  # Estimate the adjusted p-value 
  res <- as.data.frame(topTags(fit, adjust.method = correction, n = dim(gene_counts)[1]))
  res <- res[match(rownames(gene_counts), rownames(res)),]
  
  # Change columns names to plot data 
  colnames(res) <- c("logFC", "logCPM", "pvalue", "padj")
  
 
  # DEGs list
  labels_sig <- rownames(res[which(res$padj <= fdr_cutoff & abs(res$logFC) > lfc_cutoff),])
  
 
  
  
  ##############################################################################
  #                            Data Processing
  ##############################################################################
  
  
  
  ##############################################################################
  #                                 Plot
  ##############################################################################
  

  
  ##############################################################################
  #                               Save data 
  ##############################################################################
  
  
  
  # All results
  colnames(res_log2) <- paste(mt, colnames(res_log2), sep = "_")
  data <- cbind(res_df, res_log2)
  write.table(data, paste(ref, ";All_", md, "blindFALSE", "_pval", threshold,".txt", sep = ""))
  write.xlsx(data, paste(ref, ";All_", md, "blindFALSE", "_pval", threshold,".xlsx", sep = ""), overwrite = TRUE)
  
  # Differential expressed genes
  colnames(m) <- paste(mt, colnames(m), sep = "_")
  sel <- cbind(res_df, m)
  write.table(data, paste(ref, ";DEGs_", md, "blindFALSE", threshold,".txt", sep = ""))
  write.xlsx(data, paste(ref, ";DEGs_", md, "blindFALSE", threshold,".xlsx", sep = ""), overwrite = TRUE)
  
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

write.table(as.data.frame(log_data), paste(path, project, "/1_DEG_v2_DESeq_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")

