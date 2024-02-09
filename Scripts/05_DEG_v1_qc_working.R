################################################################################
#                 DIFFERENTIALLY EXPRESSED GENES SCRIPT
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
# This script is the quality control of the data set. In addition, the gene count 
# matrix is filtered by each comparison so the input data in each comparison made 
# later on is the same. This filtering of the lowly expressed genes is based on the 
# filterByExprss available in EdgeR. The design formula per each analysis is also 
# created here.


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


#-------------------------------------------------------------------------------
#                             REQUIREMENTS
#-------------------------------------------------------------------------------

## Gene count
# -------------
# 
# Matrix with the genes counts, the columns must correspond to the samples names 
# and the row to each gene. 
#
# This is the output of aligning the reads towards a reference genome. If you 
# want to work with transcripts instead of genes, another approach should be 
# followed.

## Metadata or sample information file
# ---------------------------------------
# 
# Must have the columns:
#   - Sample: The samples must be order and the samples names must be the same 
#       as the columns of the gene count matrix. This column must be named Sample
#   - Variable column the name must be specify in the code 
#   - Covariates such as RIN, pv200, Age, ... This variables must be specify in 
#     the code

## Annotation file
# ------------------
# 
# The annotation file must present at least three columns, Ensembl and Symbol gene 
# annotation and the specie. This file can be found in the following path: 
# W:\DATA_shared\Genomes_Rocky\DEG_Annotation\Annotated_Genes_20231121
# 
# This annotation information correspond to GRCh38.p14 and GRCm39 release 110.


## 0_Sample_info_XXXX.log
# -------------------------
# 
# This file will be used in all the scripts, which gives the order to the data 
# and is used in different several comparison and testings.

#-------------------------------------------------------------------------------


# References
#-------------

# Van den Berge, K., Hembach, K. M., Soneson, C., Tiberi, S., Clement, L., Love, 
# M. I., ... & Robinson, M. D. (2019). RNA sequencing data: hitchhiker's guide to 
# expression analysis. Annual Review of Biomedical Data Science, 2, 139-173.
# https://www.annualreviews.org/doi/pdf/10.1146/annurev-biodatasci-072018-021255
#
# Stark, R., Grzelak, M., & Hadfield, J. (2019). RNA sequencing: the teenage 
# years. Nature Reviews Genetics, 20(11), 631-656.
# https://www.nature.com/articles/s41576-019-0150-2 


# Reference Manuals
# ---------------------
# DESeq2: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Limma-voom: https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
# EdgeR:  https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Reference Pipelines
# ---------------------
# DESeq2: https://github.com/hbctraining/DGE_workshop/
# Limma-voom: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# EdgeR: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# Wilcoxon rank-sum test: https://rpubs.com/LiYumei/806213


################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC64"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/user/"
path <- "W:/ulazcano/"

# Date of the log file 0_Sample_info_XXXX.log
logdate <- "20231222"

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
# outliers <- c("AC1", "AC2")
outliers <- NULL


# Filtering lowly expressed genes
#
# Filtering the data can be an option but in most cases lowly expressed genes are 
# filtered. The filtering method is based on the function filterByExpr from edgeR.
# As a result. a unique matrix per each comparison is obtain and given to DESeq2, 
# limma-voom, EdgeR and Wilcoxon. However, there is a primary filtering step in 
# which we remove the genes with 0 counts.
# 
# Two filtering methods are propose based on the data type, which can be clinical or 
# cell line data sets. Both filtering includes two steps:
#
#   A. In cell line data sets, (1) genes with minimum number of counts (min_total = 15; 
#   default) or less are considered lowly expressed genes. (2) Count the number of 
#   samples per gene which overcome minimum number of counts (min_count = 10; default)
#   and if the minimum sample size is not accomplish the gene is deleted. The 
#   minimum sample size is the smallest condition size.
# 
#   B. In clinical data sets (n_large > 30), (1) genes with minimum number of counts 
#   (min_total = 15; default) or less are considered lowly expressed genes. 
#   (2) Count the number of samples per gene which overcome minimum number of counts 
#   (min_count = 10; default) and if the minimum sample size per the proportion 
#   (min_prop = 0.70; default) is not accomplish the gene is deleted. 
# 
# Minimum number of counts per gene
min_total <- 15
# Minimum number of counts per sample
min_count <- 10
# Large sample size 
n_large <- 30
# Proportion
min_prop <- 0.7

### Threshold criteria 

## Significance level
# The significance level is the threshold set to considered an adjusted p-value 
# from a gene to be significant or not. 
# We use the adjusted p-value (or false discovery rate) as the measure to consider 
# the genes significant, instead of using the p-value due to the high dimensionality 
# of our data 

# Options
# fdr_cutoff <- 0.01
fdr_cutoff <- 0.05      # Default option

## Log2 fold change threshold
# Be careful the you have to make sure that the threshold is always refereed as 
# log2 fold change and not fold change 
lfc_cutoff <- log2(1.5)         # Default option

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- "BH"


# Color list
# Option 1: Let the pipeline choose the colors
# Function named *color_palette* selects the colors for the condition
color_list <- list(Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
                   Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
                   Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))
# # Option 2: Select your own colors
# color_list <- list(trt = c(`0` = "#A6DAB0", `4` = "#C18BB7", `24` = "#D7B0B0", `48` = "#8BABD3"),
#                    Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
#                    Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
#                    Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/libraries_degs.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))

# Analysis ID 
# Created using a timestap to trace all de outputs belonging to a certain setup
analysis_ID <- format(Sys.time(), "%Y%m%d%H%M%S")

# Load log file 
logfile <- read.table(paste(path, project, "/log/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/04_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
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
var_exp <- NULL
# <- unlist(strsplit(logfile$covariance, split = ","))

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))


## Generate color list
# Automatically generate the colors for the treatment condition
# Change condition name, to make sure that fits the trt name
if(length(color_list)<4){color_list <- color_palette(color_list, trt, lvl_ord, palette = "Dark2")}else{names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")}

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Data frame with summary information from all the comparisons
# Columns are Comparison, Method, Genes, Upregulated and Downregulated
out_df <- data.frame()


# Create a workbook for pca data
exc_pca <- createWorkbook()


################################################################################
#                               LOAD FILES                       
################################################################################


# Sample information/ Metadata
data_info <- read.csv(file = paste(path, project, "/Sample_info.csv", sep = ""), header = TRUE)

# Gene count matrix
raw_counts <- read.table(file = dir_infiles, sep = "\t",  header = TRUE, stringsAsFactors = TRUE)


print(dim(raw_counts))
print(dim(data_info))

print(head(raw_counts))
print(head(data_info))



################################################################################
#                        SET WORKING DIRECTORY
################################################################################


# Create a folder for the analysis. This folder contains several folders 
# classified in Results and Figures. 

# Create output directory
dir.create(file.path(dir_out,"05_DEG_ANALYSIS"), showWarnings = FALSE)
dir_out <- paste(dir_out,"/05_DEG_ANALYSIS", sep='')
setwd(dir_out)

# Output directory with the ID
dir.create(file.path(dir_out, analysis_ID))
dir_ID <- paste(dir_out, "/", analysis_ID, sep = "")

# Quality control figures folder
dir.create(file.path(dir_ID, "QC"), showWarnings = FALSE)
dir_fig <- paste(dir_ID, "/QC", sep='')

# Results folder
dir.create(file.path(dir_ID, "Results"), showWarnings = FALSE)
dir_output <- paste(dir_ID, "/Results", sep='')


################################################################################
#                           PREPROCESSING 
################################################################################


#### Metadata information #### 
# 
# Create a sample information file with the following information
# - Sample
# - Experimental condition
# - Covariates variables adjusted with z-scores
sample_info <- data.frame(Sample = data_info$Sample)
sample_info[, trt] <- factor(data_info[[trt]], levels = lvl_ord)

# Row names as the sample name
row.names(sample_info) <- sample_info$Sample


#### Remove outliers #### 
if(is.null(outliers) == FALSE){
  sample_info <- sample_info[-which(rownames(sample_info) %in% outliers),]
  raw_counts <- raw_counts[,-which(rownames(raw_counts) %in% outliers)]
  print(paste("Remove outliers:", outliers, sep = " ")); print(ncol(raw_counts))
  project <- paste(project, outliers, sep = "_", collapse = "_")
}else{print("No sample was considered an outlier.")}


#### Add covariates effects #### 
if(!is.null(var_exp)){
  # Normal for EdgeR and limma-voom
  sample_info[, var_exp] <- data_info[, which(colnames(data_info) %in% var_exp)]
  # Z_score for DESeq2
  for (k in 1:length(var_exp)){
    if(length(unique(sample_info[, var_exp[k]]))>1){
      sample_info[, paste(var_exp[k], "_zscore", sep = "")] <- as.vector(scale(data_info[[var_exp[k]]], center = TRUE, scale = TRUE))}}
} else{print("No variation correction")}



### Verify the order of the columns are the same as the sample information file #### 
raw_counts <- raw_counts[, match(sample_info$Sample, colnames(raw_counts))]

if(sum(rownames(sample_info) == colnames(raw_counts)) != ncol(raw_counts)){
  cat("ERROR: Samples in count matrix are not ordered", "REORDER COUNT MATRIX COLUMNS", sep = "\n")
}else{print("Samples in count matrix are ordered")
  raw_counts <- raw_counts[, match(sample_info$Sample, colnames(raw_counts))]
}


#### Filter #### 
# Number of samples
n <- ncol(raw_counts)

# Filter genes with no counts based on the number of samples
# This is a common filtering step for all the comparison levels
raw_counts <-  raw_counts[(rowSums(raw_counts) != 0), ]

## Statistical summary
cat(paste("Total genes in raw count:", dim(raw_counts)[1]))
#cat(paste("Total genes in gene count:", dim(gene_counts)[1]))
print("There is NA data?")
print(is.na(raw_counts) %>% table())


## Barplot verification the quality of the gene counts
bp <- data.frame(Samples = colnames(raw_counts),
                 trt = sample_info[trt],
                 Values = as.vector(colSums(raw_counts)))

ggplot(bp, aes(x = factor(Samples, levels = sample_info$Sample), y = Values, fill = !!as.name(trt)))+
  geom_bar(stat = "identity")+
  labs(y = "Counts", x = "")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1, size = 5), legend.position = "none")
ggsave(filename = paste("00_Barplot_rawcounts_", project, "_", analysis_ID, ".pdf", sep =""), height = 4, width = 4, plot = last_plot(), path = dir_fig, bg = "white")


################################################################################
#                               COMPARISON 
################################################################################


# DESeq2 design formula
# Used to estimate the variance stabilization (VST) method proposed in DESeq2
# 
# Steps
# 1. Check there are covariates
# 2. Check the values of the covariate are not equal to avoid colinearity in 
#   the model. If their values are equal the variable is not included in the 
#   model.
# 3. Create the design formula
design_cond <- design_condition("DESeq2", trt, var_exp, metadata = sample_info)
print(design_cond)



################################################################################
#                           Data transformation
################################################################################
# 
# The methods for the differentially expressed genes used are DESeq2, EdgeR, 
# limma-voom and Wilcoxon test
# 
# Two analysis paths. 
#   1. DESeq2 analysis with NO FILERED data to try to maximize DESeq2 potential.
#   2. FILTERED data for DESeq2, EdgeR and limmavoom (Wilcoxon)
# 
# Normalized and filtered gene counts for the DEG analysis are performed in 
# following steps. 


for (i in 1:2){
  if(i == 1){
    filter_lab <- "nofilter"
    gene_counts <- raw_counts
  }else{
    # Filter gene count matrix per each comparison
    #
    # This data will be used when running DESeq2, EdgeR, limma-voom and Wilcoxon to
    # perform the test for the differentially expressed genes
    keep <- filter_genecounts(raw_counts, sample_info, trt, min_count = min_count, min_prop = min_prop, n_large = n_large, min_total = min_total)
    gene_counts <- raw_counts[keep,]
    filter_lab <- "filtered"
  }
  
  print(filter_lab)
  print(dim(gene_counts))
  
  # Final output
  df <- gene_counts
  df$Ensembl <- rownames(df)
  
  
  ##############################################################################
  #                                 DESeq2
  ##############################################################################
  
  
  # Step 1: Create DESeq object
  # ----------------------------------------------------------------------------
  #
  # DESeqDataSetFromMatrix has all the information necessary to run DESeq. 
  # - Count matrix
  # - sample_info with the following columns: condition and covariates. Data 
  #     associated to the samples valuable to perform the contrast. The row names
  #     must be the samples names
  # - Design formula with the experimental conditions and the covariates 
  dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = sample_info, design =  eval(parse(text = design_cond)))
  
  
  # Step 2: Estimate size factor
  # ----------------------------------------------------------------------------
  dds <- estimateSizeFactors(dds)
  
  
  # Step 3: Variance stabilizing transformation
  # ----------------------------------------------------------------------------
  # 
  # Variance stabilization methods in log2 scale to interpret the data
  # 
  # Choose VST for samples size group smaller than 30. Why? 
  # 
  # if you have many samples (e.g. 100s), the rlog function might take too long, and so the vst function will be 
  # a faster choice. The rlog and VST have similar properties, but the rlog requires fitting a shrinkage term for 
  # each sample and each gene which takes time. See the DESeq2 paper for more discussion on the differences 
  # (Love, Huber, and Anders 2014). The rlog function is not as sensitive as vst to the size factors, which can be 
  # an issue when size factors vary widely.
  
  # Why blind = TRUE
  # 
  # blind, for whether the transformation should be blind to the sample information specified by the design formula. 
  # When blind equals TRUE (the default), the functions will re-estimate the dispersions using only an intercept. 
  # This setting should be used in order to compare samples in a manner wholly unbiased by the information about 
  # experimental groups, for example to perform sample quality assurance.
  
  # Output
  # - Columns are samples
  # - Rows are genes
  # 
  #                    N_3_E1    N_2_E4    N_2_E3   N_2_E2   H4_2_E4   H4_2_E3   H4_2_E2   H4_2_E1
  # ENSG00000000003 11.075692 10.993586 11.209109 11.03978 11.136152 11.147595 11.118555 11.125602
  # ENSG00000000419 11.959761 11.938768 12.033994 11.84389 11.846617 11.951190 11.951552 11.888236
  # ENSG00000000457 10.476816 10.557566 10.470185 10.47396 10.568444 10.488905 10.569326 10.394093
  
  # Why VST or RLOG?
  # 
  # There are two different methods to stabilize the variance in DESeq2, vst and rlog. 
  #   The variance stabilization transformation (vst)
  #   The regularized log transformation (rlog) inherently accounts for differences in sequencing depth.
  # 
  # We decided to exclusively use VST, because the reason is given is due to computing time as you can see
  # Note on running time: if you have many samples (e.g. 100s), the rlog function might take too long, and so the vst
  # function will be a faster choice. The rlog and VST have similar properties, but the rlog requires fitting a shrinkage
  # term for each sample and each gene which takes time. See the DESeq2 paper for more discussion on the differences
  # (Love, Huber, and Anders 2014).
  # 
  # The point of these two transformations, the VST and the rlog is to remove the dependence of the variance on the mean, 
  # particularly the high variance of the logarithm of count data when the mean is low. Both VST and rlog use the 
  # experiment-wide trend of variance over mean, in order to transform the data to remove the experiment-wide trend. 
  # Note that we do not require or desire that all the genes have exactly the same variance after transformation. Indeed, 
  # in a figure below, you will see that after the transformations the genes with the same mean do not have exactly the 
  # same standard deviations, but that the experiment-wide trend has flattened. It is those genes with row variance 
  # above the trend which will allow us to cluster samples into interesting groups.
  
  # Estimate the biggest group sample size
  group_n <- max(as.vector(tabulate(sample_info[[trt]])))
  
  m_blindTRUE <- as.data.frame(assay(vst(dds, blind = TRUE)))
  m_blindFALSE <- as.data.frame(assay(vst(dds, blind = FALSE)))
  vsd_type <- "VST"
  
  
  # Step 4: Normalization
  # ----------------------------------------------------------------------------
  # 
  # The normalization method is DESeq2's median of ratios in which counts divided 
  # by sample-specific size factors determined by median ratio of gene counts 
  # relative to geometric mean per gene. 
  # 
  # Accounting for sequencing depth and RNA composition
  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  
  
  # Step 5: Stack matrix
  # ----------------------------------------------------------------------------
  # Matrix to plot the data
  m_blindTRUE_s <- stack(m_blindTRUE)
  m_blindTRUE_s_mod <- m_blindTRUE_s
  m_blindTRUE_s_mod[trt] <- rep(sample_info[[trt]], each = dim(gene_counts)[1])


  ##############################################################################
  #                               QC Plots 
  ##############################################################################
  
  
  ## COUNTS DISTRIBUTION PLOTS
  # Distribution of gene read counts per each sample
  ggplot(data = m_blindTRUE_s, aes(x = values, fill = ind))+
    geom_histogram(stat = "bin", bins = 200, position="identity", alpha = 0.3, show.legend = FALSE)+
    labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")
  ggsave(filename = paste("00_Dist_genecounts_samples_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig, bg = "white")
  
  # Distribution of gene read counts with grid per each sample
  ggplot(data = m_blindTRUE_s, aes(x = values, fill = ind))+
    geom_histogram(stat = "bin", bins = 200, show.legend = FALSE)+
    facet_wrap(~ ind, ncol = 4)+
    labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")+
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1))
  ggsave(filename = paste("00_Dist_genecounts_samples_grid_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep=""), height = 8, width = 6, plot = last_plot(), path = dir_fig, bg = "white")
  
  
  ## BOXPLOT
  ggplot(data = m_blindTRUE_s_mod, aes(x = ind, y = values, fill = !!as.name(trt)))+
    geom_boxplot(linewidth = 0.09, outlier.size = 0.5)+
    labs(x = "", y = "Log2(Counts)")+
    scale_fill_manual(values = color_list[[trt]])+
    theme(legend.position = "none", text = element_text(size = 5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), panel.grid = element_blank())
  ggsave(filename = paste("01_Distribution_boxp_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), path = dir_fig, height = 4, width = 6, plot = last_plot(), bg = "white")
  
  ## DENSITY PLOT
  # A few things to note about this figure is that there is a huge peak at exactly
  # -3.321928 which can be ignored because this is value is calculated from log2(0.1)
  # explained in the box above. This peak consists of all 0-values (inactive genes)
  # which isn’t very important to us now.
  pdf(paste(dir_fig, "/02_Density_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), height = 5, width = 5, bg = "white")
  plotDensity(m_blindTRUE, col = rep(color_list[[trt]], each = 4), xlab = "Log2(Counts)", ylab = "Density", main = "Expression Distribution")
  dev.off()
  
  
  ## MODELING COUNT DATA
  # Perform with raw gene gene counts 
  model_m <- data.frame(Mean = apply(gene_counts, 1, mean), Variance = apply(gene_counts, 1, var))
  head(model_m)
  
  ggplot(data = model_m, aes(x = Mean, y = Variance))+
    geom_point(alpha = 0.20, colour = "#696969")+
    geom_line(aes(x = Mean, y = Mean), color = "red", show.legend = FALSE)+
    scale_y_log10()+
    scale_x_log10()+
    labs(x = "Mean gene expression level (log10 scale)", y = " Gene-level variance (log10 scale)", title = "")
  ggsave(filename = paste("03_NB_fit_data_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig, bg = "white")
  
  
  ##  CORRELATION 
  #
  # Execute a Pearson correlation which accept possible NAs.
  # The acceptance of the NA should be valuable in the future in case, we accept
  # comparisons when a sample is missing a PSI value.
  #
  # The results of the correlation matrix must be aligned with the results
  # in the heatmap and PCA.
  pem <- cor(m_blindTRUE, method = "pearson", use = "na.or.complete")
  plot_h <- pheatmap(pem, color = colorRampPalette(brewer.pal(9, "Blues"))(255),
                     cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 6, fontsize_col = 6, border_color = NA, treeheight_row = 0, treeheight_col = 0)
  
  pdf(paste(dir_fig, "/05_Cor_pearson_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_h)
  dev.off()
  
  
  ## PCA PLOTS
  plot_pcas <- pca_plot(m_blindTRUE, trt, sample_info, color_list)
  
  ggsave(filename = paste("PCA_params_genecounts_", project, "_", filter_lab, "_scree", "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[1]], path = dir_fig, height = 4, width = 4, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[4]], path = dir_fig, height = 5, width = 6, bg = "white")
  addWorksheet(exc_pca, paste("R", filter_lab, sep = "_"))
  addWorksheet(exc_pca, paste("GenesPCs", filter_lab, sep = "_"))
  writeData(exc_pca, as.data.frame(plot_pcas[[5]]), sheet = paste("R", filter_lab, sep = "_"))
  writeData(exc_pca, as.data.frame(plot_pcas[[6]]), sheet = paste("GenesPCs", filter_lab, sep = "_"))
  
  
  ## HEATMAP
  plot_heatmap <- heatmap_plot(m_blindTRUE, sample_info, trt, color_list)
  
  pdf(paste(dir_fig, "/Heatmap_zscore_genecounts_", project, "_", filter_lab, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_heatmap)
  dev.off()
  
  
  ##############################################################################
  #                               EdgeR/limma-voom
  ##############################################################################
  
  if(i == 2){
    
    # Step 1: Create DGEList object
    # ----------------------------------------------------------------------------
    deg <- DGEList(counts = gene_counts, group = sample_info[,trt])
    if(!is.null(var_exp)){deg$samples[var_exp] <- sample_info[,var_exp]}
    
    
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
    
    
    # Step 3: Counts per million normalization method 
    # ----------------------------------------------------------------------------
    # 
    # In the counts per million method, counts are scaled by total number of reads
    # 
    # Accounting for sequencing depth
    cpm_counts <- cpm(deg, log = FALSE, normalized.lib.sizes = TRUE)
    cpm_counts <- as.data.frame(cpm_counts)
    
    
    
    ################################################################################
    #                               OUTPUT DATA           
    ################################################################################
    
    
    colnames(cpm_counts) <- paste("CPM_", colnames(cpm_counts), sep = "")
    df_cpm <- cbind(df, cpm_counts)
    }
  
  colnames(m_blindTRUE) <- paste("VST_", colnames(m_blindTRUE), sep = "")
  colnames(m_blindFALSE) <- paste("VST_", colnames(m_blindFALSE), sep = "")
  colnames(norm_counts) <- paste("Norm_", colnames(norm_counts), sep = "")
  
  counts_blindTRUE <- cbind(df, m_blindTRUE)
  counts_blindFALSE <- cbind(df, m_blindFALSE)
  
  counts_blindTRUE <- cbind(counts_blindTRUE, norm_counts)
  counts_blindFALSE <- cbind(counts_blindFALSE, norm_counts)
  
  
  ################################################################################
  #                               SAVE DATA           
  ################################################################################
  
  # Save DESeq
  write.table(counts_blindTRUE, paste(dir_output,"/GeneCount_", vsd_type , "_blindTRUE_", project, "_", filter_lab, "_", analysis_ID, ".txt", sep = ""))
  write.table(counts_blindFALSE, paste(dir_output,"/GeneCount_", vsd_type , "_blindFALSE_", project, "_", filter_lab, "_", analysis_ID, ".txt", sep = ""))
}

# Save EdgeR/limma-voom
write.table(df_cpm, paste(dir_output,"/GeneCount_CPM_", project, "_", filter_lab, "_", analysis_ID, ".txt", sep = ""))

# Save metadata information
write.table(sample_info, paste(dir_output, "/Metadata_", project, "_", analysis_ID, ".txt", sep = ""))

# Save PCA data
saveWorkbook(exc_pca, file =  paste(dir_output, "/PCA_data_QC_", project, ";", vsd_type, "blindTRUE_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)


################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
log_data <- c()
log_data$analysis_ID <- analysis_ID
log_data$Date <- Sys.time()
log_data$project_name <- project
log_data$Organism <- logfile$Organism
log_data$dir_out <- dir_out
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_ord, collapse =",")
log_data$Outliers <- paste(outliers, collapse = ",") 
log_data$Varexp <- paste(var_exp, collapse = ",") 
log_data$filter_lab<-filter_lab
log_data$min_count <- min_count
log_data$min_total <- min_total
log_data$n_large <- n_large
log_data$min_prop <- min_prop
log_data$fdr_cutoff <- fdr_cutoff
log_data$lfc_cutoff <- lfc_cutoff
log_data$correction <- correction
log_data$Variance <- vsd_type
log_data$contrast <- logfile$contrast
log_data$colortrt <- paste(color_list[[trt]], collapse = ",")
log_data$colorheat <- paste(color_list[["Heatmap"]], collapse = ",")
log_data$colordir<-  paste(color_list[["Direction"]], collapse = ",")
log_data$colorsh <- paste(color_list[["Shared"]], collapse = ",")


write.table(as.data.frame(log_data), paste(path, project, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), row.names = FALSE, eol = "\r")

print("DEG Analysis ID")
print(analysis_ID)


print("Quality Control Analysis completed!")

################################################################################
#                             FOLLOW THE ANALYSIS           
################################################################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# source(paste(path, project, "/Scripts/05_DEG_v2_DESeq2.R", sep = ""))
# source(paste(path, project, "/Scripts/05_DEG_v2_EdgeR.R", sep = ""))
# source(paste(path, project, "/Scripts/05_DEG_v2_limma.R", sep = ""))
# source(paste(path, project, "/Script/05_DEG_v2_wilcoxon.R", sep = ""))          #ONLY WITH BIG DATA SETS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

