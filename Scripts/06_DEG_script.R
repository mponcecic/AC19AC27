################################################################################
#                 DIFFERENTIALLY EXPRESSED GENES SCRIPT
################################################################################


# Summary 
#------------

# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory. 

# The method applied to estimate the DEGs are DESeq2, EdgeR, limma-voom
# and Wilcoxon rank-sum test. 
# 
# First, you will find a chunk of code where you can adjust all the parameters, 
# filtering, thresholds and colors needed to perform in this analysis. This 
# variables are deeply detailed in the code. Following the code, it's time to 
# load the data sets, be careful everything is matched.
# 
# Second,
# 
# Third,  
# 
# Fourth,



# Code chunk contain between the following separators can be modify and must be 
# modify to perform the correct analysis. Two chunk will be found one referring
# to the parameter selection and the order for loading the data.


# DESeq2 recommendations include all the comparison group in the DESeq object 
# and in the exploratory analysis of the data

# Does not consider interactions in the design formula, but can be added to the script

# Authors: Ivana Rondón and María Ponce


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
# The annotation file must present at least two columns being the first one the 
# Ensembl gene annotation and the second the Symbol annotation. 

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
#                         LOAD LIBRARIES
################################################################################


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(require(openxlsx))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(factoextra)) 
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(affy))        # PlotDensity
suppressPackageStartupMessages(library(dendsort))    # dist


################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC58"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Date of the log file
logdate <- ""
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load log file 
logfile <- read.table(paste(path, project_name, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/04_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
# dir_out <- paste(path, project, sep = "")   # Default option
dir_out <- paste(path, project, "/05_DEGs/", sep = "")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- logfile$condition

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- logfile$condition_order


# Contrast
# If contrast set NULL the different contrasts will be created based on the 
# lvl_order vector. Nonetheless, you can add the contrast manually following the 
# 
# DESeq2 formula:
# contrast <- list(c("Experimental Variable", "Experimental level", "Reference level"), 
#                  c("Experimental Variable", "Experimental level", "Reference level"))
contrast <- NULL


# Variance sources to include in the model
# Can be set NULL
# Mouse: RIN
# Human: RIN, dv200, Age, ... 
#
# Options
# var_exp <- c("Age", "dv200")
# var_exp <- NULL
var_exp <- c("RIN")


### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
# outliers <- c("AC1", "AC2")
outliers <- NULL

# Filtering lowly expressed genes
#
# Filtering the data can be an option but in most cases lowly expressed genes are 
# filtered. Two methods are present here to filter these genes.
#
#   A. Genes with 10 counts or less are considered lowly expressed genes and are 
#     excluded in experiments with less than 30 samples.
# 
#   B. In data sets with high dimensionality (n>30), genes with 50% of the samples 
#     with zero counts are removed. 
# 
# If filtering is set FALSE, genes with 0 counts will be removed.
# 
# This step will be performed for limma-voom and EdgeR independently.
# 
# Options:
# filter_cutoff <- FALSE
filter_cutoff <- TRUE


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
# log2 fold change and not fold change <------
lfc_cutoff <- log2(1.5)         # Default option

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- "BH"


## Data transformation considering the metadata information or not, by default 
# this is set FALSE
# blind = TRUE, unbiased why to approach data without considering the experimental 
# design 
# blind = FALSE, considers the experimental design and differences in the counts 
# might be attribute to the experimental conditions
# 
# Options:
# blind <- TRUE
blind <- FALSE        # Default option

## Shrinkage estimators
# Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq run 
# without LFC shrinkage
# 
# Options: 
# shk <- "ashr"
# shk <- "apelgm"
# shk <- "normal" 
shk <- NULL

# Color list
# Better to choose colors manual but there is a function named *color_palette* 
# which selects the colors for the condition

# Optional:
# color_list <- list(Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
#                    Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
#                    Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))
color_list <- list(trt = c(Control = "#A6DAB0", `4` = "#C18BB7", `24` = "#D7B0B0", `48` = "#8BABD3"), 
                   Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
                   Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
                   Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Data frame with summary information from all the comparisons
# Columns are Comparison, Method, Genes, Upregulated and Downregulated
out_df <- data.frame()

# Specie 
specie <- logfile$Organism


################################################################################
#                               LOAD FILES                       
################################################################################

# Reference genome
ref_genome <- read.table(paste(path, "/DEG_annotation/gene_annotation_20231106.txt", sep = ""),
                     col.names = c("Symbol", "EnsemblID", "Organism"))
genome <- ref_genome[which(ref_genome$Organism == specie), ]

# Sample information/ Metadata
data_info <- read.csv(file = paste(path, project, "/Sample_info.csv", sep = ""), header = TRUE)

# Gene count matrix
raw_counts <- read.table(file = dir_infiles, sep = "\t",  header = TRUE, stringsAsFactors = TRUE)



print(dim(raw_counts))
print(dim(data_info))
print(dim(genome))

print(head(raw_counts))
print(head(data_info))
print(head(genome))


################################################################
#                     FUNCTIONS
################################################################


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



################################################################################
#                     SET WORKING DIRECTORY
################################################################################


# Create a folder for the analysis. This folder contains several folders 
# classified in Results and Figures. 

# Create output directory
dir.create(file.path(dir_out,"04_DEG_ANALYSIS"))
dir_out <- paste(dir_out,"04_DEG_ANALYSIS", sep='')
setwd(dir_out)

# <- Filtered results ->
# Files folder
dir.create(file.path(dir_out,"Results"), showWarnings = FALSE)
dir_raw <- paste(dir_out,"/Results", sep='')
# Figures folder
dir.create(file.path(dir_out, "Figures"), showWarnings = FALSE)
dir_fig_raw <- paste(dir_out, "/Figures", sep='')



################################################################################
#                           PREPROCESSING 
################################################################################


#### Create contrast list #### 
if(is.null(contrast)== TRUE){contrast <- create_contrast(trt, lvl_ord)}


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
}else{print("No sample was considered an outlier.")}


#### Add covariates effects #### 
if(is.null(var_exp) == FALSE){
  # Normal for EdgeR and limma-voom
  sample_info[, var_exp] <- data_info[, which(colnames(data_info) %in% var_exp)]
  # Z_score for DESeq2
  sample_info[, paste(var_exp, "_zscore", sep = "")] <- as.vector(scale(data_info[, which(colnames(data_info) %in% var_exp)], center = TRUE, scale = TRUE))
} else{print("No variation correction")}


### Verify the order of the columns are the same as the sample information file #### 
gene_counts <- raw_counts[, match(sample_info$Sample, colnames(raw_counts))]

if(sum(rownames(sample_info) == colnames(gene_counts)) != ncol(gene_counts)){
  cat("ERROR: Samples in count matrix are not ordered", "REORDER COUNT MATRIX COLUMNS", sep = "\n")
}else(print("Samples in count matrix are ordered"))
### Add modification and verify que esta todo bien 

#### Filter #### 
# Number of samples
n <- ncol(gene_counts)

# Filter lowly expressed genes based on the number of samples 
if(filter_cutoff == TRUE){
  if(n <= 30){
    # Filter genes with 10 count or less
    gene_counts <- gene_counts[which(rowSums(gene_counts) >= 10),]
    label <- "filtering_10"
  } else {
    # Filter genes with 0 counts in the 50% of the samples
    gene_counts <- gene_counts[which(rowSums(gene_counts > 0) >= round(0.5*n)),]
    label <- "filtering_50"
  }
  # Filtering set FALSE, we only filter the gene with 0 counts
} else {
  # Filter genes with 10 count or less
  gene_counts <- gene_counts[which(rowSums(gene_counts) != 0),]
  label <- "filtering_0"}


a <- c("1", "1", "1", rep("2", 5))

which(table(a)>3)
which(table(sample_info[[trt]])>30)

print(label)
print(paste("Total genes in raw count: \n", dim(raw_counts)[1]))
print(paste("Total genes in gene count: \n", dim(gene_counts)[1]))



#### Annotation #### 

# Annotate all the genes in the matrix with the Symbol annotation

# CAUTION: Symbol id present more than one Ensembl identifier. 
# I propose to perform the analysis using the Ensembl identifier, if not the 
# mean value of the genes with the same Symbol id should be performed.

annot <- genome[which(genome$GeneID %in% rownames(gene_counts)),]
gene_names <- data.frame(GeneID = annot$GeneID, Symbol = annot$Symbol)


#### Design formula ####
#
# Two different design formulas are created for DESeq2 and EdgeR/limma-voom
# DESeq2 design formula
design_cond = ifelse(is.null(var_exp) == TRUE, paste("~", trt, sep = " "), paste("~", paste(var_exp, "_zscore", " +", sep = "", collapse = " "), trt, sep = " "))

# EdgeR and limma-voom 
# Build a design/model matrix 
design_cond2 = ifelse(is.null(var_exp) == TRUE, paste("~ 0 + ", trt, sep = " "), paste("~ 0 +", paste(var_exp, "+", sep = " ", collapse = ""), trt, sep = " "))
m_model = model.matrix(as.formula(design_cond2), sample_info)




################################################################################
#                             DESeq2
################################################################################


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
                              colData = sample_info,                       # Input: Metadata
                              design = eval(parse(text = design_cond)),    # Design matrix
                              tidy = FALSE,                                # Default; Option: TRUE = First column of count data is the row names for the count matrix
                              ignoreRank = FALSE                           # Default: Reserved for DEXSeq developers
)


# Step 2: Run DESeq
# ------------------------------------------------------------------------------
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
pdf(file = paste(dir_fig_raw, "/Dispersion_DESeq2_", project, ".pdf", sep = ""), width = 6, height = 5)
plotDispEsts(dds,
             # ymin,                          # Lower bound for points on the plot, points beyond == triangles 
             genecol = "black",               # Gene-wise dispersion estimates
             fitcol = "red",                  # Fitted estimates
             finalcol = "dodgerblue",         # Final estimates used for testing
             legend = TRUE, 
             log = "xy",                      
             cex = 0.45,                      # Change the size of symbols/text relative to the default size
             # type = "",                     # type: p, l, b, c, ... 
             # main = "",
             # xlab = "", 
             # ylab = ""
)
dev.off()


# Step 3: Normalization 
# ------------------------------------------------------------------------------
# 
# This step is also performed using EdgeR and limma-voom but the method is 
# different. Here we use the mean ratios method and in the other cases they use
# TMM nethod (Trimmed Mean of M-values).
# 
# The median of ratios is obtained by counts divided by sample-specific size 
# factors determined by median ratio of gene counts relative to geometric mean 
# per gene. 
# This method corrects the sequencing detph and the RNA composition. 
# Can be used for DE analysis but not for within sample comparisons

dds <- estimateSizeFactors(dds)


# Step 4: Transform reads counts and aggregate 
# ------------------------------------------------------------------------------
#
# For the future exploratory analysis, we must transform the data to interpretate 
# the data. First, we normalized the count matrix and afterwards we scale the data. 
# 
# There are several options to scale the data, DESeq2 offers rlog() and vst() options, 
# some people also use rlog2(count +1)
#   - The rlog() is log2(count data) which minimizes the differences between the samples for 
# rows with small counts. 
#   - The vst() is the variance stabilization method from the fitted dispersion-mean 
# relation(s) and then transforms the count data (normalized by division by the size factors 
# or normalization factors)

# Normalized data and scale data 
dds_norm <- log(counts(dds, normalized = TRUE) + 1, base = 2)
dds_norm <- data.frame(dds_norm)
dds_norm$GeneID <- rownames(dds_norm)
print(dim(dds_norm))


if(blind == FALSE){
  # VST transformation
  dds_vst <- assay(vst(dds, blind = FALSE))
  # rlog transformation
  dds_rlog <- assay(rlog(dds, blind = FALSE))
} else {
  # VST transformation
  dds_vst <- assay(vst(dds, blind = TRUE))
  # rlog transformation
  dds_rlog <- assay(rlog(dds, blind = TRUE))
}

# log 2 transformation
dds_log2 <- log(assay(dds), base = 2)


#############################################################################################################
# Alternative to regularized log2 when the sample size is bigger.
# rlog takes longer than vst, that's the main reason why tis will be an option

# if(ncol(gene_counts)>=30){
#   print("Regularized log2 transformation")
#   m <- assay(rlog(dds, blind = blind))
# } else {
#   print("Variance stabilizing transformation")
#     m <- assay(vst(dds, blind = blind))}
#############################################################################################################



# Step 5: Annotated and aggregate data 
# ------------------------------------------------------------------------------
#
# Some genes are presented in different regions of the genome and present an ENSMBL 
# identifier unique per each position. Nonetheless, the gene (Symbol) is refereed 
# as the same being the reason the data must be aggregated together. 
#
# The mean is performed between gene counts with different ENSMBL identifiers 
# shared by the same Symbol identifier.
# 
# The annotation is performed using the gene information file from the the genome 
# used to mapped the data. 

# Annotated the data by merging the normalized counts with the gene info
dds_norm_mod <- merge(dds_norm, gene_names, by = "GeneID")

# Aggregate the based on Symbol id and performed mean 
dds_norm_mod1 <- aggregate(dds_norm_mod, by = list(dds_norm_mod$Symbol), mean)
print(dim(dds_norm))

# Symbol name as row names
rownames(dds_norm_mod1) <- dds_norm_mod1$Group.1

# Remove unnecessary columns 
dds_norm_mod1 <- dds_norm_mod1[, -which(colnames(dds_norm_mod1) %in% c("Group.1", "Symbol", "GeneID"))]



################################################################################
#                             EdgeR
################################################################################


# Step 1: Create DGEList object
# ----------------------------------------------------------------------------
#
# DGEList contains all the information necessary to run EdgeR and limma-voom 
# - Count matrix
# - Metadata with the following columns: Sample, Time and RIN. Data associated 
#       to the samples valuable to perform the contrast

# The options of the SummarizedExperiment can be applied in this object.
deg <- DGEList(counts = gene_counts, group = sample_info[,trt])
deg$samples[var_exp] <- sample_info[,var_exp]


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

# Biological coefficient variation per each gene 
pdf(file = paste(dir_fig_raw, "/Dispersion_BCV_EdgeR_", project,".pdf", sep = ""), height = 5 , width = 5)
plotBCV(degde)
dev.off()

# MA plot 
pdf(file = paste(dir_fig_raw, "/MA_plot_EdgeR_", project,".pdf", sep = ""), height = 5 , width = 5)  
plotMD(degde, column = 1); abline(h = 0, col  ="red", lty = 2, lwd = 2)
dev.off()



################################################################################
#                             limma-voom
################################################################################


# Step 3: Run limma-voom 
# ----------------------------------------------------------------------------
#
# Run voom transformation 
# Transform count data to log2-counts per million (logCPM), estimate the mean-variance 
# relationship and use this to compute appropriate observation-level weights. The 
# data are then ready for linear modelling.


## Voom transformation 
dlim <- voom(deg, m_model, plot = FALSE, save.plot = TRUE)

# MA plot
plot_dt <- data.frame(x = as.vector(dlim$voom.xy[[1]]), y = as.vector(dlim$voom.xy[[2]]))
ggplot(plot_dt, aes(x = x, y = y))+
  geom_point()+
  geom_smooth(se = FALSE, color = "red")+
  labs(x = dlim$voom.xy[[3]], y = dlim$voom.xy[[4]])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste(dir_fig_raw, "/MA_plot_limma-voom_", project,".pdf", sep = ""), height = 4, width = 6, plot = last_plot())

# Model fit 
dlim <- lmFit(dlim, m_model, ndups = NULL, spacing = NULL, block = NULL, correlation, weights = NULL, method = "ls")



################################################################################
#                           Output matrix
################################################################################


# This matrix contains:
# - Filtered gene counts with Symbol and Ensembl annotation
# - Normalized and scaled counts
# - Variance stabilized transformation
# - Regularized log2 transformation 

# Normalization and log2(counts+1)
dds_norm1 <- dds_norm

# Rename columns
colnames(dds_norm1) <- paste("Normalized", colnames(dds_norm1), sep = "_")
colnames(dds_vst) <- paste("VST", colnames(dds_vst), sep = "_")
colnames(dds_rlog) <- paste("rlog", colnames(dds_rlog), sep = "_")
colnames(dds_log2) <- paste("log2", colnames(dds_log2), sep = "_")

# Merge all together
data <- cbind(gene_counts, dds_norm1 ,dds_vst, dds_rlog, dds_log2)
data$GeneID <- rownames(data)
data <- merge(data, gene_names, by = "GeneID")
data <- data  %>% select(c("GeneID", "Symbol", everything()))


#################################################################################################################################
#                             EXPLORATORY ANALYSIS
#################################################################################################################################


################################################################################
#                             RAW COUNTS                       
################################################################################

# Stack raw count matrix
raw_st <- stack(raw_counts)


## COUNTS DISTRIBUTION PLOTS
# Distribution of raw read counts, after filtering, per samples
ggplot(data = raw_st, aes(x = values, fill = ind))+
  geom_histogram(stat = "bin", bins = 200, position="identity", alpha = 0.3, show.legend = FALSE)+
  xlim(-10, 2000)+
  labs(x = "Counts", y = "Genes", title = "Raw Expression")
ggsave(filename = paste("00_Dist_rawcounts_samples_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)

# Distribution of raw read counts, after filtering, per each samples
ggplot(data = raw_st, aes(x = values, fill = ind))+
  geom_histogram(stat = "bin", bins = 200, show.legend = FALSE)+
  facet_wrap(~ ind, ncol = 4)+
  xlim(-10, 2000)+
  labs(x = "Counts", y = "Genes", title = "Raw Expression")+
  theme(strip.background = element_blank())
ggsave(filename = paste("00_Dist_rawcounts_samples_grid_", project, ".pdf", sep=""), height = 8, width = 6, plot = last_plot(), path = dir_fig_raw)


## MODELING COUNT DATA
#
# Verify data is suited to used the Negtive Binomial model instead of Poisson distribution
# To fit the NB model the standard deviation > mean
# To fit the Poisson the standard deviation == mean

model_m <- data.frame(Mean = apply(raw_counts, 1, mean), Variance = apply(raw_counts, 1, var))
head(model_m)

ggplot(data = model_m, aes(x = Mean, y = Variance))+
  geom_point(alpha = 0.20, colour = "#696969")+
  geom_line(aes(x = Mean, y = Mean), color = "red", show.legend = FALSE)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Mean gene expression level (log10 scale)", y = " Gene-level variance (log10 scale)", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste("00_NB_fit_data_rawcounts_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)


################################################################################
#                              GENE COUNTS
################################################################################


## Statistical summary

print("Dimension raw gene count matrix")
print(dim(raw_counts))

print("Dimension gene count matrix")
print(dim(gene_counts))

print("Summary gene count matrix")
print(summary(gene_counts))

print("There is NA data?")
print(is.na(gene_counts) %>% table())


# Stack matrix
st <- stack(gene_counts)

# Log2 transformation
colnames(dds_log2) <- sample_info$Sample
st_rg <- as.data.frame(stack(dds_log2))
st_rg[trt] <- rep(lvl_ord, each = dim(gene_counts)[1]*length(lvl_ord))
st_rg$col <- as.character(st_rg$col)

# Barplot data
bp <- data.frame(Samples = colnames(gene_counts),
                 trt = sample_info[trt],
                 Values = as.vector(colSums(gene_counts)))


# Log 2 transform data for plotting it
print("Regularized log 2 transformation")
m <- dds_rlog
colnames(m) <- sample_info$Sample
print(head(m))


## COUNTS DISTRIBUTION PLOTS
# Distribution of gene read counts, after filtering, per each sample
ggplot(data = st, aes(x = values, fill = ind))+
  geom_histogram(stat = "bin", bins = 200, position="identity", alpha = 0.3, show.legend = FALSE)+
  xlim(-10, 2000)+
  labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")
ggsave(filename = paste("01_Dist_genecounts_samples_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)

# Distribution of gene read counts, after filtering, with grid per each sample
ggplot(data = st, aes(x = values, fill = ind))+
  geom_histogram(stat = "bin", bins = 200, show.legend = FALSE)+
  facet_wrap(~ ind, ncol = 4)+
  xlim(-10, 2000)+
  labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")+
  theme(strip.background = element_blank())
ggsave(filename = paste("01_Dist_genecounts_samples_grid_", project, ".pdf", sep=""), height = 8, width = 6, plot = last_plot(), path = dir_fig_raw)


## MODELING COUNT DATA
#
# Verify data is suited to used the Negtive Binomial model instead of Poisson distribution
# To fit the NB model the standard deviation > mean
# To fit the Poisson the standard deviation == mean

model_m <- data.frame(Mean = apply(gene_counts, 1, mean), Variance = apply(gene_counts, 1, var))
head(model_m)

ggplot(data = model_m, aes(x = Mean, y = Variance))+
  geom_point(alpha = 0.20, colour = "#696969")+
  geom_line(aes(x = Mean, y = Mean), color = "red", show.legend = FALSE)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Mean gene expression level (log10 scale)", y = " Gene-level variance (log10 scale)", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste("01_NB_fit_data_genecounts_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)


## BARPLOT
# Total gene counts per sample
ggplot(bp, aes(x = factor(Samples, levels = sample_info$Sample), y = Values, fill = !!as.name(trt)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Values, fontface = "bold", vjust = -0.5), size = 1.5,
            position = position_dodge(width = 1), inherit.aes = TRUE)+
  labs(y = "Counts", x = "")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1, size = 5),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste("02_Distribution_bar_genecounts_", project, ".pdf", sep =""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)

ggplot(bp, aes(x = factor(Samples, levels = sample_info$Sample), y = Values, fill = !!as.name(trt)))+
  geom_bar(stat = "identity")+
  labs(y = "Counts", x = "")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1, size = 5),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste("02_Distribution_bar_genecounts_wo_", project, ".pdf", sep =""), height = 4, width = 4, plot = last_plot(), path = dir_fig_raw)


## BOXPLOT
# log 2 transform gene count distribution to visualize possible outliers
ggplot(data = st_rg, aes(x = factor(col, levels = sample_info$Sample), y = value, fill = !!as.name(trt)))+
  geom_boxplot(linewidth = 0.09)+
  labs(x = "", y = "Log2(Counts)")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), panel.grid = element_blank())
ggsave(filename = paste("03_Distribution_boxp_genecounts_",project,".pdf", sep = ""), path = dir_fig_raw, height = 4, width = 4, plot = last_plot())


## DENSITY PLOT
#
# log 2 transform gene count
#
# A few things to note about this figure is that there is a huge peak at exactly
# -3.321928 which can be ignored because this is value is calculated from log2(0.1)
# explained in the box above. This peak consists of all 0-values (inactive genes)
# which isn’t very important to us now.

pdf(paste(dir_fig_raw, "/04_Density_genecounts_", project,".pdf", sep = ""), height = 5, width = 5)
plotDensity(dds_log2, col = rep(color_list[[trt]], each = 4),
            xlab = "Log2(Counts)", ylab = "Density", main = "Expression Distribution")
dev.off()


## PCA PLOTS
#
# To perform the PCA analysis we need a matrix with the PSI values per
# each sample. The samples must be place in the rows and the events in the
# columns. This matrix is centered and scale by using prcomp().
#
# A data frame with the original data and the metadata must be used to plot.
# This data is used to increase interpretation power.
#
# The principal components which explains the variability of the data should
# at least sum 65% of the variance of the data.

# PCA matrix
m_t <- t(m)
print(dim(m_t))
m_pca <- prcomp(m_t, scale. = TRUE)

# Save output
pca1 <- m_pca$x
pca2 <- m_pca$rotation

# Scree plot
# Percentage of variances explained by each principal component
pca_scree <- fviz_eig(m_pca,choice = "variance", addlabels = TRUE, ggtheme = theme_classic())
ggsave(filename = paste("PCA_Screeplot_", project, ".pdf", sep =""),
       plot = pca_scree, path = dir_fig_raw, height = 4, width = 6)

# Variable graphs
# Contributions of the events in each principal component
# Positive correlated variables point to the same side of the plot
# Negative correlated variables point to opposite sides of the graph.
pca_var <- fviz_pca_var(m_pca, col.var = "contrib", gradient.cols = c("blue", "yellow", "red"))+
  theme_classic()
ggsave(filename = paste("PCA_VarContribution_", project, ".pdf", sep =""),
       plot = pca_var, path = dir_fig_raw, height = 4, width = 6)

# PC1 vs PC2
pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1,2),
                         geom.ind = "text", repel = TRUE, labelsize = 4,
                         col.ind = sample_info[[trt]],
                         addEllipses = TRUE, ellipse.level = 0.95,
                         title = "")+
  scale_color_manual(values = color_list[[trt]])+
  scale_fill_manual(values = color_list[[trt]])+
  theme(legend.position = "none")

# PC1 vs PC3
pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1,3),
                         geom.ind = "text", repel = TRUE, labelsize = 4,
                         col.ind = sample_info[[trt]],
                         addEllipses = TRUE, ellipse.level = 0.95,
                         legend.title = "Treatment", title = "")+
  scale_color_manual(values = color_list[[trt]])+
  scale_fill_manual(values = color_list[[trt]])+
  theme(legend.position = "none")

# PC1 vs PC4
pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1,4),
                         geom.ind = "text", repel = TRUE, labelsize = 4,
                         col.ind = sample_info[[trt]],
                         addEllipses = TRUE, ellipse.level = 0.95,
                         legend.title = "Treatment", title = "")+
  scale_color_manual(values = color_list[[trt]])+
  scale_fill_manual(values = color_list[[trt]])+
  theme(legend.position = "none")

# Save PCA plots
ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_gene_counts_", project, ".pdf", sep =""),
       plot = pca_1vs2, path = dir_fig_raw,
       height = 5, width = 6)
ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_gene_counts_", project, ".pdf", sep =""),
       plot = pca_1vs3, path = dir_fig_raw,
       height = 5, width = 6)
ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_gene_counts_", project, ".pdf", sep =""),
       plot = pca_1vs4, path = dir_fig_raw,
       height = 5, width = 6)



## DISTANCE MATRIX
#
# Execute a Euclidean distance to performed sample-to-sample analysis
#
# The results of the distance matrix must be aligned with the results
# in the heatmap and PCA.

eum <- as.matrix(dist(t(m)))

pheatmap(eum,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         border_color = NA,
         treeheight_row = 0,
         filename = paste(dir_fig_raw, "/Cor_euclidean_genecounts_", project, ".pdf", sep = ""),
         height = 4, width = 4)


##  CORRELATION 
#
# Execute a Pearson correlation which accept possible NAs.
# The acceptance of the NA should be valuable in the future in case, we accept
# comparisons when a sample is missing a PSI value.
#
# The results of the correlation matrix must be aligned with the results
# in the heatmap and PCA.

# Pearson Correlation Matrix
pem <- cor(m, method = "pearson", use = "na.or.complete")

# Plot
pheatmap(pem,
         color = colorRampPalette(brewer.pal(9, "Blues"))(255),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 0,
         filename = paste(dir_fig_raw, "/Cor_pearson_genecounts_", project, ".pdf", sep = ""),
         height = 4, width = 4)


## HEATMAP
#
# Clustering of all the differentially expressed genes with different layers
# of information such as the treatment and the event type.
# The heatmap function needs a
#   - Matrix data (m): Samples as columns and genes as rows
#   - Annotation samples (sample_col): Treatment corresponding to each sample
#   - Annotation genes (sample_row): Event type corresponding to each genes

# Column information
# Sample columns information
# Should only contain the samples of the tissue followed by the treatment
sample_col <-  sample_info %>% dplyr::select(all_of(trt))
rownames(sample_col) <- sample_info$Sample
sample_col[trt] <- factor(sample_info[[trt]])

pheatmap(m,
         scale = "row",
         color = color_list$Heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sample_col,
         annotation_colors = color_list,
         show_rownames = FALSE,
         border_color = NA,
         treeheight_row = 0,
         filename = paste(dir_fig_raw, "/Heatmap_Zscore_genecounts_", project, ".pdf", sep = "")
)

pheatmap(m,
         scale = "row",
         color = color_list$Heatmap,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = sample_col,
         annotation_colors = color_list,
         show_rownames = FALSE,
         border_color = NA,
         treeheight_row = 0,
         filename = paste(dir_fig_raw, "/Heatmap_Zscore_genecounts_force_", project, ".pdf", sep = "")
)


## Library size estimation factor
# Save a plot with the library estimation size performed by DESeq2 and EdgR/limma-voom
# DESeq2 use the mean of the ratios method
# EdgeR use the trimmed mean of M-values  method

lib_size <- as.data.frame(cbind(dds$sizeFactor, deg$samples$norm.factors))
colnames(lib_size) <- c("DESeq2", "EgdeR/limma")
norm_fact <- stack(lib_size)
norm_fact$Sample <- sample_info$Sample

ggplot(norm_fact, aes(x = Sample, y = values, color = ind))+
  geom_point(position = position_dodge(0.3), size = 2)+
  geom_line(mapping = aes(group = ind))+
  labs(y = "Factors estimation")+
  scale_color_manual(values = c("#87CEEB","#228B22"))+
  theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank(), legend.position = "top")
ggsave(filename = paste("Size_factors_", project, ".pdf", sep = ""), height = 5, width = 6, plot = last_plot(), path = dir_fig_raw)


#################################################################################################################################
#                             COMPARISONS
#################################################################################################################################


for(i in 1:length(contrast)){
  
  # Contrast 
  name <- paste(contrast[[i]][2], "vs", contrast[[i]][3], sep = "")
  print(name)
  
  # Comparison levels
  comp_lvl <- c(contrast[[i]][2],  contrast[[i]][3])
  
  # Select the contrast levels
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% comp_lvl)]
  
  # Metadata
  metadata <- sample_info[which(sample_info[[trt]] %in% comp_lvl),]
  metadata[,trt] <- factor(metadata[,trt])
  
  # Save the results per each comparison
  res_all <- data.frame()
  
  # Save DEGs per each comparison
  deg_data <- data.frame()
  
  
  ###############################################################
  #                Create comparison folder
  ###############################################################
  
  
  # Create a folder with the comparison name. This folder contains 
  # several folders classified in Results and Figures.
  dir.create(file.path(dir_out, name), showWarnings = FALSE)
  dir_outfolder <- paste(dir_out, "/", name, sep = "")
  
  # <- Final results ->
  # Files folder
  dir.create(file.path(dir_outfolder ,"Results"), showWarnings = FALSE)
  dir_outfiles <- paste(dir_outfolder ,"/Results",sep='')
  # Figures folder
  dir.create(file.path(dir_outfolder , "Shared"), showWarnings = FALSE)
  dir_sh <- paste(dir_outfolder ,"/Shared", sep='')
  
  
  loop <- ifelse(length(metadata$Sample) >= 30, 4, 3)
  
  for(j in 1:loop){
    if (j == 1){
      
      ##########################################################################
      #                             Run DESeq2
      ##########################################################################
      
      
      analysis <- "DESeq2"
      print(analysis)
      
      # Figures folder
      dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
      dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
      
      
      # Log2 fold change result table for an specific comparison
      resl <- results(object = dds,                      # DESeqDataSet  
                      contrast = contrast[[i]],          # Constrast 
                      test = "Wald",                     # default; Wald Test
                      minmu = 0.5,                       # default; lower bound on the estimate count while fitting the GLM
                      
                      ## Independent filtering  
                      alpha = 0.1,              # default = 0.1
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
      
      
      if(is.null(shk) != TRUE){
        shk_res <- lfcShrink(dds = dds, res = resl, contrast = contrast[[i]], type = shk, quiet = TRUE)
        
        pdf(paste(dir_fig, "/01_MA_plot_lfcshk_", analysis, "_", project, "_", name ,".pdf", sep = ""), height = 4, width = 5)
        DESeq2::plotMA(shk_res)
        dev.off()
        
        # Change columns names to plot data 
        colnames(shk_res) <- paste("shrk", shk, colnames(shk_res), sep = "_")
        
        shk_res <- shk_res[match(rownames(gene_counts), rownames(shk_res)),]
        res_all <- cbind(res_all, shk_res)
        
      }else{print("Shrinkage estimators won't be estimated")}
      
      
    } else if (j == 2) {
      
      
      ##########################################################################
      #                     Run EdgeR
      ##########################################################################
      
      
      analysis <- "EdgeR"
      print(analysis)
      
      # Figures folder
      dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
      dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
      
      
      control <- contrast[[i]][3]
      experimental <- contrast[[i]][2]
      
      # Perform contrast
      fit <- exactTest(degde, pair = c(control, experimental), dispersion = "auto", rejection.region = "doubletail")
      
      # Estimate the adjusted p-value 
      res <- as.data.frame(topTags(fit, adjust.method = "BH", n = dim(gene_counts)[1]))
      res <- res[match(rownames(gene_counts), rownames(res)),]
      colnames(res) <- paste(analysis, colnames(res), sep = "_")
      res_all <- cbind(res_all, as.data.frame(res))
      
      # Change columns names to plot data 
      colnames(res) <- c("logFC", "logCPM", "pvalue", "padj")
      
      labels_sig <- rownames(res[which(res$padj <= fdr_cutoff & abs(res$logFC) > lfc_cutoff),])
      
      # MA plot 
      # In the x axis, average logCPM(counts per million) but in DESeq2 is mean of normalized counts
      pdf(file = paste(dir_fig, "/00_MA_plot_", analysis, "_", project,".pdf", sep = ""), height = 5 , width = 5)  
      plotSmear(fit, de.tags = labels_sig, ylab = "log fold change"); abline(h = 0, col  ="black", lty = 1, lwd = 2)
      dev.off()
      
      
    } else if (j == 3) {
      
      
      ##########################################################################
      #                     Run limma-voom
      ##########################################################################
      
      
      analysis <- "limma-voom"
      print(analysis)
      
      # Figures folder
      dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
      dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
      
      # Create contrast
      control <- paste(trt, contrast[[i]][3], sep = "")
      experimental <- paste(trt, contrast[[i]][2], sep = "")
      cont <- makeContrasts(paste(experimental, "-", control, sep = ""), levels = colnames(coef(dlim)))
      
      # Estimate contrast per each gene
      tmp <- contrasts.fit(dlim, cont)
      tmp <- eBayes(tmp)
      
      # Extract results
      res <- topTable(tmp, adjust.method = correction, n = dim(gene_counts)[1])
      res <- res[match(rownames(gene_counts), rownames(res)),]
      colnames(res) <- paste(analysis, colnames(res), sep = "_")
      res_all <- cbind(res_all, as.data.frame(res))
      
      # Change columns names to plot data 
      colnames(res) <- c("logFC", "MeanExp", "t", "pvalue", "padj", "B")
      
      # MA plot 
      plot_dt <- res %>% mutate(Color = "No")
      plot_dt$Color[which(res$padj <= fdr_cutoff & abs(res$logFC) > lfc_cutoff)] <- "Color"
      
      ggplot(plot_dt, aes(x = MeanExp, y = logFC, color = Color))+
        geom_point()+
        geom_hline(yintercept = 0, col = "darkgrey", linetype = "solid", size = 2)+
        labs(x = "log2 CPM", y = "log fold change")+
        scale_color_manual(values = c("blue", "black"),
                           labels = c("No", "Color"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
      ggsave(filename = paste(dir_fig, "/00_MA_plot_", analysis, "_", project, "_", name, ".pdf", sep = ""), height = 4, width = 6, plot = last_plot())
      
    } else {
      
      
      ##########################################################################
      #                   Run Wilcoxon rank-sum test
      ##########################################################################
      
      
      analysis <- "Wilcoxon"
      print(analysis)
      
      # Figures folder
      dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
      dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
      
      
      # Perform TMM normalization and transform into conts per million (CPM). This step
      # it is recommended by researches because Wilcoxon is not a regression based 
      # model and thus it cannot adjust for possible cofunding factors
      count_cpm <- as.data.frame(cpm(deg))
      count_cpm <- count_cpm[, metadata$Sample]
      
      # Perform Wilcoxon test
      res <- wilcoxon_test(count_cpm, correction, metadata, trt, contrast)
      
      # Extract results
      resl <- res
      colnames(resl) <- paste(analysis, colnames(res), sep = "_")
      res_all <- cbind(res_all, resl)
      
      # Completar y hacer una funcion con esto 
      # pseudo_ref <- sqrt(rowSums(count_cpm))
      # ratios <- count_cpm/pseudo_ref
      # median <- apply(ratios, 2, function(x) median(x[which(!is.na(x))]))
      
    }
    
    ############################################################################
    #                       DATA PROCESSING
    ############################################################################
    
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
    
    # List of DEGs per each analysis method used
    t <- df[, c("GeneID", "logFC", "padj", "Direction")] %>% mutate(Group = analysis)
    deg_data <- rbind(deg_data, t)
    
    
    # Data frame per comparison an analysis
    if(analysis == "DESeq2" & !is.null(shk)){
      data_df <- cbind(res_df, shk_res); data_df <- cbind(data_df, gene_counts)
    }else {data_df <- cbind(res_df, gene_counts)}
    data_df <- data_df %>% select(GeneID, Symbol, DEG, Direction, everything())
    
    # Save data per each comparison 
    write.table(data_df, file = paste(dir_outfiles, "/", name, "_", project , "_", analysis, "_", label, "_", threshold, "_info_table.txt" , sep =""))
    
    
    ###############################################################
    #                 Statistical summary
    ###############################################################
    
    
    a <- data.frame(Comparison = name, Method = analysis, Genes = length(!is.na(res_df$padj)), DEGs = length(which(res_df$DEG == "YES")), Upregulated = length(which(res_df$Direction == "Upregulated")), Downregulated = length(which(res_df$Direction == "Downregulated")))
    out_df <- rbind(out_df, a)
    
    
    ############################################################################
    #                       VISUALIZATION
    ############################################################################
    
    
    ### All genes
    # Histogram FC distribution
    A <- ggplot(data = res_df, aes(x = logFC)) +
      geom_histogram( fill = "#6696CC", color = "black")+
      labs(x = "log2FC", y = "Counts")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # Histogram p-value adjusted distribution
    B <- ggplot(data = res_df, aes(x = padj)) +
      geom_histogram( fill = "#6696CC", color = "black")+
      labs(x = "Adjusted p-value", y = "Counts")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    ### Differentially expressed genes
    # Histogram FC distribution 
    C <- ggplot(data = df, aes(x = logFC)) +
      geom_histogram( fill = "#6696CC", color = "black")+
      labs(x = "log2FC", y = "Counts")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # Histogram p-value adjusted distribution
    D <- ggplot(data = df, aes(x = padj)) +
      geom_histogram( fill = "#6696CC", color = "black")+
      labs(x = "Adjusted p-value", y = "Counts")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    fig <- ggarrange(A, B, C, D, ncol = 2, nrow = 2, widths = 10, heights = 10)
    ggsave(filename = paste("00_Validation_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), width = 10, height = 10, plot = fig, path = dir_fig)
    
    
    
    ######################### VOLCANO PLOT #########################
    
    # Possibility of adding labels from interesting gene ids
    cont_res <- res_df
    # cont_res$labels <- ifelse(cont_res$Symbol %in% head(cont_res[order(abs(cont_res$logFC)), "Symbol"], 10), cont_res$Symbol, NA)
    
    # DEGs
    # ggplot(data = cont_res, aes(x = logFC, y = -log10(padj), col = DEG, label = labels))+
    ggplot(data = cont_res, aes(x = logFC, y = -log10(padj), col = DEG))+
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
      geom_point(size = 2, alpha = 0.6)+
      scale_color_manual(values = c("grey", "blue"),
                         labels = c("Not significant", "Significative"))+
      labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)")+
      # geom_text_repel(size=2.5, max.overlaps = Inf)+
      theme(text = element_text(size = 8),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none")
    ggsave(paste("Volcano_DEGs_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 5, width = 5, plot = last_plot(), path = dir_fig)
    
    ggplot(data = res, mapping = aes(x = pvalue))+ 
      geom_histogram(aes(y = ..density..), colour = "black", fill = "white")+
      geom_density(alpha = .2, fill = "blue", linewidth = 0.7)
    
    # DEGs with up and down information
    # ggplot(data = cont_res, aes(x = logFC, y = -log10(padj), col = Direction, label = labels))+
    ggplot(data = cont_res, aes(x = logFC, y = -log10(padj), col = Direction))+
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      geom_hline(yintercept = -log10(fdr_cutoff), col = "gray", linetype = 'dashed')+
      geom_point(size = 2, alpha = 0.6)+
      scale_color_manual(values = as.vector(color_list$Direction))+
      labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "")+
      # geom_text_repel(size=2.5, max.overlaps = Inf)+
      theme(text = element_text(size = 8),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none")
    ggsave(paste("Volcano_DEGs_Direction_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 5, width = 5, plot = last_plot(), path = dir_fig)
    
    
    ######################### WATERFALL PLOT #########################
    
    dt_wf <- df[which(!(duplicated(df$GeneID))),]
    dt_wf$GeneID <- factor(dt_wf$GeneID, levels = dt_wf$GeneID[order(dt_wf$logFC, decreasing = FALSE)])
    
    ggplot(dt_wf, aes(x = GeneID, y = logFC, fill = Direction)) +
      geom_bar(stat = "identity")+
      scale_fill_manual(values = as.vector(color_l$Direction)[-2])+
      xlab("Differentially expressed genes")+
      ylab("Log2 Fold Change")+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
    ggsave(paste("Waterfall_DEGs_",analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), height = 5, width = 5, plot = last_plot(), path = dir_fig)
    
    
    ## PCA PLOTS
    #
    # To perform the PCA analysis we need a matrix with the PSI values per
    # each sample. The samples must be place in the rows and the events in the
    # columns. This matrix is centered and scale by using prcomp().
    #
    # A data frame with the original data and the metadata must be used to plot.
    # This data is used to increase interpretation power.
    #
    # The principal components which explains the variability of the data should
    # at least sum 65% of the variance of the data.
    
    # PCA matrix
    m_t <- t(m)
    print(dim(m_t))
    m_pca <- prcomp(m_t, scale. = TRUE)
    
    # Save output
    pca1 <- m_pca$x
    pca2 <- m_pca$rotation
    
    # Scree plot
    # Percentage of variances explained by each principal component
    pca_scree <- fviz_eig(m_pca,choice = "variance", addlabels = TRUE, ggtheme = theme_classic(), main = "Screeplot")
    
    # Variable graphs
    # Contributions of the events in each principal component
    # Positive correlated variables point to the same side of the plot
    # Negative correlated variables point to opposite sides of the graph.
    pca_var <- fviz_pca_var(m_pca, col.var = "contrib", gradient.cols = c("blue", "yellow", "red"))+
      theme_classic()+ 
      labs(title = "Variance Contribution")
    
    fig <- ggarrange(pca_scree, pca_var, ncol = 2, nrow = 1, widths = 10, heights = 4)
    ggsave(filename = paste("PCA_params_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""), width = 10, height = 4, plot = fig, path = dir_fig)
    
    
    # PC1 vs PC2
    pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1,2),
                             geom.ind = "text", repel = TRUE, labelsize = 4,
                             col.ind = metadata[[trt]],
                             addEllipses = TRUE, ellipse.level = 0.95,
                             title = "")+
      scale_color_manual(values = color_l[[trt]])+
      scale_fill_manual(values = color_l[[trt]])+
      theme(legend.position = "none")
    
    # PC1 vs PC3
    pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1,3),
                             geom.ind = "text", repel = TRUE, labelsize = 4,
                             col.ind = metadata[[trt]],
                             addEllipses = TRUE, ellipse.level = 0.95,
                             legend.title = "Treatment", title = "")+
      scale_color_manual(values = color_l[[trt]])+
      scale_fill_manual(values = color_l[[trt]])+
      theme(legend.position = "none")
    
    # PC1 vs PC4
    pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1,4),
                             geom.ind = "text", repel = TRUE, labelsize = 4,
                             col.ind = metadata[[trt]],
                             addEllipses = TRUE, ellipse.level = 0.95,
                             legend.title = "Treatment", title = "")+
      scale_color_manual(values = color_l[[trt]])+
      scale_fill_manual(values = color_l[[trt]])+
      theme(legend.position = "none")
    
    # Save PCA plots
    ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""),
           plot = pca_1vs2, path = dir_fig,
           height = 5, width = 6)
    ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""),
           plot = pca_1vs3, path = dir_fig,
           height = 5, width = 6)
    ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep =""),
           plot = pca_1vs4, path = dir_fig,
           height = 5, width = 6)
    
    
    
    ## DISTANCE MATRIX
    #
    # Execute a Euclidean distance to performed sample-to-sample analysis
    #
    # The results of the distance matrix must be aligned with the results
    # in the heatmap and PCA.
    
    eum <- as.matrix(dist(t(m)))
    
    pheatmap(eum,
             color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 6,
             fontsize_col = 6,
             border_color = NA,
             treeheight_row = 0,
             filename = paste(dir_fig, "/Euclidean_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""),
             height = 4, width = 4)
    
    
    ##  CORRELATION 
    #
    # Execute a Pearson correlation which accept possible NAs.
    # The acceptance of the NA should be valuable in the future in case, we accept
    # comparisons when a sample is missing a PSI value.
    #
    # The results of the correlation matrix must be aligned with the results
    # in the heatmap and PCA.
    
    # Pearson Correlation Matrix
    pem <- cor(m, method = "pearson", use = "na.or.complete")
    
    # Plot
    pheatmap(pem,
             color = colorRampPalette(brewer.pal(9, "Blues"))(255),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 6,
             fontsize_col = 6,
             border_color = NA,
             treeheight_row = 0,
             treeheight_col = 0,
             filename = paste(dir_fig, "/Pearson_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = ""),
             height = 4, width = 4)
    
    
    ## HEATMAP
    #
    # Clustering of all the differentially expressed genes with different layers
    # of information such as the treatment and the event type.
    # The heatmap function needs a
    #   - Matrix data (m): Samples as columns and genes as rows
    #   - Annotation samples (sample_col): Treatment corresponding to each sample
    #   - Annotation genes (sample_row): Event type corresponding to each genes
    
    # Column information
    # Sample columns information
    # Should only contain the samples of the tissue followed by the treatment
    sample_col <-  metadata %>% dplyr::select(all_of(trt))
    rownames(sample_col) <- metadata$Sample
    sample_col[trt] <- factor(metadata[[trt]])
    
    # # Row information
    # # Should only contain the differential events of the tissue
    # # Info: Name (Gen_Event) and event type
    # sample_row <- df_sel %>%
    #   distinct(NAME, .keep_all = TRUE) %>%
    #   select(NAME, which(colnames(df_sel) %in% paste("PSI_", metadata$Sample, sep = "")), EVENT_TYPE)
    # sample_row <- data.frame(COMPLEX = sample_row$EVENT_TYPE)
    # rownames(sample_row) <- df_sel$NAME
    
    
    # Choosing heatmap font size
    if (dim(m)[1]>=60 & dim(m)[1]<200){font_row = 2} else if (dim(m)[1]>=200){font_row = 1} else if (dim(m)[1]<=60 & dim(m)[1]>=20){font_row = 2} else {font_row = 5}
    if (dim(m)[2] <= 10){font_col = 7} else if (dim(m)[2] > 10 & dim(m)[2] < 20){font_col = 5} else {font_col = 2}
    
    
    pheatmap(m,
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
             clustering_callback = callback,   # Sorting dendrogram
             filename = paste(dir_fig, "/Heatmap_Zscore_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = "")
    )
    
    pheatmap(m,
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
             clustering_callback = callback,   # Sorting dendrogram
             filename = paste(dir_fig, "/Heatmap_Zscore_nonames_", analysis, "_", threshold, "_", name, "_", project, ".pdf", sep = "")
    )
    
    
  }
  
  
  
  #########################################################################################################################
  #                       SHARED GENES BETWEEN METHODS
  #########################################################################################################################
  
  
  ##############################################################################
  #                         SHARED EVENTS
  ##############################################################################
  
  
  # List of shared events
  deg_list <- split(deg_data$GeneID, deg_data$Group)
  
  # List of the duplicated events
  dup_tab <- as.data.frame(table(deg_data$GeneID))
  dup_gen <- dup_tab[dup_tab$Freq > 1,]
  
  
  # Duplicate events per tissue
  #   - Times repeated per group
  #   - Groups in which appear
  #   - Direction of the event (up/down)
  #   - Congruence: Direction the event was the same or not in 
  #     both groups
  duplicates <- deg_data[which(deg_data$GeneID %in% dup_gen$Var1),] %>% 
    group_by(GeneID) %>%
    summarise(N = n_distinct(Group),
              Group = paste(unique(Group), collapse = ", "),
              Direction = paste(Direction, collapse = ", ")) %>%
    filter(N > 1)
  
  # Congruence add to the duplicates
  duplicates$Congruence <- sapply(duplicates$Direction, check_congruence)
  
  
  ########### VENN DIAGRAM ###########
  
  # Selecting colors
  if(length(names(deg_list))==3) {ven_col <- color_list$Shared[-4]} else {ven_col <- color_list$Shared}
  
  # Visualization of the shared genes among the different methods used
  ggvenn(deg_list,
         digits = 2,                                    # Two decimals
         fill_color = ven_col,       # Fill color
         fill_alpha = 0.5,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 3,                           # Default 
         text_size = 2.5                              # Default
  )
  ggsave(filename = paste("00_Venn_Analysis_", threshold, "_", name, "_", project, ".pdf", sep =""), path = dir_sh, height = 4, width = 4, plot = last_plot())
  
  
  ########### UPSET PLOT ###########
  # 
  # Visualization of the shared genes among the different methods used
  upset1 <- upset(fromList(deg_list),                       
                  order.by = "freq",                          # Order by frequency
                  nsets = length(deg_list),                  # Number of sets
                  keep.order = TRUE,                          
                  mainbar.y.label = "Shared genes",       # Y axis label
                  sets.x.label = "Total genes",                    # X axis label
                  main.bar.color = "#2E8B57",      # Y axis bar color
                  sets.bar.color = "#483D8B",      # X axis bar color
                  matrix.color =  "black"        # Matrix dots color
  )     
  
  pdf(file = paste(dir_sh, "/00_Upset_Analysis_", threshold, "_", name, "_", project,".pdf", sep = ""), height = 5, width = 10, bg = "white")
  print(upset1)
  dev.off()
  
  
  ##############################################################################
  #                           UNIQUE EVENTS
  ##############################################################################
  
  
  # List of unique events
  unq_tab <- dup_tab[dup_tab$Freq == 1,]
  
  # Data frame of unique events
  unq_events <- deg_data[which(deg_data$GeneID %in% unq_tab$Var1),]
  
  
  ##############################################################################
  #                           UNIQUE EVENTS
  ##############################################################################
  
  
  s_all <- stack(res_all)
  s_all <- separate(s_all, col = ind, into = c("Group", "Parameter"), sep = "_")
  
  s <- s_all[s_all$Parameter %in% c("log2FoldChange", "logFC", "pvalue", "padj", "PValue", "FDR", "P.Value", "adj.P.Val"),]
  s$Parameter[which(s$Parameter == "log2FoldChange")] = "logFC"
  s$Parameter[which(s$Parameter %in% c("PValue", "P.Value"))] = "pvalue"
  s$Parameter[which(s$Parameter %in% c("adj.P.Val", "FDR"))] = "padj"
  s$GeneID <- rownames(res_all)
  
  
  
  for (l in 2:length(names(deg_list))){
    a <- "DESeq2"
    b <- names(deg_list)[l]
    
    subset_dt <- s[which(s$Group %in% c(a,b)),]
    sub_dt <- subset_dt %>% mutate(Param = paste(Group, Parameter, sep = "_")) %>% select(-c(Group, Parameter, GeneID))
    sub_dt <- unstack(sub_dt)
    colnames(sub_dt) <- c("A", "B", "C", "A2", "B2", "C2")
    sub_dt$GeneID <- unique(subset_dt$GeneID)
    head(sub_dt)
    
    
    A <- ggplot(sub_dt, aes(x = -log10(B), y = -log10(B2), color = GeneID))+
      geom_point(size = 1, alpha = 0.5, show.legend = FALSE) +
      labs(x = paste(a, "pvalue", sep = "_"), y = paste(b, "pvalue", sep = "_"))+
      theme(text = element_text(size = 13,  family = "serif"), axis.text = element_text(size = 12), 
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
    
    B <- ggplot(sub_dt, aes(x = log10(C), y = -log10(C2), color = GeneID))+
      geom_point(size = 1, alpha = 0.5, show.legend = FALSE)+
      labs(x = paste(a, "padj", sep = "_"), y = paste(b, "adj", sep = "_"))+
      theme(text = element_text(size = 13,  family = "serif"), axis.text = element_text(size = 12), 
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))    
    
    C <- ggplot(sub_dt, aes(x = A, y = A2 , color = GeneID))+
      geom_point(size = 1, alpha = 0.5, show.legend = FALSE)+
      labs(x = paste(a, "logFC", sep = "_"), y = paste(b, "logFC", sep = "_"))+
      theme(text = element_text(size = 13,  family = "serif"), axis.text = element_text(size = 12), 
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
    
    
    fig <- ggarrange(A, B, C, ncol = 2, nrow = 2, widths = 10, heights = 10)
    ggsave(filename = paste("01_", a, "vs", b, "_", threshold, "_", name, "_", project, ".pdf", sep =""), width = 10, height = 10, plot = fig, path = dir_sh)
    
    
  }
  
  
  ################################################################################
  #                            SAVE DATA
  ################################################################################
  
  # Add information rows
  dupli <- duplicates %>% select(-c("N", "Direction", "Congruence")) 
  
  res_all$GeneID <- rownames(res_all)
  res_all$Symbol <- gene_names$Symbol
  res_all <- merge(res_all, dupli, by = "GeneID", all = TRUE)
  res_all <- res_all %>% select(GeneID, Symbol, Group, everything())
  
  write.table(res_all, paste(dir_outfiles, "/DEGs_", threshold, "_", name, "_", project, ".txt", sep = ""), row.names = FALSE)
  write.xlsx(res_all, paste(dir_outfiles, "/DEGs_", threshold, "_", name, "_", project, ".xlsx", sep = ""), overwrite = TRUE)
  
  write.csv(duplicates, paste(dir_outfiles, "/Shared_DEGs_", threshold, "_", name, "_", project, ".csv", sep = ""), row.names = FALSE)
  write.csv(unq_events, paste(dir_outfiles, "/Unique_DEGs_", threshold, "_", name, "_", project, ".csv", sep = ""), row.names = FALSE)
  
}


# Save file 
write.table(data, paste(dir_raw, "/Count_matrix_norm_VST_rlog_log2_", label, "_", project, ".txt", sep = ""), row.names = FALSE)

# Save library size factors, normalization methods
write.table(lib_size, paste(dir_raw, "/Library_size_factors_",project,".txt", sep = ""), row.names = TRUE)

# Save the summary per each comparison and method
write.csv(out_df, paste(dir_raw, "/Comparison_summary_", threshold, "_", project, ".csv", sep = ""), row.names = FALSE)

# Save R session information
sessionInfo() %>% capture.output(file = paste("session_info_", format(Sys.time(), "%Y%m%d.%H%M%S"), ".txt", sep = ""))

