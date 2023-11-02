################################################################################
#                 DIFFERENTIALLY EXPRESSED GENES SCRIPT
################################################################################

# Folders
# Input:
# Output:





################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Project name
project <- "AC58"

# Pathway to the folders and files
# Can be your personal folder in BigData
path <- "W:/mponce/"

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/03_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
# dir_out <- paste(path, project, sep = "")   # Default option
dir_out <- paste(path, project, sep = "")


# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- "Time"

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- c("Control", "4", "24", "48")


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
# Better to choose the manual 
color_list <- list(trt = c(Control = "#A6DAB0", `4` = "#C18BB7", `24` = "#D7B0B0", `48` = "#8BABD3"), 
                   Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
                   Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
                   Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Data frame with summary information from all the comparisons
# Columns are Comparison, Method, Genes, Upregulated and Downregulated
out_df <- data.frame()



################################################################################
#                               LOAD FILES                       
################################################################################


# Reference genome
# Used in the mapping with STAR which was refdata-gex-GRCh38-2020-A
genome <- read.table(paste(path, project, "/05_DEGs/03_STAR/Human_101/geneInfo.tab", sep = ""), skip = 1, 
                     col.names = c("GeneID", "Symbol", "Biotype"))

# Sample information/ Metadata
data_info <- read.csv(file = paste(path, project, "/Sample_info.csv", sep = ""), header = TRUE)

# Gene count matrix
raw_counts <- read.table(file = dir_infiles, sep = "\t",  header = TRUE, stringsAsFactors = TRUE)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


print(dim(raw_counts))
print(dim(data_info))
print(dim(genome))

print(head(raw_counts))
print(head(data_info))
print(head(genome))



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



#######################################################################
#                            LOG FILE                        
#######################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project_name <- project
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_ord, collapse = ",")
# log_data$path <- dir_out
# log_data$pathRocky <- path_rocky
# log_data$filedir <- files
# log_data$filedirRocky <- files_rocky

# Color list
log_data$ColorTrt <- paste(color_list[[trt]], collapse = ", ")
log_data$ColorHeatmap <- paste(color_list$Heatmap, collapse = ", ")
log_data$ColorDirection <- paste(color_list$Direction, collapse = ", ")
log_data$ColorShared <- paste(color_list$Shared, collapse = ", ")



write.table(as.data.frame(log_data), paste(dir_out, "/DEGs_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")




