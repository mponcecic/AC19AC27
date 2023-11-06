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
logfile <- read.table(paste(path, project, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/03_STAR/RawCounts_", project,".txt", sep = "")

# Output directory
# dir_out <- paste(path, project, sep = "")   # Default option
dir_out <- paste(path, project, sep = "")

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
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))


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

# 
# ## Data transformation considering the metadata information or not, by default 
# # this is set FALSE
# # blind = TRUE, unbiased why to approach data without considering the experimental 
# # design 
# # blind = FALSE, considers the experimental design and differences in the counts 
# # might be attribute to the experimental conditions
# # 
# # Options:
# # blind <- TRUE
# blind <- FALSE        # Default option
# 
# ## Shrinkage estimators
# # Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq run 
# # without LFC shrinkage
# # 
# # Options: 
# # shk <- "ashr"
# # shk <- "apelgm"
# # shk <- "normal" 
# shk <- NULL

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
# specie = "Human"
specie <- logfile$Organism


################################################################################
#                               LOAD FILES                       
################################################################################

# Reference genome
ref_genome <- read.table(paste(path, "/DEG_annotation/gene_annotation_20231106.txt", sep = ""),
                         col.names = c("Symbol", "GeneID", "Organism"))
genome <- ref_genome[which(ref_genome$Organism == specie), -3]

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


################################################################################
#                        SET WORKING DIRECTORY
################################################################################


# Create a folder for the analysis. This folder contains several folders 
# classified in Results and Figures. 

# Create output directory
dir.create(file.path(dir_out,"04_DEG_ANALYSIS"))
dir_out <- paste(dir_out,"/04_DEG_ANALYSIS", sep='')
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
}else{print("Samples in count matrix are ordered")
  gene_counts <- raw_counts[, match(sample_info$Sample, colnames(raw_counts))]
}
  

#### Filter #### 
# Number of samples
n <- ncol(gene_counts)

# Filter genes with no counts based on the number of samples
# This is a common filtering step for all the comparison levels
gene_counts <- gene_counts[which(rowSums(gene_counts) != 0),]
label <- "filtering_0"

print(label)
cat(paste("Total genes in raw count:", dim(raw_counts)[1]))
cat(paste("Total genes in gene count:", dim(gene_counts)[1]))



#### Annotation #### 

# Annotate all the genes in the matrix with the Symbol annotation

# CAUTION: Symbol id present more than one Ensembl identifier. 
# I propose to perform the analysis using the Ensembl identifier, if not the 
# mean value of the genes with the same Symbol id should be performed.
annot <- genome[which(genome$GeneID %in% rownames(gene_counts)),]
gene_names <- data.frame(GeneID = annot$GeneID, Symbol = annot$Symbol)




# Save gene counts
write.table(gene_counts, paste(dir_raw,"/GeneCount_filtered_", project,".txt", sep = ""))
# Save annotation data 
write.table(gene_names, paste(dir_raw,"/Annotation_", project,".txt", sep = ""))




################################################################################
#                           COMPARISON 
################################################################################

#### Create contrast list #### 
if(is.null(contrast)== TRUE){contrast <- create_contrast(trt, lvl_ord)}


for (i in 1:length(contrast)){
  
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
  
  #### Design formula ####
  #
  # Two different design formulas are created for DESeq2 and EdgeR/limma-voom
  # DESeq2 design formula
  design_cond = ifelse(is.null(var_exp) == TRUE, paste("~", trt, sep = " "), paste("~", paste(var_exp, "_zscore", " +", sep = "", collapse = " "), trt, sep = " "))

  # EdgeR and limma-voom
  # Build a design/model matrix
  design_cond2 = ifelse(is.null(var_exp) == TRUE, paste("~ 0 + ", trt, sep = " "), paste("~ 0 +", paste(var_exp, "+", sep = " ", collapse = ""), trt, sep = " "))
  m_model = model.matrix(as.formula(design_cond2), metadata)

  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  # Create a folder for the analysis. This folder contains several folders 
  # classified in Results and Figures. 
  
  # Create output directory
  dir.create(file.path(dir_out, name))
  dir_outfolder <- paste(dir_out,"/", name, sep='')
  setwd(dir_outfolder)
  
  
  # Files folder
  dir.create(file.path(dir_outfolder,"Results"), showWarnings = FALSE)
  dir_output <- paste(dir_outfolder,"/Results", sep='')
  # Quality control figures folder
  dir.create(file.path(dir_outfolder, "Figures"), showWarnings = FALSE)
  dir_fig_qc <- paste(dir_outfolder, "/Figures", sep='')
  
  
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
log_data$dir_out <- dir_raw
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

write.table(as.data.frame(log_data), paste(path, project, "/1_DEG_qc_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")


################################################################################
#                       KEEP THE ANALYSIS GOING            
################################################################################

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(paste(path, project, "/Scripts/06_DEG_v2_DESeq2.R", sep = ""))
source(paste(path, project, "/Scripts/06_DEG_v2_EdgeR.R", sep = ""))
source(paste(path, project, "/Scripts/06_DEG_v2_limma.R", sep = ""))
# source(paste(path, project, "/Script/06_DEG_v2_wilcoxon.R", sep = ""))          #ONLY WITH BIG DATA SETS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


