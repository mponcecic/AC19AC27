################################################################################
#                           PARAMETER SELECTION                       
################################################################################


# RIGTH NOW ALL THE VARIABLES HERE ARE BASED ON THE DEG SCRIPT



# Project name
project <- "AC58"

# Pathway to the folders and files
# Can be your personal folder in BigData
path <- "W:/mponce/"


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

