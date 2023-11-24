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
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Date of the log file 0_Sample_info_XXXX.log
logdate <- "20231121"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/log/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory. Raw gene counts  
dir_infiles <- paste(path, project, "/04_STAR/RawCounts_", project,".txt", sep = "")

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
# filtered. The filtering method is based on the function filterByExpr from edgeR.
# As a result. a unique matrix per each comparison is obtain and given to DESeq2, 
# limma-voom, EdgeR and Wilcoxon. However, there is a primary filterin step in 
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
# log2 fold change and not fold change <------
lfc_cutoff <- log2(1.5)         # Default option

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- "BH"


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

# Automatically generate the colors for the treatment condition
if(length(color_list)<4){color_list <- color_palette(color_list, trt, lvl_ord, palette = "Dark2")}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  

# Genome annotation path
anot_path <- sub(pattern = "/.*", replacement = "", path)
if (anot_path == ""){anot_path <- "/vols/GPArkaitz_bigdata"}

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
ref_genome <- read.table(paste(anot_path, "/DATA_shared/Genomes_Rocky/DEG_Annotation/Annotated_Genes_20231121.txt", sep = ""), header = TRUE)
genome <- ref_genome[which(ref_genome$Specie == specie), -3]
genome$Name <- paste(genome$Symbol, genome$Ensembl, sep = "_")
print(dim(genome))


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
dir.create(file.path(dir_out,"05_DEG_ANALYSIS"))
dir_out <- paste(dir_out,"/05_DEG_ANALYSIS", sep='')
setwd(dir_out)



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
  project <- paste(project, outliers, sep = "_", collapse = "_")
}else{print("No sample was considered an outlier.")}


#### Add covariates effects #### 
if(is.null(var_exp) == FALSE){
  # Normal for EdgeR and limma-voom
  sample_info[, var_exp] <- data_info[, which(colnames(data_info) %in% var_exp)]
  # Z_score for DESeq2
  for (k in 1:length(var_exp)){
    if(length(unique(sample_info[, var_exp[k]]))>1){
      sample_info[, paste(var_exp[k], "_zscore", sep = "")] <- as.vector(scale(data_info[[var_exp[k]]], center = TRUE, scale = TRUE))}}
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
annot <- genome[which(genome$Ensembl %in% rownames(gene_counts)),]
gene_names <- annot[match(rownames(gene_counts), annot$Ensembl), -3]


## Barplot verification the quality of the gene counts
bp <- data.frame(Samples = colnames(gene_counts),
                 trt = sample_info[trt],
                 Values = as.vector(colSums(gene_counts)))

ggplot(bp, aes(x = factor(Samples, levels = sample_info$Sample), y = Values, fill = !!as.name(trt)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Values, fontface = "bold", vjust = -0.5), size = 1.5,
            position = position_dodge(width = 1), inherit.aes = TRUE)+
  labs(y = "Counts", x = "")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1, size = 5), legend.position = "none")
ggsave(filename = paste("Barplot_genecounts_", project, ".pdf", sep =""), height = 4, width = 4, plot = last_plot(), path = dir_out, bg = "white")

ggplot(bp, aes(x = factor(Samples, levels = sample_info$Sample), y = Values, fill = !!as.name(trt)))+
  geom_bar(stat = "identity")+
  labs(y = "Counts", x = "")+
  scale_fill_manual(values = color_list[[trt]])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1, size = 5), legend.position = "none")
ggsave(filename = paste("Barplot_genecounts_wo_", project, ".pdf", sep =""), height = 4, width = 4, plot = last_plot(), path = dir_out, bg = "white")




################################################################################
#                               SAVE DATA 
################################################################################


# Save gene counts
write.table(gene_counts, paste(dir_out,"/GeneCount_filtered_", project,".txt", sep = ""))

# Save metadata information
write.table(sample_info, paste(dir_out, "/Metadata_", project,".txt", sep = ""))

# Save annotation data 
write.table(gene_names, paste(dir_out,"/Annotation_", project,".txt", sep = ""))


################################################################################
#                               COMPARISON 
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
  
  # Metadata
  metadata <- sample_info[which(sample_info[[trt]] %in% comp_lvl),]
  metadata[,trt] <- factor(metadata[,trt])
  
  # Gene counts
  counts_comp <- gene_counts[, metadata$Sample]
  
  # DESeq2 design formula
  # Used to estimate the variance stabilization (VST) method proposed in DESeq2
  # 
  # Steps
  # 1. Check there are covariates
  # 2. Check the values of the covariate are not equal to avoid colinearity in 
  #   the model. If their values are equal the variable is not included in the 
  #   model.
  # 3. Create the design formula
  if(is.null(var_exp) == FALSE){
    var_design <- NULL
    for (k in 1:length(var_exp)){if(length(unique(metadata[, var_exp[k]]))>1){var_design <- c(var_design, var_exp[k])}}
    design_cond = ifelse(is.null(var_design) == FALSE, paste("~", paste(var_design, "_zscore", " +", sep = "", collapse = " "), trt, sep = " "), paste("~", trt, sep = " "))
  } else{design_cond <- paste("~", trt, sep = " ")}
  
  print(design_cond)
  
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  
  # Create a folder for the analysis. This folder contains several folders 
  # classified in Results and Figures. 
  
  # Create output directory
  dir.create(file.path(dir_out, name), showWarnings = FALSE)
  dir_outfolder <- paste(dir_out,"/", name, sep='')
  setwd(dir_outfolder)
  
  # Files folder
  dir.create(file.path(dir_outfolder,"Results"), showWarnings = FALSE)
  dir_output <- paste(dir_outfolder,"/Results", sep='')
  # Quality control figures folder
  dir.create(file.path(dir_outfolder, "Figures"), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder, "/Figures", sep='')
  
  
  
  ##############################################################################
  #                               Filtering 
  ##############################################################################
  
  ## Filter gene count matrix per each comparison
  # 
  # This data will be used when running DESeq2, EdgeR, limma-voom and Wilcoxon to
  # perform the test for the differentially expressed genes
  keep <- filter_genecounts(counts_comp, metadata, trt, min_count = 10, min_prop = 0.7, n_large = 30, min_total = 10)
  df <- counts_comp[keep,]
  
  
  ##############################################################################
  #                           Data transformation 
  ##############################################################################
  
  # Step 1: Create DESeq object
  # ----------------------------------------------------------------------------
  #
  # DESeqDataSetFromMatrix has all the information necessary to run DESeq. 
  # - Count matrix
  # - Metadata with the following columns: condition and covariates. Data 
  #     associated to the samples valuable to perform the contrast. The row names
  #     must be the samples names
  # - Design formula with the experimental conditions and the covariates 
  dds <- DESeqDataSetFromMatrix(countData = df, colData = metadata[,-1], design =  eval(parse(text = design_cond)))
  
  
  # Step 2: Variance stabilizing transformation
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
  
  
  # Estimate the biggest group sample size
  group_n <- max(as.vector(tabulate(metadata[[trt]])))
  
  if(group_n < 30){
    m <- as.data.frame(assay(vst(dds, blind = FALSE)))
    md <- "VST"
  }else{
    m <- as.data.frame(assay(rlog(dds, blind = FALSE)))
    md <- "RLOG"}
  
  # Step 3: Stack matrix
  # ----------------------------------------------------------------------------
  # 
  # Modify the matrix to easily plot it
  m_s <- stack(m)
  m_s_mod <- m_s
  m_s_mod[trt] <- rep(metadata[[trt]], each = dim(df)[1])
  
  
  
  ##############################################################################
  #                           Exploratory Analysis 
  ##############################################################################
  
  
  ## Statistical summary
  
  print("Dimension raw gene count matrix")
  print(dim(counts_comp))
  
  print("Dimension gene count matrix")
  print(dim(df))
  
  print("Summary gene count matrix")
  print(summary(df))
  
  print("There is NA data?")
  print(is.na(df) %>% table())
  
  
   
  ##############################################################################
  #                                 Plots 
  ##############################################################################
  
  
  ## COUNTS DISTRIBUTION PLOTS
  # Distribution of gene read counts per each sample
  ggplot(data = m_s, aes(x = values, fill = ind))+
    geom_histogram(stat = "bin", bins = 200, position="identity", alpha = 0.3, show.legend = FALSE)+
    labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")
  ggsave(filename = paste("00_Dist_genecounts_samples_", name, "_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig, bg = "white")
  
  # Distribution of gene read counts with grid per each sample
  ggplot(data = m_s, aes(x = values, fill = ind))+
    geom_histogram(stat = "bin", bins = 200, show.legend = FALSE)+
    facet_wrap(~ ind, ncol = 4)+
    labs(x = "Counts", y = "Genes", title = "Gene Counts Expression")+
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.01, hjust = 1))
  ggsave(filename = paste("00_Dist_genecounts_samples_grid_", name, "_", project, ".pdf", sep=""), height = 8, width = 6, plot = last_plot(), path = dir_fig, bg = "white")
  
  
  ## BOXPLOT
  ggplot(data = m_s_mod, aes(x = ind, y = values, fill = !!as.name(trt)))+
    geom_boxplot(linewidth = 0.09)+
    labs(x = "", y = "Log2(Counts)")+
    scale_fill_manual(values = color_l[[trt]])+
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), panel.grid = element_blank())
  ggsave(filename = paste("01_Distribution_boxp_genecounts_", name, "_", project,".pdf", sep = ""), path = dir_fig, height = 4, width = 4, plot = last_plot(), bg = "white")

  
  ## DENSITY PLOT
  # A few things to note about this figure is that there is a huge peak at exactly
  # -3.321928 which can be ignored because this is value is calculated from log2(0.1)
  # explained in the box above. This peak consists of all 0-values (inactive genes)
  # which isn’t very important to us now.
  
  pdf(paste(dir_fig, "/02_Density_genecounts_", name, "_", project,".pdf", sep = ""), height = 5, width = 5, bg = "white")
  plotDensity(m, col = rep(color_l[[trt]], each = 4), xlab = "Log2(Counts)", ylab = "Density", main = "Expression Distribution")
  dev.off()
  
  
  ## MODELING COUNT DATA
  # Perform with raw gene gene counts 
  model_m <- data.frame(Mean = apply(df, 1, mean), Variance = apply(df, 1, var))
  head(model_m)
  
  ggplot(data = model_m, aes(x = Mean, y = Variance))+
    geom_point(alpha = 0.20, colour = "#696969")+
    geom_line(aes(x = Mean, y = Mean), color = "red", show.legend = FALSE)+
    scale_y_log10()+
    scale_x_log10()+
    labs(x = "Mean gene expression level (log10 scale)", y = " Gene-level variance (log10 scale)", title = "")
  ggsave(filename = paste("03_NB_fit_data_genecounts_", name, "_", project, ".pdf", sep=""), height = 4, width = 4, plot = last_plot(), path = dir_fig, bg = "white")
  
  
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
  
  pdf(paste(dir_fig, "/05_Cor_pearson_genecounts_", name, "_", project, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_h)
  dev.off()

  
  ## PCA PLOTS
  plot_pcas <- pca_plot(m, trt, metadata, color_l)

  ggsave(filename = paste("PCA_params_genecounts_", name, "_", project, "_scree.pdf", sep = ""), plot = plot_pcas[[1]], path = dir_fig, height = 4, width = 4, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_genecounts_", name, "_", project, ".pdf", sep = ""), plot = plot_pcas[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_genecounts_", name, "_", project, ".pdf", sep = ""), plot = plot_pcas[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_genecounts_", name, "_", project, ".pdf", sep = ""), plot = plot_pcas[[4]], path = dir_fig, height = 5, width = 6, bg = "white")

  ## HEATMAP
  plot_heatmap <- heatmap_plot(m, metadata, trt, color_l)
  
  pdf(paste(dir_fig, "/Heatmap_zscore_genecounts_", name, "_", project, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
  print(plot_heatmap)
  dev.off()
  

  
  ##############################################################################
  #                               Save data 
  ##############################################################################
  
  
  # Save transform data with blind = FALSE
  # Why blind = FALSE
  if(group_n < 30){dds_trs <- assay(vst(dds, blind = FALSE))}else{dds_trs <- assay(rlog(dds, blind = FALSE))}
  write.table(dds_trs, paste(dir_output,"/GeneCount_", md , "_blindFALSE_", name, "_", project, ".txt", sep = ""))
  
  # Save filtered gene counts per comparison 
  write.table(df, paste(dir_output,"/GeneCount_", name, "_", project,".txt", sep = ""))
  
  # Save metadata
  write.table(metadata, paste(dir_output,"/Metadata_", name, "_", project,".txt", sep = ""))
  
  # Save QC file information
  logdate <- format(Sys.time(), "%Y%m%d")
  sum_contrast <- c()
  sum_contrast$Date <- Sys.time()
  sum_contrast$Contrast <- name
  sum_contrast$Genes <- dim(counts_comp)[1]
  sum_contrast$GenesFiltered <- dim(df)[1] 
  sum_contrast$design <- design_cond
  sum_contrast$Transformation <- md
  write.table(sum_contrast, paste(dir_output,"/QC_result_", name, "_", project,".txt", sep = ""))
  
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
log_data$dir_out <- dir_out
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_ord, collapse =",")
log_data$Outliers <- paste(outliers, collapse = ",") 
log_data$Varexp <- paste(var_exp, collapse = ",") 
log_data$min_count <- min_count
log_data$fdr_cutoff <- fdr_cutoff
log_data$lfc_cutoff <- lfc_cutoff
log_data$correction <- correction
log_data$contrast <- paste(unlist(contrast), collapse = ",")
log_data$colortrt <- paste(color_list[[1]], collapse = ",")
log_data$colorheat <- paste(color_list[[2]], collapse = ",")
log_data$colordir<-  paste(color_list[[3]], collapse = ",")
log_data$colorsh <- paste(color_list[[4]], collapse = ",")

write.table(as.data.frame(log_data), paste(path, project, "/log/1_DEG_qc_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



################################################################################
#                             FOLLOW THE ANALYSIS           
################################################################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# source(paste(path, project, "/Scripts/06_DEG_v2_DESeq2.R", sep = ""))
# source(paste(path, project, "/Scripts/06_DEG_v2_EdgeR.R", sep = ""))
# source(paste(path, project, "/Scripts/06_DEG_v2_limma.R", sep = ""))
# source(paste(path, project, "/Script/06_DEG_v2_wilcoxon.R", sep = ""))          #ONLY WITH BIG DATA SETS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


