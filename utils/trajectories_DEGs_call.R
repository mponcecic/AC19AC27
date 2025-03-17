##############################################################################
#                           Create directory
##############################################################################

# Create output folder 
dir.create(file.path(dir_out, analysis), showWarnings = FALSE)
dir_output <- paste(dir_out, "/", analysis, sep = "")

# File directory
dir_in <- paste(dir_wd, "/05_DEG_ANALYSIS/", analysis_ID2, "/Results/", sep = "")
dir.create(file.path(dir_out, "PCA"))
dir_pcas <- paste(dir_out, "/PCA", sep = "") 


##############################################################################
#                           Set parameters
##############################################################################

# Load log file 
logfile <- read.table(paste(dir_wd, "/log/5_DEG_qc_", analysis_ID2, ".log", sep = ""), header = TRUE)

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}

# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- logfile$condition

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- as.numeric(unlist(str_split(logfile$condition_order, pattern = ",")))

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

# Variance transformation method
md <-  logfile$Variance

# Threshold criteria 
# Significance level
# Default. pvalue < 0.05 
fdr_cutoff <- logfile$fdr_cutoff
# Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$lfc_cutoff

# Threshold label
threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")


# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")))
names(color_list) <- c(trt, "Heatmap")
names(color_list[[trt]]) <- lvl_ord

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)


################################################################################
#                               Load data 
################################################################################


# Load sample information 
sample_info <- read.table(file = paste(dir_in, "Metadata_", project, "_",  analysis_ID2, ".txt", sep = ""), header = TRUE)
rownames(sample_info) <- sample_info$Sample
sample_info <- sample_info[, c("Sample", trt)]
sample_info[[trt]] <- factor(sample_info[[trt]], levels = lvl_ord)


# Phenotypical data
pheno <- data.frame(trt = unique(sample_info[[trt]]))
colnames(pheno) <- trt
rownames(pheno) <- pheno[[trt]]

# Load selected counts  

# From the developers of Mfuzz: "For RNA-seq, you could input log-transformed 
# normalized cpm and then use the standardise function. If you are looking specifically 
# for peaks, you can spare the log transformation, since it reduces the dynamic 
# range... TPM would also be OK, as TPM and CPM only differ by the inverse length 
# of the transcript since Mfuzz::standardise() divides through the standard 
# deviation for each gene, this factor is taken out again."

# User Manual RNfuzzyApp
# Several methods exist for read count normalization and as we implemented 3 different DEA
# methods, several normalization methods are available. If you choose the TCC method, DESeq2
# TMM will be proposed as normalization methods. If DESeq2 [6] is chosen, DESeq2 will be applied
# and finally if edgeR [7] is chosen, TMM RLE, upperquartile will be offered. Once you launch
# normalization and DEG analysis (Figure 7 a), you will receive as results the read counts
# normalized across all samples (Figure 7 b), the results table, with the DE analysis results across all
# conditions (Figure 7 c) and the DEG table, with only those genes that are significantly differentially

if(analysis == "DESeq2" | analysis == "DESeq2_NoFilter"){
  # Load data: log-transformed median of ratios normalized gene counts
  mt <- "Norm_"
  data <- read.delim(file = paste(dir_in, analysis, "_", project, ";Selected_", md, "blindFALSE_", threshold, "_", analysis_ID2, ".txt", sep = ""), header = TRUE, sep = " ")
} else{
  # Load data: normalized counts (CPM)
  mt <- "CPM_"
  data <- read.delim(file = paste(dir_in, analysis, "_", project, ";Selected_CPM_", threshold, "_", analysis_ID2, ".txt", sep = ""), header = TRUE, sep = " ")
}

# Remove column with the comparison name
df <- data %>% distinct(Name, .keep_all = TRUE) %>% select(-1)
df <- df[, -which(colnames(df) %in% c("DEG", "Direction", "logFC", "padj", "shrklogFC", "log2FC", "shrklog2FC", "MeanExp", "lfcSE", "stat", "pvalue"))]

################################################################################
#                             Process data 
################################################################################


# Matrix with the normalized DE gene counts 
sel_counts <- df[,paste(mt, sample_info$Sample, sep = "")]
colnames(sel_counts) <- gsub(mt, "", colnames(sel_counts))
sel_counts <- log2(sel_counts+1)


################################################################################
#                             PCA data 
################################################################################

exc_pca <- createWorkbook()

pcas <- pca_plot(sel_counts, trt, sample_info, color_list)

# Save PCA plots
ggsave(filename = paste("PCA_params_",  analysis, "_Selected_", project, "_", threshold,  "_scree_", analysis_ID, ".pdf", sep = ""), plot = pcas[[1]], path = dir_pcas, height = 4, width = 4, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_", analysis, "_Selected_", project, "_", threshold, analysis_ID, ".pdf", sep = ""), plot = pcas[[2]], path = dir_pcas, height = 5, width = 6, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_", analysis, "_Selected_", project, "_", threshold, analysis_ID, ".pdf", sep = ""), plot = pcas[[3]], path = dir_pcas, height = 5, width = 6, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_", analysis, "_Selected_", project, "_", threshold, analysis_ID, ".pdf", sep = ""), plot = pcas[[4]], path = dir_pcas, height = 5, width = 6, bg = "white")
addWorksheet(exc_pca, "Rotation")
addWorksheet(exc_pca, "GenesPCs")  
writeData(exc_pca, as.data.frame(pcas[[5]]), sheet = "Rotation")
writeData(exc_pca, as.data.frame(pcas[[6]]), sheet = "GenesPCs")

# Save PCA data
saveWorkbook(exc_pca, file =  paste(dir_pcas, "/PCA_data_", analysis, "_", project, ";", md, "blindFALSE_", threshold, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)


