##############################################################################
#                           Create directory
##############################################################################

# Create output folder 
dir.create(file.path(dir_out, analysis), showWarnings = FALSE)
dir_output <- paste(dir_out, "/", analysis, sep = "")

# File directory
dir_in <- paste(dir_wd, "/04_SPLICING_ANALYSIS", analysis_ID2, "/Results/", sep = "")
dir.create(file.path(dir_out, "PCA"))
dir_pcas <- paste(dir_out, "/PCA", sep = "") 


##############################################################################
#                           Set parameters
##############################################################################

# Load log file 
logfile <- read.table(paste(dir_wd, "/log/4_AS_qc_", analysis_ID2, ".log", sep = ""), header = TRUE)

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


# Threshold criteria 
# Significance level
# Default. pvalue < 0.05 
pval_cutoff <- logfile$pval_cutoff
# Log2 fold change threshold
# Default. log2(1.5)
DPSI_cutoff <- logfile$DPSI_cutoff

# Threshold label
threshold <- threshold <- paste("DPSI_", DPSI_cutoff ,"_pval_", pval_cutoff, sep = "")


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
sample_info <- read.table(file = paste(dir_wd, "/Sample_info.csv", sep = ""), header = TRUE, sep = ",")
rownames(sample_info) <- sample_info$Sample
sample_info <- sample_info[, c("Sample", trt)]
sample_info[[trt]] <- factor(sample_info[[trt]], levels = lvl_ord)

if(!is.null(outliers)){sample_info <- sample_info[!which(sample_info$Sample %in% outliers),] }

# Phenotypical data
pheno <- data.frame(trt = unique(sample_info[[trt]]))
colnames(pheno) <- trt
rownames(pheno) <- pheno[[trt]]


# Load selected events: PSI values 
df <- read.delim(file = paste(dir_in, analysis, "_Selected_", project, "_", threshold, "_", analysis_ID, "_Primers.txt", sep = ""), header = TRUE, sep = " ")


################################################################################
#                             Process data 
################################################################################


# Remove specific columns per comparison
df <- df[, -which(colnames(df) %in% c("Significative", "PVAL_WILCOXON_GRPA_VS_GRPB", "ABS_DPSI", "PVAL_WILCOXON", "PSI_GRPB", "PSI_GRPA", "MEDIANPSI_GRPB", "DPSI_GRPA_MINUS_GRPB", "PSIMARGIN_BETWEEN_GRPA_GRPB", "PSIRANGE_GRPA", "AVGPSI_GRPA", 
                                      "MEDIANPSI_GRPA", "MINPSI_GRPB", "MAXPSI_GRPB", "PSIRANGE_GRPB", "AVGPSI_GRPB", "MINPSI_ALL", "MAXPSI_ALL", "PSIRANGE_All", "MINPSI_GRPA", "MAXPSI_GRPA"))]
# Select Exon Skipping events
df <- df[which(df$COMPLEX %in% c("S", "C1", "C2", "C3", "MIC", "ANN")),]

# Remove column with the comparison name
df <- df %>% distinct(Name, .keep_all = TRUE) %>% select(-1)

# Matrix with the normalized DE gene counts 
sel_counts <- df[,paste("PSI_", sample_info$Sample, sep = "")]
colnames(sel_counts) <- gsub("PSI_", "", colnames(sel_counts))



################################################################################
#                             PCA data 
################################################################################

exc_pca <- createWorkbook()

pcas <- pca_plot(sel_counts, trt, sample_info, color_list)

# Save PCA plots
ggsave(filename = paste("PCA_params_",  analysis, "_Selected_", project, "_", threshold,  "_scree.pdf", sep = ""), plot = pcas[[1]], path = dir_pcas, height = 4, width = 4, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_", analysis, "_Selected_", project, "_", threshold, "_EX.pdf", sep = ""), plot = pcas[[2]], path = dir_pcas, height = 5, width = 6, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_", analysis, "_Selected_", project, "_", threshold, "_EX.pdf", sep = ""), plot = pcas[[3]], path = dir_pcas, height = 5, width = 6, bg = "white")
ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_", analysis, "_Selected_", project, "_", threshold, "_EX.pdf", sep = ""), plot = pcas[[4]], path = dir_pcas, height = 5, width = 6, bg = "white")
addWorksheet(exc_pca, "Rotation")
addWorksheet(exc_pca, "GenesPCs")  
writeData(exc_pca, as.data.frame(pcas[[5]]), sheet = "Rotation")
writeData(exc_pca, as.data.frame(pcas[[6]]), sheet = "GenesPCs")

# Save PCA data
saveWorkbook(exc_pca, file =  paste(dir_pcas, "/PCA_data_", analysis, "_Selected_EX_", project, "_", threshold, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)

