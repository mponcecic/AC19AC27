

################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC19AC27"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/user/"
path <- "W:/mponce/"

# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID <- "20240801092404"

# Select analysis/method performed
analysis <- "DESeq2"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# General output directory
dir_out <- paste(path, project, sep = "")

# Load libraries
source(paste(path, project, "/utils/libraries_degs.R", sep = ""))
# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))

# Load log file 
logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)

# Output directory
dir_infiles <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID,"/Results", sep = "")

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


# Significance level
fdr_cutoff <- logfile$fdr_cutoff

# Log2 fold change threshold
lfc_cutoff <- round(logfile$lfc_cutoff, 2)

# Organism
organism <- logfile$Organism

# Threshold label
theshold <- paste("padj_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep = "")

# Variance stabilization
md <- logfile$Variance

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))


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

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){
  outliers <- NULL
} else {
  outliers <-  paste(logfile$Outliers, collapse = ",")
  project <- paste(project, paste(outliers, collapse = "_"), sep = "_")}

#Label of filtering
filter_lab <- logfile$filter_lab

################################################################################
#                         Load data  
################################################################################

# File name 
load_file <- paste(analysis, "_", project, ";All_", md, "blindFALSE_", theshold, sep = "")

# Load universe
bg <- read.table(paste(dir_infiles, "/", load_file, "_", analysis_ID, ".txt", sep = ""), header = TRUE)
  
# Load metadata file 
sample_info <- read.table(paste(dir_infiles, "/Metadata_", project, "_", analysis_ID, ".txt", sep = ""))



################################################################################
#                         noDOX vs DOX24H 
################################################################################

# Comparison
comparison <- "noDOXvsDOX24H"
lvl_ord <- c("noDOX", "DOX24H")

# Figure directory
dir_fig <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID,"/", comparison, "/DESeq2/", sep = "")

# Color list
color_l <- color_list
color_l$Condition <- color_l$Condition[names(color_l$Condition) %in% lvl_ord] 

# Process data 
metadata <- sample_info[sample_info$Condition %in% lvl_ord, ]
data <- bg[which(bg$Comparison == comparison & bg$DEG == "YES"), c("Name", paste("VST_", metadata$Sample, sep = ""))]
colnames(data)[-1] <- metadata$Sample
rownames(data) <- data$Name
m <- data[,-1]
print(head(m))

plot_hm <- heatmap_plot(m, metadata = metadata, trt, color_l, callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")})

pdf(paste(dir_fig, "/Heatmap_zscore_DESeq2_noDOXvsDOX24H_AC19AC27_redone_20240801092404.pdf", sep = ""), height = 4, width = 4, bg = "white")
print(plot_hm)
dev.off()



################################################################################
#                         noDOX vs DOX42H 
################################################################################

# Comparison
comparison <- "noDOXvsDOX42H"
lvl_ord <- c("noDOX", "DOX42H")

# Figure directory
dir_fig <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID,"/", comparison, "/DESeq2/", sep = "")

# Color list
color_l <- color_list
color_l$Condition <- color_l$Condition[names(color_l$Condition) %in% lvl_ord] 

# Process data 
metadata <- sample_info[sample_info$Condition %in% lvl_ord, ]
data <- bg[which(bg$Comparison == comparison & bg$DEG == "YES"), c("Name", paste("VST_", metadata$Sample, sep = ""))]
colnames(data)[-1] <- metadata$Sample
rownames(data) <- data$Name
m <- data[,-1]
print(head(m))

plot_hm <- heatmap_plot(m, metadata = metadata, trt, color_l, callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")})

pdf(paste(dir_fig, "/Heatmap_zscore_DESeq2_noDOXvsDOX42H_AC19AC27_redone_20240801092404.pdf", sep = ""), height = 4, width = 4, bg = "white")
print(plot_hm)
dev.off()