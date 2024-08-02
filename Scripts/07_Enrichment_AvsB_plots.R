
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

# Comparison
comparison <- "noDOXvsDOX24H"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------



################################################################################
#                         SET WORKNING DIRECTORY
################################################################################


# General output directory
dir_out <- paste(path, project, sep = "")

# Directory files
dir_files <- paste(dir_out, "/07_ENRICHMENT/", analysis, "_", analysis_ID, "/Results/", sep = "")

# Directory figures
dir_fig <- paste(dir_out, "/07_ENRICHMENT/", analysis, "_", analysis_ID, "/Figures/", sep = "")
setwd(dir_fig)


################################################################################
#                               LOAD FILE 
################################################################################


source(paste(dir_out, "/utils/libraries_degs.R", sep = ""))

# Theme ggplot2  
theme_splicing <- theme_bw()+ 
  theme(axis.text.x = element_text(color = "black", size = 6), axis.text.y = element_text(color = "black", size = 6), axis.ticks = element_line(color = "black", linewidth = 0.3), 
        legend.title = element_text(size = 8), legend.text = element_text(size = 6), panel.grid = element_line(size = 0.3))
theme_set(theme_splicing)



################################################################################
#                               LOAD FILE 
################################################################################


# Load log file 
logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)

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

# Organism
organism <- logfile$Organism

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))

# Threshold criteria 
# Significance level
# Default. pvalue < 0.05 
pval_cutoff <- logfile$pval_cutoff
# Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$DPSI_cutoff

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <- paste(unlist(str_split(logfile$Outliers, ",")), collapse = "_")}

# Percentaje
if(is.null(logfile$Percentage)){perc <- NULL}else{perc <- logfile$Percentage}

# Threshold label
threshold <- paste("padj_", pval_cutoff, "_log2FC_", lfc_cutoff, "_Annotated_", analysis_ID, sep = "")

# Load annotated file
data <- read.xlsx(paste(dir_files, list.files(dir_files), sep = ""), sheet = comparison)
print(dim(data))

# select only significative pathways
data <- data[data$si == TRUE, ]

################################################################################
#                               PLOT 
################################################################################

# Select the firs 15 pathways
data <- data[order(data$p_value, decreasing = TRUE),]
dt <- data
path_order <- dt$term_name

# Dotplot Precision 
ggplot(data = dt, aes(x = p_value, y = factor(term_name, levels = path_order), color = p_value, size = precision))+
  geom_point()+
  scale_color_gradient(name = "p-value", high = "blue", low = "red")+
  scale_size(name = "Precision")+ 
  labs(x = "", y = "")
ggsave(filename = paste("00_Dotplot_", project, "_", comparison, "_", threshold, ".pdf", sep = ""), height = 4, width = 4.5, plot = last_plot())


# Dotplot Intersection size
ggplot(data = dt, aes(x = p_value, y = factor(term_name, levels = path_order), color = p_value, size = intersection_size))+
  geom_point()+
  scale_color_gradient(name = "p-value", high = "blue", low = "red")+
  scale_size(name = "Intersection size")+ 
  labs(x = "", y = "")
ggsave(filename = paste("01_Dotplot_", project, "_", comparison, "_", threshold, ".pdf", sep = ""), height = 4, width = 4.5, plot = last_plot())


