
################################################################################
#                       SHARED GENES BETWEEN METHODS
################################################################################


# Summary
# ---------

# Compare the differentially expressed genes among the different methods. 




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
logdate <- "20231110"

# Select the methods you want to compare
# analysis_list <- c("DESeq2", "EdgeR", "limma-voom", "Wilcoxon")
analysis_list <- c("DESeq2", "EdgeR", "limma-voom")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/1_DEG_qc_", logdate, ".log", sep = ""), header = TRUE)


# Output directory
# dir_out <- paste(path, project, sep = "")   # Default option
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
# var_exp <- NULL
var_exp <- c("RIN")

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}


# Filtering lowly expressed genes
# Default: TRUE
filter_cutoff <- logfile$filter_cutoff


### Threshold criteria 

## Significance level
# Default. pvalue < 0.05 
fdr_cutoff <- logfile$fdr_cutoff

## Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$lfc_cutoff

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- logfile$correction

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")),
                   Direction = unlist(str_split(logfile$colordir, pattern = ",")),
                   Shared = unlist(str_split(logfile$colorsh, pattern = ",")))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")
names(color_list[[trt]]) <- lvl_ord
names(color_list[["Direction"]]) <- c("Downregulated", "Not significant", "Upregulated")
names(color_list[["Shared"]]) <- c("DESeq2", "EdgeR", "limma-voom", "Wilcoxon")

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Specie 
# specie = "Human"
specie <- logfile$Organism 



##############################################################################
#                         SHARED EVENTS
##############################################################################



# Annotation for all the results 
ref <- paste(analysis, "_", name, "_", project, sep = "")

paste(dir_output, "/", ref, ";All_", md, "blindFALSE_", threshold,".txt", sep = "")




################################################################################
#                               COMPARISONS
################################################################################


for (i in 1:length(contrast)){
  
  
  # Contrast 
  name <- paste(contrast[[i]][2], "vs", contrast[[i]][3], sep = "")
  print(name)
  
  # Contrast levels 
  control <- contrast[[i]][3]
  experimental <- contrast[[i]][2]
  
  # Select the contrast levels and method colors
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% c(experimental, control))]
  color_l[["Shared"]] <- color_l[["Shared"]][which(names(color_l[["Shared"]]) %in% analysis_list)]
  
  # Data frame with all the methods information 
  data <- c() 
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  
  # Load output directory
  dir_outfolder <- paste(dir_out, "/05_DEG_ANALYSIS/", name, sep='')
  setwd(dir_outfolder)
  
  # Result folder
  dir_output <- paste(dir_outfolder, "/Results", sep='')
  
  # Method comparison figures folder
  dir.create(file.path(dir_outfolder , "Method_comparison"), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder ,"/Method_comparison", sep='')
  
  
  ##############################################################################
  #                               Load Data
  ##############################################################################
  
  
  # Load sample information per comparison
  metadata <- read.table(paste(dir_output, "/Metadata_", name, "_", project, ".txt", sep = ""))
  metadata[,trt] <- factor(metadata[,trt], levels = c(control, experimental))
  
  # Load differentially expressed genes results
  for (j in 1:length(analysis_list)) {
    
    # Select analysis
    analysis <- analysis_list[j]
    
    # Results annotation 
    ref <- paste(analysis, "_", name, "_", project, sep = "")
    
    # Load results
    genes <- read.table(paste(dir_output, "/", ref, ";All_", md, "blindFALSE_", threshold,".txt", sep = ""))
    genes <- genes %>% mutate(Method = analysis) %>% select(GeneID, Method,  DEG, logFC, padj, Direction, metadata$Sample)
    
    data <- rbind(data, genes)
  }
  
  # Differential expressed genes
  data_de <- data[which(data$DEG == "YES"),]
  
  
  
  ##############################################################################
  #                         Shared Genes
  ##############################################################################
  
  
  # List of shared events
  genes_list <- split(data$GeneID, data$Method)
  de_list <- split(data_de$GeneID, data_de$Method)
  
  
  # List of the duplicated events
  dup_tab <- as.data.frame(table(data_de$GeneID))
  dup_gen <- dup_tab[dup_tab$Freq > 1,]
  
  
  # Duplicate events per tissue
  #   - Times repeated per group
  #   - Groups in which appear
  #   - Direction of the event (up/down)
  #   - Congruence: Direction the event was the same or not in 
  #     both groups
  duplicates <- data_de[which(data_de$GeneID %in% dup_gen$Var1),] %>% 
    group_by(GeneID) %>%
    summarise(N = n_distinct(Method),
              Method = paste(unique(Method), collapse = ", "),
              Direction = paste(Direction, collapse = ", ")) %>%
    filter(N > 1)
  
  # Congruence add to the duplicates
  duplicates$Congruence <- sapply(duplicates$Direction, check_congruence)
  
  
  
  ## VENN DIAGRAM
  # Visualization of the shared genes among the different methods used
  ggvenn(de_list,
         digits = 2,                                  # Two decimals
         fill_color = color_l,                        # Fill color
         fill_alpha = 0.5,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 3,                           # Default 
         text_size = 2.5                              # Default
  )
  ggsave(filename = paste("Venn_diagram_", name, "_", project, ".pdf", sep =""), path = dir_sh, height = 4, width = 4, plot = last_plot(), bg = "white")
  
  
  ## UPSET PLOT 
  # Visualization of the shared genes among the different methods used
  upset1 <- upset(fromList(de_list),                       
                  order.by = "freq",                          # Order by frequency
                  nsets = length(de_list),                  # Number of sets
                  keep.order = TRUE,                          
                  mainbar.y.label = "Shared genes",       # Y axis label
                  sets.x.label = "Total genes",                    # X axis label
                  main.bar.color = "#2E8B57",      # Y axis bar color
                  sets.bar.color = "#483D8B",      # X axis bar color
                  matrix.color =  "black"        # Matrix dots color
  )     
  pdf(file = paste(dir_sh, "/Upset_plot_", name, "_", project,".pdf", sep = ""), height = 5, width = 10, bg = "white")
  print(upset1)
  dev.off()
  
  
  
  ##############################################################################
  #                           Unique Genes
  ##############################################################################
  
  
  # List of unique events
  unq_tab <- dup_tab[dup_tab$Freq == 1,]
  
  # Data frame of unique events
  unq_events <- deg_data[which(data_de$GeneID %in% unq_tab$Var1),]
  
}



################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- dir_out







