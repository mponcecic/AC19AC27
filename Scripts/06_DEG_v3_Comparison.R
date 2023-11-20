
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
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
theme_set(theme_DEGs)

# Specie 
# specie = "Human"
specie <- logfile$Organism 

# Summary table
tab_cols <- c("Comparison", "Method", "Genes", "DEGs", "Upregulated", "Downregulated")
sum_tab <- data.frame(matrix(nrow = 0, ncol = length(tab_cols))) 

# Create a list with the different possible method comparisons
comparison <- list()
for (p in 2:length(analysis_list)) {for (t in 1:(p - 1)) {comparison[length(comparison) + 1] <- list(c(analysis_list[t], analysis_list[p]))}}



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
  # Data frame for the correlation plot 
  tab_cor <- c() 
  
  
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
    genes <- genes %>% mutate(Method = analysis) %>% select(GeneID, Method, DEG, logFC, padj, Direction, metadata$Sample)
    data <- rbind(data, genes)

    # Results summary table 
    sum_tab <- rbind(sum_tab, c(name, analysis, nrow(genes), length(genes$DEG == "YES"), length(genes$Direction == "Upregulated"), length(genes$Direction == "Downregulated")))
    
    # Data frame for the correlation plot
    gen_tab <- genes %>% mutate(paste(log_FC, analysis, sep ="_") = logFC) %>% select(GeneID, paste(log_FC, analysis, sep ="_"), paste(padj, analysis, sep ="_"))
    tab_cor <- cbind(tab_cor, gen_tab)
  }
  
  # Rename columns name 
  colnames(sum_tab) <- tab_cols
  
  # Differential expressed genes
  data_de <- data[which(data$DEG == "YES"),]
  
  
  
  ## SCATTERPLOT
  # Select the method comparison and create a scatterplot
  for (k in 1:length(comparison)) {
    comp <- comparison[[k]]
    parm <- c(paste("logFC", comp, sep = "_"), paste("padj", comp, sep = "_"))
    
    # Create a data frame with the following columns
    # - logFC_analysis1
    # - logFC_analysis2
    # - padj_analysis1
    # - padj_analysis2
    # - Significance: The significance of the genes is known, possible options: 
    #   analysis1, analysis2, analysis1_analysis2 and "-" for the non-significant
    #   The significance is based on the logFC and adjusted p-value
    m <- tab_cor[, which(colnames(tab_cor %in% parm))]
    m <- m[match(colnames(m), parm),]
    m$sig <- "-"
    colnames(m) <- c("A", "B", "C", "D", "Significance")
    m$Significance[which(abs(A) > lfc_cutoff & C <= fdr_cutoff)] <-  comp[1]
    m$Significance[which(abs(B) > lfc_cutoff & D <= fdr_cutoff)] <-  comp[2]
    m$Significance[which(abs(A) > lfc_cutoff & C <= fdr_cutoff & abs(B) > lfc_cutoff & D <= fdr_cutoff)] <-  paste(comp, collapse = "_")
    
    # Scatterplot
    ggplot(m, aes(x = A, y = B))+
      geom_point(alpha = 0.4)
    ggsave(paste(dir_fig, "/Corr_", paste(comp, collapse = "_vs_"), "_",  name, "_", project, ".pdf"), plot = last_plot(), height = 4, width = 4, bg = "white")
    
    # Scatterplot with colors
    ggplot(m, aes(x = A, y = B, color = Significance))+
      geom_point(alpha = 0.4)+
      labs(x = comp[1], y = comp[2], title = "log2 Fold Change")+ 
      theme(legend.position = "bottom", legend.box = "horizontal")
    ggsave(paste(dir_fig, "/Corr_", paste(comp, collapse = "_vs_"), "_sig_", name, "_", project, ".pdf"), plot = last_plot(), height = 4, width = 4, bg = "white")
    
  }
  
  
  
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
  
  
  
  ##############################################################################
  #                           Unique Genes
  ##############################################################################
  
  
  # List of unique events
  unq_tab <- dup_tab[dup_tab$Freq == 1,]
  
  # Data frame of unique events
  unq_events <- deg_data[which(data_de$GeneID %in% unq_tab$Var1),]
  
  
  
  ##############################################################################
  #                                 PLOTS
  ##############################################################################
  
  
  ## DENSITY PLOT
  # Representation of the adjusted p-value and log2 fold-change for allS
  A <- ggplot(data, aes(x = padj, fill = Method)) +
    geom_density(alpha = 0.5)+
    scale_fill_manual(values = color_l$Shared)+
    labs(x = "adjusted p-value", y = "Counts")
  B <- ggplot(data, aes(x = logFC, fill = Method)) +
    geom_density(alpha = 0.5)+
    scale_fill_manual(values = color_l$Shared)+
    labs(x = "log2FC", y = "Counts")
  fig <- ggarrange(A, B, ncol = 2, nrow = 1, widths = 10, heights = 5)
  ggsave(filename = paste("00_Histogram_verif_", ref, ".pdf", sep = ""), path = dir_fig, plot = fig, width = 4, height = 2, bg = "white")
  
  
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
  
  
}



# Summary table 
write.csv(sum_tab, paste(dir_outfolder, "/Comparison_result_table_", project,".csv", sep = ""), row.names = FALSE)



################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- dir_out







