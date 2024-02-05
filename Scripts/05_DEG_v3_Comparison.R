
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
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/user/"
path <- "W:/user/"

# Date of the log file 5_DEG_qc_xxxxx.log
logdate <- "20231204"

# Select the methods you want to compare
# analysis_list <- c("DESeq2", "EdgeR", "limma-voom", "Wilcoxon")
analysis_list <- c("EdgeR", "limma-voom")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/libraries_degs.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file
logfile <- read.table(paste(path, project, "/log/5_DEG_qc_", logdate, ".log", sep = ""), header = TRUE)


# Main directory
# dir_out <- paste(path, project, sep = "")   # Default option
dir_out <- paste(path, project, "/05_DEG_ANALYSIS", sep = "")

# Output directory 
dir_output <- paste(dir_out, "/Results", sep = "")

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
var_exp <- logfile$covariance

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

# Create a workbook
exc <- createWorkbook()


################################################################################
#                               LOAD DATA 
################################################################################


## Load metadata file 
sample_info <- read.table(paste(dir_output, "/Metadata_", project, ".txt", sep = ""))



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
  comp_lvl <- c(control, experimental)
  
  # Select the contrast levels and method colors
  color_l <- color_list
  color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% c(experimental, control))]
  color_l[["Shared"]] <- color_l[["Shared"]][which(names(color_l[["Shared"]]) %in% analysis_list)]
  
  # Threshold label
  threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")
  
  # Data frame with all the methods information 
  data <- c() 
  
  # Data frame which will gather the results per each comparison and all the 
  # methods used. In addition, it's used for the correlation plots of the logFC 
  # and adjusted p-values
  result_tab <- c() 
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  # Load input directory
  dir_outfolder <- paste(dir_out, "/", name, sep='')
  dir_infiles <- paste(dir_outfolder, "/Results", sep = "")
  
  # Method comparison figures folder
  dir.create(file.path(dir_outfolder , "Comparison"), showWarnings = FALSE)
  dir_fig <- paste(dir_outfolder ,"/Comparison", sep='')
  
  
  ##############################################################################
  #                               Load Data
  ##############################################################################
  
  
  # Load sample information per comparison
  metadata <- sample_info[which(sample_info[[trt]] %in% comp_lvl),]
  metadata[,trt] <- factor(metadata[,trt], levels = c(control, experimental))
  
  # Load QC_result_name_project.txt
  file_trs <- read.table(paste(dir_output,"/QC_result_", project,".txt", sep = ""), header = TRUE)
  md <- file_trs$Transformation
  
  # Load differentially expressed genes results
  for (j in 1:length(analysis_list)) {
    
    # Select analysis
    analysis <- analysis_list[j]
    
    # Reference
    ref <- paste(analysis, name ,project, sep = "_")
    
    # Load results
    # Select the statistic columns per each analysis
    if(analysis == "DESeq2" | analysis == "DESeq2_NoFilter"){
      df <- read.table(paste(dir_infiles, "/", ref, ";All_", md, "blindFALSE_", threshold,".txt", sep = ""), header = TRUE)
      col_nam <- c("logFC", "padj", "shrklogFC", "MeanExp", "lfcSE", "stat", "pvalue")
    } else if(analysis == "EdgeR"){
      df <- read.table(paste(dir_infiles, "/", ref, ";All_CPM_", threshold,".txt", sep = ""), header = TRUE)
      col_nam <- c("logFC", "padj", "logCPM", "pvalue")
    } else if(analysis == "limma-voom"){
      col_nam <- c("logFC", "pvalue", "logCPM", "t", "padj", "B")
      df <- read.table(paste(dir_infiles, "/", ref, ";All_CPM_", threshold,".txt", sep = ""), header = TRUE)
    } else {
      col_nam <- c("logFC", "pvalue", "padj")}

     
    genes <- df %>% mutate(Method = analysis) %>% select(Name, Symbol, Ensembl, Method, DEG, logFC, padj, Direction, metadata$Sample)
    data <- rbind(data, genes)

    # Results summary table 
    sum_tab <- rbind(sum_tab, c(name, analysis, nrow(genes), length(which(genes$DEG == "YES")), length(which(genes$Direction == "Upregulated")), length(which(genes$Direction == "Downregulated"))))
    
    # Rename statistics with the method used
    colnames(df)[which(colnames(df) %in% col_nam)] <- paste(col_nam, analysis, sep = "_")
    
    if(j==1){
      # Data frame for the correlation plot
      gen_tab <- df %>% select(Name, Symbol, Ensembl, paste(col_nam, analysis, sep ="_"), everything())
      gen_tab[analysis] <- ifelse(gen_tab$DEG == "YES", "X", "")
      result_tab <- gen_tab
    } else {
      # Data frame for the correlation plot
      gen_tab <- df %>% select(Name, paste(col_nam, analysis, sep ="_"), DEG)
      gen_tab[analysis] <- ifelse(gen_tab$DEG == "YES", "X", "")
      gen_tab <- gen_tab[,-which(colnames(gen_tab) %in% "DEG")]
      result_tab <- merge(result_tab, gen_tab)}
  }
  
  # Remove DEG column and direction column
  result_tab <- result_tab[,-which(colnames(result_tab) %in% c("DEG", "Direction"))]
  
  
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
    m <- result_tab[,which(colnames(result_tab) %in% parm)]
    m <- m[,match(colnames(m), parm)]
    m$sig <- "-"
    colnames(m) <- c("A", "B", "C", "D", "Significance")
    m$Significance[which(abs(m$A) > lfc_cutoff & m$C <= fdr_cutoff)] <-  comp[1]
    m$Significance[which(abs(m$B) > lfc_cutoff & m$D <= fdr_cutoff)] <-  comp[2]
    m$Significance[which(abs(m$A) > lfc_cutoff & m$C <= fdr_cutoff & abs(m$B) > lfc_cutoff & m$D <= fdr_cutoff)] <-  paste(comp, collapse = "_")
    
    # Scatterplot
    ggplot(m, aes(x = A, y = B))+
      geom_point(alpha = 0.4)+
      labs(x = comp[1], y = comp[2], title = "log2 Fold Change")
    ggsave(paste(dir_fig, "/Corr_", paste(comp, collapse = "_vs_"), "_",  name, "_", project, ".pdf", sep = ""), plot = last_plot(), height = 4, width = 4, bg = "white")
    
    # Scatterplot with colors
    ggplot(m, aes(x = A, y = B, color = Significance))+
      geom_point(alpha = 0.4)+
      labs(x = comp[1], y = comp[2], title = "log2 Fold Change", color = "")+ 
      theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 6))
    ggsave(paste(dir_fig, "/Corr_Sig_", paste(comp, collapse = "_vs_"), "_", name, "_", project, ".pdf", sep = ""), plot = last_plot(), height = 4, width = 4, bg = "white")
    
  }
  
  
  ##############################################################################
  #                         Shared Genes
  ##############################################################################
  
  
  # List of shared events
  genes_list <- split(data$Ensembl, data$Method)
  de_list <- split(data_de$Ensembl, data_de$Method)
  
  
  # List of the duplicated events
  dup_tab <- as.data.frame(table(data_de$Ensembl))
  dup_gen <- dup_tab[dup_tab$Freq > 1,]
  
  
  # Duplicate events per tissue
  #   - Times repeated per group
  #   - Groups in which appear
  #   - Direction of the event (up/down)
  #   - Congruence: Direction the event was the same or not in 
  #     both groups
  duplicates <- data_de[which(data_de$Ensembl %in% dup_gen$Var1),] %>% 
    group_by(Ensembl) %>%
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
  unq_events <- data_de[which(data_de$Ensembl %in% unq_tab$Var1),]
  
  
  
  ##############################################################################
  #                                 PLOTS
  ##############################################################################
  
  
  ## DENSITY PLOT
  # Representation of the adjusted p-value and log2 fold-change for allS
  A <- ggplot(data, aes(x = padj, fill = Method)) +
    geom_density(alpha = 0.5)+
    scale_fill_manual(values = color_l$Shared)+
    labs(x = "adjusted p-value", y = "Counts")+
    theme(legend.position = "none")
  B <- ggplot(data, aes(x = logFC, fill = Method)) +
    geom_density(alpha = 0.5, )+
    scale_fill_manual(values = color_l$Shared)+
    labs(x = "log2FC", y = "Counts")+
    theme(legend.position = "none")
  fig <- ggarrange(A, B, ncol = 2, nrow = 1, widths = 10, heights = 5)
  ggsave(filename = paste("00_Histogram_verif_", paste(analysis_list, collapse = "_"), ".pdf", sep = ""), path = dir_fig, plot = fig, width = 4, height = 2, bg = "white")
  
  
  ## VENN DIAGRAM
  # Visualization of the shared genes among the different methods used
  ggvenn(de_list,
         digits = 2,                                  # Two decimals
         fill_color = as.vector(color_l[["Shared"]]),                        # Fill color
         fill_alpha = 0.5,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 3,                           # Default 
         text_size = 2.5                              # Default
  )
  ggsave(filename = paste("Venn_diagram_", paste(analysis_list, collapse = "_"), "_", name, "_", project, ".pdf", sep =""), path = dir_fig, height = 4, width = 4, plot = last_plot(), bg = "white")
  
  
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
  pdf(file = paste(dir_fig, "/Upset_plot_", paste(analysis_list, collapse = "_"), "_", name, "_", project,".pdf", sep = ""), height = 5, width = 10, bg = "white")
  print(upset1)
  dev.off()
  
  
  
  ##############################################################################
  #                                 SAVE DATA  
  ##############################################################################
  
  # Order output
  result_tab <- result_tab %>% select(Name, Symbol, Ensembl, analysis_list, everything())
  # write.table(result_tab, paste(dir_infiles, "/", name, ";", paste(analysis_list, collapse = "_"), "_", threshold, ".txt", sep = ""))
  # write.xlsx(result_tab, paste(dir_infiles, "/", name, ";", paste(analysis_list, collapse = "_"), "_", threshold, ".xlsx", sep = ""), overwrite = TRUE)
  
  # Save data in the workbook
  addWorksheet(exc, name)
  writeData(exc, result_tab, sheet = name)
  
  }

# Save workbook
saveWorkbook(exc, file =  paste(dir_output, "/", project, ";", paste(analysis_list, collapse = "_"), "_", threshold, ".xlsx", sep = ""), overwrite = TRUE)

# Summary table 
write.csv(sum_tab, paste(dir_output, "/Comparison_result_table_", project, ".csv", sep = ""), row.names = FALSE)


################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- dir_out
log_data$AnalysisComp <- analysis_list

write.table(result_tab, paste(path, project, "/log/5_DEG_v3_", paste(analysis_list, collapse = "_"), "_", logdate, ".log", sep = ""))







