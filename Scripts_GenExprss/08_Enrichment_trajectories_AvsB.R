################################################################################
#                ENRICHMENT ANALYSIS FOR TRAJECTORIES
################################################################################

# Date: Mon May 19 09:31:19 2025
# Author: Maria Ponce

# Summary
# ---------------------

# Enrichment analysis per cluster in which we use the whole transcriptome as the 
# background/universe per each cluster. This script can only be run in DEGs to run 
# the code in Splicing data some changes must be done.
# 
# In addition, you have a few parameters to modify:
#   - mem_list: Vector of the membership values to use a cutoff. Only use this 
#     option when running FCM (Mfuzz) data.


# Input folder: W:/mponce/project/06_TRAJECTORIES/analysis_ID_traj/Results/
# Output folder: W:/mponce/project/07_ENRICHMENT/analysis_ID_enrich/



################################################################################
#                             LOAD LIBRARIES 
################################################################################


library(openxlsx)
library(dplyr)
library(gprofiler2)



################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC19AC27"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "W:/mponce/"

# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID_degs <- "20240801092404"
analysis_ID_traj <- "20250319124539"

# Select analysis/method performed
analysis <- "DESeq2"

# Method used for trajectories
clustering <- "K-means"

# Optimal number of clusters
c <- 5

# Fuzziness coefficient
m <- NULL 

# Membership/Centroid distance 
# Filter out the events/genes whit a membership value below this value
# mem_list <- c(0.8, 0.9, 0.95)
mem_list <- 0
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# General output directory
dir_out <- paste(path, project, sep = "")

# Analysis ID 
# Created using a timestap to trace all de outputs belonging to a certain setup
analysis_ID <- format(Sys.time(), "%Y%m%d%H%M%S")



################################################################################
#                               SPLICING 
################################################################################

if(analysis == "Splicing"){
  
  
  ################################################################################
  #                                 DEGS 
  ################################################################################
  
} else {
  # Load log file 
  logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID_degs, ".log", sep = ""), header = TRUE)
  
  # Input directory 
  dir_trajectories <- paste(dir_out, "/06_TRAJECTORIES/", analysis_ID_traj, "/", analysis, "/", clustering, "/Results", sep = "")
  
  # Output directory
  dir_infiles <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID_degs, "/", "Results", sep = "")
  
  # Significance level
  fdr_cutoff <- logfile$fdr_cutoff
  
  # Log2 fold change threshold
  lfc_cutoff <- round(logfile$lfc_cutoff, 2)
  
  # Threshold label
  theshold <- paste("padj_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep = "")
  
  # Variance stabilization
  md <- logfile$Variance
  
  # File name
  if(clustering == "Mfuzz"){
    load_file <- paste(clustering, "_", analysis, "_", project, ";Norm_", theshold, "_m_", m,"_c_", c, "_", analysis_ID_traj, sep = "")
  } else {load_file <- paste(clustering, "_", analysis, "_", project, ";Norm_", theshold, "_c_", c, "_", analysis_ID_traj, sep = "")}
  
  
  
  ################################################################################
  #                         Load data  
  ################################################################################
  
  # Load universe
  bg <- read.table(paste(dir_infiles, "/", analysis, "_", project, ";All_", md, "blindFALSE_", theshold, "_", analysis_ID_degs, ".txt",  sep = ""), header = TRUE)
  
  # Select unique Ensembl ids
  universe <- unique(bg$Ensembl)
  
}


##############################################################################
#                         Create working directories
##############################################################################

# Create a folder for the enrichment analysis results. This folder contains 
# several folders classified in Results and Figures. 

# Create enrichment directory
dir.create(file.path(dir_out , "07_ENRICHMENT"), showWarnings = FALSE)
dir_output <- paste(dir_out, "/07_ENRICHMENT", sep='')
dir.create(file.path(dir_output, analysis_ID), showWarnings = FALSE)
dir_output <- paste(dir_output, "/", analysis_ID, sep='')

# Load output directory
dir.create(file.path(dir_output, analysis), showWarnings = FALSE)
dir_outfolder <- paste(dir_output, "/", analysis, sep='')
setwd(dir_outfolder)

# Save files with comparison results separate
dir.create(file.path(dir_outfolder , "Results"), showWarnings = FALSE)
dir_res <- paste(dir_outfolder, "/Results", sep='')

# Save figures
dir.create(file.path(dir_outfolder , "Figures"), showWarnings = FALSE)
dir_fig <- paste(dir_outfolder, "/Figures", sep='')



################################################################################
#                         ENRICHMENT ANALYSIS
################################################################################



for (j in 1:length(mem_list)){
  
  # Membership cutoff
  mem_cutoff <- mem_list[j]
  
  # Create Workbook
  exc <- createWorkbook()
  
  # Data frame with 
  sumy <- c()
  
  for (i in 1:c) {
    sheet <- i
    data <- read.xlsx(paste(dir_trajectories, "/", load_file, ".xlsx", sep = ""), sheet = sheet)
    
    if(clustering == "Mfuzz"){
      idx <- which(data$Membership >= mem_cutoff)
      genes <- list(unique(data$Ensembl[idx]))
      a <- c(length(idx), dim(data)[1])
      sumy <- rbind(sumy, a)
      
    } else {genes <- list(unique(data$Ensembl))}
    
    
    
    
    # Enrichment analysis
    enrich_dt <- gost(genes, 
                      organism = "hsapiens",                   # Organism
                      ordered_query = FALSE,
                      significant = TRUE,                      # All results not only statistically significant
                      measure_underrepresentation = FALSE,     # Over-representation or under-representation (TRUE)
                      user_threshold = 0.05,                   # Significance level
                      correction_method = "fdr",               # Multi-testing correction, options "bonferroni" or "fdr"
                      domain_scope = "custom",                 # Domain/background set; Genes with at least one annotation
                      custom_bg = universe,
                      numeric_ns = "", 
                      sources = NULL, 
                      as_short_link = FALSE,
                      highlight = FALSE,
                      multi_query = FALSE, 
                      evcodes = TRUE                      #
    )
    
    # Manhattan plot
    plot_gs <- gostplot(enrich_dt, interactive = FALSE, capped = FALSE)
    # publish_gostplot(p = plot_gs, filename = paste(dir_fig, "/Enrichment_plot_", analysis, "_", theshold, "_c_", c, "_", analysis_ID , ".pdf", sep = ""), width = 8, height = 6)
    pdf(paste(dir_fig, "/Enrichment_plot_", analysis, "_", theshold, "_c_", i, "_m_", mem_cutoff, "_", analysis_ID , ".pdf", sep = ""), width = 8, height = 6, bg = "white")
    print(plot_gs)
    dev.off()
    
    results <- enrich_dt$result[order(enrich_dt$result$p_value),]
    results <- results %>% select(term_name, source, term_id, p_value, significant, query, query_size, intersection_size, term_size, effective_domain_size, intersection, precision, recall, source_order, parents)

    addWorksheet(exc, sheet)
    writeData(exc, sheet, results)
  }
  
  saveWorkbook(exc, file = paste(dir_res, "/", load_file, "_Annotated_m_", mem_cutoff, ".xlsx", sep = ""), overwrite = TRUE)

  if(clustering == "Mfuzz"){
    colnames(sumy) <- c(paste("Membership>=", sep = ""), "Total")
    rownames(sumy) <- 1:c
    # write.csv(sumy, paste(dir_res, "/", load_file, "_Summary.csv", sep = ""), row.names = TRUE)
    write.csv(sumy, paste(dir_res, "/", load_file, "_Summary_m_", mem_cutoff,".csv", sep = ""), row.names = TRUE)

  } else{}
  
  
}

################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- analysis_ID
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- dir_out
log_data$Analysis <- analysis
log_data$project_name <- project
log_data$Organism <- logfile$Organism
log_data$condition <- logfile$condition
log_data$condition_order <- logfile$condition_order
log_data$Cluster <- c
if(analysis == "Splicing"){
  log_data$analysis_ID_degs <- analysis_ID_degs
  log_data$dpsi_cutoff <- dpsi_cutoff
  log_data$pval_cutoff <- pval_cutoff
} else { 
  log_data$analysis_ID_splicing <- analysis_ID_degs
  log_data$fdr_cutoff <- fdr_cutoff
  log_data$lfc_cutoff <- lfc_cutoff
  log_data$Variance <- md
}
log_data$mem_cutoff <- paste(mem_cutoff, collapse = ",")
log_data$analysis_ID_traj <- analysis_ID_traj

write.table(as.data.frame(log_data), paste(path, project, "/log/7_Enrichment_", analysis, "_", analysis_ID, ".log", sep = ""), row.names = FALSE, eol = "\r")




