################################################################################
#                   ENRICHMENT ANALYSIS FOR TRAJECTORIES 
################################################################################

################################################################################
#                             LOAD LIBRARIES 
################################################################################


library(openxlsx)
library(dplyr)
library(gprofiler2)
library(stringr)



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


################################################################################
#                               SPLICING 
################################################################################

if(analysis == "Splicing"){
  
  # Load log file 
  logfile <- read.table(paste(dir_wd, "/log/4_AS_qc_", logdate, ".log", sep = ""), header = TRUE)
  
  # Input directory 
  dir_trajectories <- paste(dir_out, "/06_TRAJECTORIES/", analysis, "/", clustering, "/Results", sep = "")
  
  # Output directory
  dir_infiles <- paste(dir_out, "/04_SPLICING_ANALYSIS/Results", sep = "")
  
  # Select the identification for the analysis
  label_gen <- "GENE"
  
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
  
  # File name
  if(clustering == "FCM"){
    load_file <- paste(clustering, "_", analysis, "_", project, ";", theshold, "_m_", m,"_c_", c, sep = "")
  } else {load_file <- paste(clustering, "_", analysis, "_", project, ";", theshold, "_c_", c, sep = "")}
  
  
  ################################################################################
  #                                 DEGS 
  ################################################################################
  
} else {
  # Load log file 
  logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)
  
  # Output directory
  dir_infiles <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID,"/Results", sep = "")
  
  # Select the identification for the analysis
  label_gen <- "Ensembl"
  
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

# Load output directory
dir.create(file.path(dir_output, paste(analysis, analysis_ID, sep = "_")), showWarnings = FALSE)
dir_outfolder <- paste(dir_output, "/", analysis, "_", analysis_ID, sep='')
setwd(dir_outfolder)

# Save files with comparison results separate
dir.create(file.path(dir_outfolder , "Results"), showWarnings = FALSE)
dir_res <- paste(dir_outfolder, "/Results", sep='')

# Figures folder
dir.create(file.path(dir_outfolder , "Figures"), showWarnings = FALSE)
dir_fig <- paste(dir_outfolder ,"/Figures", sep='')


# Save files with comparison results separate
dir.create(file.path(dir_outfolder , "GSEA"), showWarnings = FALSE)
dir_gsea <- paste(dir_outfolder, "/GSEA", sep='')




################################################################################
#                          
################################################################################


# Create Workbook
exc <- createWorkbook()

# Assembly
if(organism == "Human"){assembly <- "hsapiens"}else{assembly <- "mmusculus"}


for (i in 1:length(contrast)) {
  
  # Contrast 
  name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(name)
  
  # Select the DEGs
  genes <- bg[which(bg$Comparison == name & bg$DEG == "YES"), "Ensembl"]
  print(length(genes))
  
  # Enrichment analysis
  enrich_dt <- gost(genes,
                    organism = assembly,  # Organism
                    evcodes = TRUE ,      # Save gene names
                    ordered_query = FALSE,
                    significant = FALSE,                     # All results not only statistically significant
                    measure_underrepresentation = FALSE,     # Over-representation or under-representation (TRUE)
                    user_threshold = 0.05,                   # Significance level
                    correction_method = "fdr",               # Multi-testing correction, options "bonferroni" or "fdr"
                    domain_scope = "custom",                 # Domain/background set; Genes with at least one annotation
                    custom_bg = universe,
                    numeric_ns = "",
                    sources = NULL,
                    as_short_link = FALSE,
                    highlight = FALSE,
                    multi_query = FALSE
  )

  # Manhattan plot
  plot_gs <- gostplot(enrich_dt, interactive = FALSE, capped = FALSE)

  results <- enrich_dt$result[order(enrich_dt$result$p_value),]
  results <- results %>% select(term_name, source, term_id, p_value, significant, query, query_size, intersection_size, term_size, effective_domain_size, intersection, precision, recall, source_order, parents)


  addWorksheet(exc, name)
  writeData(exc, name, results)
}

saveWorkbook(exc, file = paste(dir_res, "/", load_file, "_Annotated_", analysis_ID,".xlsx", sep = ""))


################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
log_data <- c()
log_data$Date <- Sys.time()
log_data$analysis_ID <- analysis_ID
log_data$Directory <- dir_out
log_data$Analysis <- analysis
log_data$project_name <- project
log_data$Organism <- organism
log_data$condition <- logfile$condition
log_data$condition_order <- logfile$condition_order
log_data$Outliers <- paste(outliers, collapse = ",") 
log_data$Varexp <- logfile$Varexp
log_data$filter_lab <- filter_lab
log_data$min_count <- logfile$min_count
log_data$min_prop <- logfile$min_prop
log_data$n_large <- logfile$n_large
log_data$min_total <- logfile$min_total
log_data$zscore <- logfile$zscore
log_data$fdr_cutoff <- fdr_cutoff
log_data$lfc_cutoff <- lfc_cutoff
log_data$correction <- logfile$correction
log_data$Variance <- logfile$vsd_type
log_data$zscore <- logfile$zscore


write.table(as.data.frame(log_data), paste(path, project, "/log/7_Enrichment_", analysis, "_", analysis_ID, ".log", sep = ""), row.names = FALSE, eol = "\r")




#####################################################################################################
#                               GSEA DATA PREPARATION
#####################################################################################################


