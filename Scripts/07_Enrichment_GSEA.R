
################################################################################
#                           GSEA preparation files
################################################################################

################################################################################
#                             LOAD LIBRARIES 
################################################################################


library(dplyr)
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
  
}

# Load for GSEA files creation
metadata <- read.csv(paste(dir_out, "/Sample_info.csv", sep = ""))
metadata <- metadata[match(metadata$Sample, rev(metadata$Sample)),]


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

# GSEA output directory
dir_gsea <- paste(dir_outfolder, "/GSEA", sep='')


##############################################################################
#                         Save files
##############################################################################

for (i in 1:length(contrast)) {
  
  # Levels
  experimental <- contrast[[i]][3]  
  control <- contrast[[i]][2]
  
  # Contrast 
  name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(name)
  
  # Sample info
  sample_info <- metadata[which(metadata$Condition %in% c(control, experimental)),]
  
  # Select all the genes in the analysis
  subgenes <- bg[which(bg$Comparison == name), ]
  subgenes <- subgenes[!is.na(subgenes$Ensembl),]

  # Create expression data matrix
  dt <- subgenes[, paste("Norm_", sample_info$Sample, sep = "")]
  colnames(dt) <- gsub("Norm_", "", colnames(dt))
  dt <- dt[, match(sample_info$Sample, colnames(dt))]
  dt <- dt %>% mutate(NAME = subgenes$Ensembl, DESCRIPTION = "na") %>% select(NAME, DESCRIPTION, everything())
  write.table(dt, file = paste(dir_gsea, "/GSEA_", name, "_expressiondata_", analysis_ID, ".txt", sep = ""), row.names = F, quote = FALSE, sep = "\t")
  
  ### PHENO DATA
  phenodata <- c()
  phenodata$first <- paste(dim(sample_info)[1], 2, 1, collapse = " ")
  phenodata$second <- paste("#" , experimental, control, collapse = " ")
  phenodata$third <- paste(sample_info$Condition, collapse = " ")
  write.table(t(data.frame(phenodata)), file = paste(dir_gsea, "/GSEA_", name, "_phenodata_", analysis_ID, ".cls", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}





