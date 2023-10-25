
################################################################################
#                               METADATA                          
################################################################################


# Summary
# --------
#
# Generate a sample information file also refered as metadata. The sample 
# information file is used since the beginning of the analysis to generate 
# different scripts and in the differential expression analysis which is based 
# in the project information. Allowing the execution of the scripts with minor 
# adjustments.
# 
# This file must be MANUALLY adjusted for the different projects. 
# 
# Content
#   - Named Sample_info.csv
#   - Must have at least, the following variables:
#     - Sample: Variable with the samples names 
#     - Treatment: Specific treatment per sample. Must be a factor variable.



################################################################################
#                               LOAD LIBRARIES                           
################################################################################


suppressMessages(library(stringr)) 
suppressMessages(library(pdftools)) 
suppressMessages(library(tidyverse)) 

################################################################################
#                                 LOAD DATA                         
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC58"

# Pathway to the folders and files
# Can be your personal folder in BigData
path <- "W:/mponce/"

# Input directory. Raw gene counts  
dir_in <- "W:/DATA_shared/AC-58_TotalRNAseq/"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Output directory
dir_out <- paste(path, project, sep = "")

# Set working directory
setwd(dir_out)

# List of samples names
samples <- unique(gsub("_1.fastq.gz","", list.files(path = paste(dir_in, "FASTQs/", sep = ""), pattern = "_1.fastq.gz")))

file <- list.files(path = dir_in, pattern = "Library_Preparation")


file <- pdf_text(paste(dir_in, file, sep = ""))


# https://www.youtube.com/watch?v=bJH-S2iaxNo&ab_channel=DataCentricInc.
clean_tab <- function(x, last_sample){
  x <- map(x, ~ str_split(.x, "\\n") %>% unlist())
  x <- reduce(x, c)
  
  tab_start <- str_which(tolower(x), "library id")[1]
  tab_end <- str_which(x, last_sample)[1]
  
  tab_end <- tab_end[min(which(tab_end > tab_start))]
  
  tab <- x[(tab_start):(tab_end)]
  tab <- str_replace_all(tab, "\\s{2,}", "|")
  tab <- tab[-(1:4)]

  tab_df <- as.data.frame(do.call(rbind, strsplit(tab, "|", fixed = TRUE)))[,-1]
  
  colnames(tab_df) <- c("LibraryID", "GAPID", "SampleID", "Sample", "Organism", "Group", "QUBIT", "RIN")
  
  res <- tab_df %>% select(Sample, Organism, RIN)
  
  return(res)
}

clean_tab(file, "AC-58_L16")


################################################################################
#                             SAMPLE INFORMATION                           
################################################################################


# Number of samples
n <- length(samples)

# Sample names
sample_names <- str_sort(samples, numeric = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample information
sample_info <- data.frame(
  Sample = sample_names,
  Time = c(rep("4", 4), rep("24", 4), rep("48", 4), rep("Control", 4)),
  RIN = c("10","10","10","10","10","10","10","10","10","10","10","9.9","10","10","10","10"),
  Organism = rep(x = "Homo sapiens", n),
  Ind = c("E1","E2","E3","E4","E1","E2","E3","E4","E4","E1","E2","E3","E2","E3","E4","E1"))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------



#######################################################################
#                             SAVE DATA                         
#######################################################################

# Save sample information file
write.csv(sample_info, file = paste(dir_out, "/Sample_info.csv", sep = ""), row.names = FALSE)
