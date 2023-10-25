##################################################
#              MULTIFASTQC REPORT                #
##################################################

# Summary
#---------

# Create a report with the result of the quality 
# control performed in each fastq trimmed file per 
# sample of AC62 project.

# RUN LOCALLY IN RSTUDIO


####################################################
#                   LOAD LIBRARIES             
####################################################


library(fastqcr)


####################################################
#           SET DIRECTORIES AND LOAD FILES             
####################################################


### General Project ###
# Project Name 
project_name <- "AC58"

# Output and input directories SSH 
dir_infiles <- paste("W:/mponce/", project_name, "/05_DEGs/02_FASTQC", sep = "" )
dir_outfiles <- dir_infiles
setwd(dir_outfiles)


### FASTQCs ###
# List fastqc.gz files
samples <- list.files(path = dir_infiles, pattern = "_fastqc.zip")
# Filter sample name 
samples_names <- gsub("_fastqc.zip", "", samples)


####################################################
#               MULTIFASTQC REPORT               
####################################################

qc_report(qc.path = dir_infiles, result.file = paste(dir_outfiles, "/",  project, "_MultiReport", sep=""), 
          interpret = TRUE)

