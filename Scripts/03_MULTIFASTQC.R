################################################################################
#                         MULTIFASTQC REPORT                
################################################################################

# Summary
#---------

# Summary report with the fastqc or quality control from the fastq files.

# Folder
# Input: Project_folder/02_FASTQC
# Output: Project_folder/02_FASTQC


################################################################################
#                             FASTQCR VERSION
################################################################################

# Fastqcr version is the latest in Febrary 2023
# Fastqcr v0.1.3
# Link: https://cran.r-project.org/web/packages/fastqcr/



################################################################################
#                   SET DIRECTORIES AND LOAD FILES             
################################################################################


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



################################################################################
#                     LOAD LIBRARIES AND FUNCTIONS                           
################################################################################


# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))



################################################################################
#                       MULTIFASTQC REPORT               
################################################################################


qc_report(qc.path = dir_infiles, result.file = paste(dir_outfiles, "/",  project, "_MultiReport", sep=""), 
          interpret = TRUE)

