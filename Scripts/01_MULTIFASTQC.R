################################################################################
#                         MULTIFASTQC REPORT                
################################################################################

# Summary
#---------

# Summary report with the fastqc or quality control from the fastq files.

# Folder
# Input: Project_folder/01_FASTQC
# Output: Project_folder/01_FASTQC


################################################################################
#                             FASTQCR VERSION
################################################################################

# Fastqcr version is the latest in Febrary 2023
# Fastqcr v0.1.3
# Link: https://cran.r-project.org/web/packages/fastqcr/



################################################################################
#                   SET DIRECTORIES AND LOAD FILES             
################################################################################


#-------------------------------------------------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project <- "XXX"


# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Output and input directories SSH 
dir_infiles <- paste(path, project, "/01_FASTQC", sep = "" )
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


qc_report(qc.path = dir_infiles, result.file = paste(dir_outfiles, "/",  project, "_MultiReport", sep = ""), 
          interpret = TRUE)

