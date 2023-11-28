
################################################################################
#                               METADATA                          
################################################################################


# Summary
# --------
#
# Generate a sample information file also refereed as metadata. The sample 
# information file is used since the beginning of the analysis to generate 
# different scripts and in the differential expression analysis which is based 
# in the project information. Allowing the execution of the scripts with minor 
# adjustments.
# 
# This file must be MANUALLY adjusted for the different projects. However, part 
# of the information is rescued from the pdf file with the library preparation 
# information (MUST HAVE THE PATTERN Library_Preparation). You should look for 
# the table with the information and save the first and last sample from the 
# table.
# 
# Content
#   - Named Sample_info.csv
#   - Must have at least, the following variables:
#     - Sample: Variable with the samples names. MUST NAME Sample
#     - Treatment: Specific treatment per sample. Must be a factor variable, can 
#       be named as you prefer but must be consistent for the following analysis


# Folder
# Input: W:/DATA_shared/Sequencing_data_folder{ACXX}/
# Output: W:/PersonalFolder/Project folder


################################################################################
#                                 LOAD DATA                         
################################################################################


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "XXX"


# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"


# Local directory with Git folders 
local_dir <- "C:/Users/mponce/CIC bioGUNE/Arkaitz group - Maria Ponce/Projects/00_DEG_Reference"


# Input directory
# Must be the folder where the library preparation pdf is found and folder 
# with the fastq are included 
files_rocky <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-XX_TotalRNAseq/"
files <- "W:/DATA_shared/AC-XX_TotalRNAseq/"
# Select input directory
dir_in <- files
  
# First and last sample found in the Library Preparation pdf table which is 
# usually in page 5.
# This is a key step to generate a data frame based on the table from the pdf
first_sample <- "AC-XX_XX"
last_sample <- "AC-XX_XX"


# Condition
trt <- "XXX"

# Contrast level order 
lvl_order <- c("A", "B", "C", "D")

# Organism
# org <- "Mouse"
org <- "Human"

# Read length
read <- 151


# Contrast
# If contrast set NULL the different contrasts will be created based on the 
# lvl_order vector. Nonetheless, you can add the contrast manually following the 
# 
# DESeq2 formula:
# contrast <- list(c("Condition", "Experimental level", "Reference level"), 
#                  c("Condition", "Experimental level", "Reference level"))
contrast <- NULL
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


################################################################################
#                         PROJECT CONFIGURATION
################################################################################



## Create folders BigData folders
# Create project folder in BigData
dir.create(file.path(path, project))
dir_out <- paste(path, project, sep = "")
setwd(dir_out)

# Create utils folder in BigData
dir.create(file.path(dir_out, "utils"))
dir_utils <- paste(dir_out, "/utils", sep = "")

# Create log file folder
dir.create(file.path(dir_out, "log"))
dir_log <- paste(dir_out, "/log", sep = "")


## Utils folder transfer files 
# utils folder in local directory
local_dir <- paste(local_dir, "/utils", sep = "")

# List file names from utils
file_list <- list.files(local_dir)

# Copy the files
file.copy(from = file.path(local_dir, file_list), to = dir_utils)



################################################################################
#                     LOAD LIBRARIES AND FUNCTIONS                           
################################################################################


# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))
# Load functions 
source(paste(path, project, "/utils/function_pdf_to_tab.R", sep = ""))
source(paste(path, project, "/utils/functions_degs.R", sep = ""))



################################################################################
#                         SAMPLE INFORMATION FILE
################################################################################


# List of samples names
samples <- unique(gsub("_1.fastq.gz","", list.files(path = paste(dir_in, "FASTQs/", sep = ""), pattern = "_1.fastq.gz")))

# PDF file with the library preparation report
# This file should be in the project DATA_shared folder 
file <- pdf_text(paste(dir_in, list.files(path = dir_in, pattern = "Library_Preparation"), sep = ""))


################################################################################
#                             SAMPLE INFORMATION                           
################################################################################


# Run the function
data <- pdf_to_tab(file, first_sample, last_sample)
print(head(data))

# Number of samples
n <- length(samples)

# Verification
ifelse(nrow(data) == n, print("The number of samples from the data table in the pdf match the samples in the folder"), paste("ERROR: Samples do not match", "CHANGE THE NUMBER OF ROWS REMOVED IN THE FUNCTION", sep = "\n."))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Select manually the columns because theorder is not the same in all the pdfs
# Usually select 
#  - Sample names compatible with the file names
#  - Organism
#  - RIN values
data <- data %>% select(V5, V6, V9) %>% arrange(desc(V5))
colnames(data) <- c("Sample", "Organism", "RIN")

# Add the condition variable
# Be careful adding the comparison levels
data[trt] <- factor(c(rep("Control", 4), rep("4", 4), rep("48", 4), rep("24", 4)), levels = lvl_order, ordered = TRUE)

# Optional
# Some sequencing projects migth include information from different tissues, projects, ...
# You can add this information as a new variable or generate two separet csv files
data <- data %>% mutate(Ind = sub(".*E", "E", data$Sample)) %>% arrange(Time)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Save data as sample information 
sample_info <- data


################################################################################
#                           CREATE CONTRAST                         
################################################################################


if(is.null(contrast)== TRUE){contrast <- create_contrast(trt, lvl_order)}


#######################################################################
#                             SAVE DATA                         
#######################################################################


# Save sample information file
write.csv(sample_info, file = paste("Sample_info.csv", sep = ""), row.names = FALSE)



#######################################################################
#                            LOG FILE                        
#######################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project_name <- project
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_order, collapse =",")
log_data$path <- dir_out
log_data$filedir <- dir_in
log_data$filedirRocky <- files_rocky
log_data$Organism <- org
log_data$read <- read
log_data$contrasts <- paste(unlist(contrast), collapse = ",")

write.table(as.data.frame(log_data), paste(dir_log, "/0_Sample_info_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



