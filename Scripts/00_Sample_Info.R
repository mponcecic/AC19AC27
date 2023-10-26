
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
# information. 
# 
# Content
#   - Named Sample_info.csv
#   - Must have at least, the following variables:
#     - Sample: Variable with the samples names. MUST NAME Sample
#     - Treatment: Specific treatment per sample. Must be a factor variable, can 
#       be named as you prefer but must be consisting for the following analysis



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

# Input directory
# Must be the folder where a the library preparation pdf is found and folder 
# with the fastq are included 
dir_in <- "W:/DATA_shared/AC-58_TotalRNAseq/"

# Last sample found in the Library Preparation pdf
last_sample <- "AC-58_L16"

# Condition
trt <- "Time"
  
# Contrast level order 
lvl_order <- c("Control", "4", "24", "48")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Output directory
dir_out <- paste(path, project, sep = "")

# Set working directory
setwd(dir_out)

# List of samples names
samples <- unique(gsub("_1.fastq.gz","", list.files(path = paste(dir_in, "FASTQs/", sep = ""), pattern = "_1.fastq.gz")))

# PDF file with the library preparation report
# This file should be in the project DATA_shared folder 
file <- pdf_text(paste(dir_in, list.files(path = dir_in, pattern = "Library_Preparation"), sep = ""))


################################################################################
#                             FUNCTION                           
################################################################################


# From pdf to data frame
pdf_to_tab <- function(x, last_sample){
  
  # Description
  # 
  # Transform the table from the pdf into a data frame to select the variables of 
  # interest
  # This can not be completely automatized because the tables differ from projects
  
  # Create a character vector using \n as the row separator
  x <- map(x, ~ str_split(.x, "\\n") %>% unlist())
  x <- reduce(x, c)
  
  # Select the beginning of the table
  # Must be always the same
  if(is.na(str_which(x, "Kibrary ID")[1])){
    tab_start <- str_which(x, "GAP ID")[1]
  } else {tab_start <- str_which(x, "Library ID")[1]}
  
  # Select the row of the end
  tab_end <- str_which(x, last_sample)[1]
  
  # Select the rows of interest
  tab <- x[(tab_start):(tab_end)]
  
  # Change column separator form spaces to "|" and remove the first four rows
  tab <- str_replace_all(tab, "\\s{2,}", "|")
  tab <- tab[-(1:4)]
  
  # Create a data frame of the table using "|" as a separator 
  # Remove the first column beacuse it's empty
  tab_df <- as.data.frame(do.call(rbind, strsplit(tab, "|", fixed = TRUE)))[,-1] 
  
  return(tab_df)
}


################################################################################
#                             SAMPLE INFORMATION                           
################################################################################


# Run the function
data <- pdf_to_tab(file, last_sample)
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


#######################################################################
#                             SAVE DATA                         
#######################################################################

# Save sample information file
write.csv(sample_info, file = paste(dir_out, "/Sample_info.csv", sep = ""), row.names = FALSE)
