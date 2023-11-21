
################################################################################
#                       TRIMMING CONTROL SUMMARY                           
################################################################################

# Summary
# --------

# After runnig cutadapt in Rocky, we want to verify all the jobs were successful 
# before we keep going with the analysis.  


################################################################################
#                         LOAD FILES AND DATA                           
################################################################################


# -----------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project_name <- "AC58"

# File path
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"
# -----------------------------------------------------------------------------------------------------------------

# Output and input directories SSH 
input_dir <- paste(path, project_name, "/02_TRIMMED", sep = "" )
output_dir <- paste(path, project_name, "/VERIFICATION", sep = "")

# Set working directory
setwd(input_dir)

# List fastq.zip files
samples <- list.files(pattern = "_trmd.fastq.gz")
# Filter sample name 
samples_names <- gsub("_trmd.fastq.gz", "", samples)

# Relevant information from the .out file
pattern <- c("State","CPU Efficiency", "Job Wall-clock", "Memory Efficiency")

# Final matrix
def_file <- data.frame()


for (i in 1:length(samples_names)){
  # Sample name
  file_name <- paste(samples_names[i], "_Trimmed.out", sep="")
  # Sample *.fastqc.zip size file
  file_size <- file.size(paste(getwd(),"/", samples[i], sep=""))
  
  # Read .out file
  # out_file <- read.fwf(file_name)
  # out_data <- data.frame(read.delim(file = file_name, header = FALSE))
  out_data <- readLines(file_name)
  
  # Select presenting the information of interest
  out_data <- data.frame(x = out_data[grepl(paste(pattern, collapse = "|"), out_data)])
  
  # Split the strings into two columns using ": " as separator
  #   First column. Variable of interest corresponding to pattern
  #   Second column. Value of the variable
  out_file <- data.frame(do.call("rbind",strsplit(as.character(out_data[,1]), ": ", fixed = TRUE)))
  
  ## Key step for the final output
  # Transpose matrix
  # Output
  #   CPU Efficiency  Job Wall-clock  Memory Efficiency
  #      Value 1          Value 2           Value 3
  out_data <- setNames(data.frame(t(out_file[,-1])), out_file[,1])
  
  # Aggregate all the information together of the corresponding file 
  out_data <- cbind(Sample=paste(samples_names[i]), Size_R1=file_size, out_data)
  
  # Aggregate all the information from every file
  # Output 
  # 
  #   Sample  CPU Efficiency  Job Wall-clock  Memory Efficiency
  #   CAS1_1   Value 1          Value 2           Value 3
  #     .      Value 1          Value 2           Value 3
  #     .      Value 1          Value 2           Value 3
  #   CAS24_2   Value 1          Value 2           Value 3
  def_file <- rbind(def_file, out_data)
}


################################################################################
#                   SAVED DATA: Performance               
################################################################################


write.csv(def_file,file = paste(output_dir ,"/02_Trimmed_Fastq_file.csv", sep = ""), row.names=FALSE)
