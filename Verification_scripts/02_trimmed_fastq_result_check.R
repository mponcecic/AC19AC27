
################################################################################
#                             TRIMMED FASTQ SUMMARY                            
################################################################################

# Summary
# ---------

# After runnig cutadapt in Rocky, we want to study how the fastq file and the 
# reads were affected. We select the total number of reads processed, the reads 
# with adapter content, reads too short, quality trimmed and the total reads 
# written.


################################################################################
#                         LOAD FILES AND DATA                           
################################################################################


# -----------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project_name <- "XXX"

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
samples <- list.files(pattern = "_1_trmd.fastq.gz")
# Filter sample name 
samples_names <- gsub("_1_trmd.fastq.gz", "", samples)

# Relevant information from the .out file
pattern <- c("Total read pairs processed", "Read 1 with adapter", "Read 2 with adapter", 
             "Pairs that were too short", "Pairs written (passing filters)", 
             "Quality-trimmed", "Total written (filtered)")


# -- Adaptador --
# Total read pairs processed
# Read 1 with adapter
# Read 2 with adapter
# 
# -- seq menores 30 --
# Pairs that were too short
# Pairs written (passing filters)
# 
# -- seq calidad mayor 10 --
# Quality-trimmed
# Total written (filtered)


# Final matrix
def_file <- data.frame()


################################################################################
#                            PROCESS                           
################################################################################


for (i in 1:length(samples_names)){
  # Sample name
  file_name <- paste(samples_names[i], "_Trimmed.out", sep="")
  
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
  out_data1 <- setNames(data.frame(t(out_file[,-1])), out_file[,1])
  # out_data2 <- out_data1[,-c(3,7)]
  colnames(out_data1) <- c("Total reads pairs processed", "Read 1 with adapter Truseq", "Read 2 with adapter Truseq", "Pairs < 30bps", "Quality-trimmed < 10")
  
  # Aggregate all the information together of the corresponding file 
  out_data3 <- cbind(Sample=paste(samples_names[i]), out_data1)
  
  # Aggregate all the information from every file
  # Output 
  # 
  #   Sample  Total reads pairs processed  Read 1 with adapter Truseq     ...
  #   CAS1_1            Value 1                     Value 2               ...     
  #     .               Value 1                     Value 2               ...     
  #     .               Value 1                     Value 2               ...     
  #   CAS24_2           Value 1                     Value 2               ...     
  def_file <- rbind(def_file, out_data3)
}

################################################################################
#                   SAVED DATA: Performance               
################################################################################


write.csv(def_file,file = paste(output_dir ,"/", project, "_Trimming_check.csv", sep = ""), row.names=FALSE)
