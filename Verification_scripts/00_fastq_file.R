################################################################################
#                        FASTQ FILES SUMMARY                           
################################################################################

# Summary
# ---------

# File size corresponding to each sample read (fastq file)


################################################################################
#                         LOAD FILES AND DATA                           
################################################################################


# Input directory
input_dir <- "W:/DATA_shared/AC-58_TotalRNAseq/FASTQs/"
# Output directory
output_dir <- "W:/mponce/AC58/VERIFICATION"

# Set working directory
setwd(input_dir)

# List fastq.gz files
samples <- list.files(pattern = ".fastq.gz")
# Filter sample name 
samples_names <- gsub(".fastq.gz", "", samples)

# Final matrix
def_file <- data.frame()


for (i in 1:length(samples_names)){
  # Sample name
  file_name <- samples_names[i]
  # Sample *.fastq.gz size file
  file_size <- file.size(paste(getwd(),"/", samples[i], sep=""))
  
  # Aggregate all the information together of the corresponding file 
  out_data <- data.frame(cbind(Sample=file_name,Size=file_size))
  
  # Aggregate all the information from every file
  # Output 
  # 
  #   Sample    Size 
  #   CAS1_1   Value 1          
  #     .      Value 1         
  #     .      Value 1         
  #   CAS24_2  Value 1
  def_file <- rbind(def_file, out_data)
}


################################################################################
#                   SAVED DATA: Performance               
################################################################################

# Export to *.csv file on the personal folder
write.csv(def_file, file = paste(output_dir,"/00_Fastq.csv", sep = ""), row.names = FALSE)

