
################################################################################
#                       STAR CONTROL SUMMARY                           
################################################################################

# Summary
# ---------

# After runnig STAR in Rocky, we want to know which percentage of the reads were 
# uniquely mapped to the reference genome. This parameter will give information 
# of how good the mapping was in among all the sample.  


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

# Input directory
input_dir <- paste(path, project_name,"/04_STAR", sep = "")

# Output directory
dir.create(file.path(paste(path, project_name, sep = ""),"VERIFICATION"))
output_dir <- paste(path, project_name,"/VERIFICATION", sep ="")


# Set working directory
setwd(input_dir)

# List fastq.zip files
samples <- gsub(".out", "", list.files(pattern = ".out"))
# Filter sample name 
samples_names <- gsub("_STAR", "", samples)

# Relevant information from the .out file
pattern <- c("Number of input reads", "Uniquely mapped reads %")

# Final matrix
def_file <- data.frame()


################################################################################
#                            PROCESS                           
################################################################################


for (i in 1:length(samples_names)){
  # Sample name
  setwd(as.character(samples[i]))
  
  file_name <- "Log.final.out"
  
  # Read .out file
  out_data <- readLines(file_name)
  
  # Select presenting the information of interest
  out_data <- data.frame(x = out_data[grepl(paste(pattern, collapse = " | "), out_data)])
  
  # Split the strings into two columns using " |\t" as separator
  #   First column. Variable of interest corresponding to pattern
  #   Second column. Value of the variable
  out_file <- data.frame(do.call("rbind",strsplit(as.character(out_data[,1]), "|\t", fixed = TRUE)))
  
  ## Key step for the final output
  # Transpose matrix
  out_data1 <- setNames(data.frame(t(out_file[,-1])), out_file[,1])
  
  # Aggregate all the information together of the corresponding file 
  out_data1 <- cbind(Sample=paste(samples_names[i]), out_data1)
  
  # Aggregate all the information from every file
  def_file <- rbind(def_file, out_data1)
  
  # Change directory to STAR_AC35
  setwd("../")
}



#Load trimming check
if(timming == TRUE){trim_file <- read.csv(file = paste(output_dir, "/", project_name,"_Trimming_check.csv", sep = ""), header = TRUE)
colnames(trim_file) <- gsub("\\.", " ", colnames(trim_file))
colnames(trim_file) <- gsub("   ", " < ", colnames(trim_file))

# Save together the results of the trimming and the alignment
final_check <- merge(trim_file,def_file, by = "Sample")
write.csv(final_check, file = paste(output_dir,"/", project_name,"_Trimmed_STAR_check.csv", sep = ""), row.names=FALSE)

} else {final_check <- def_file}





################################################################################
#                   SAVED DATA: Performance               
################################################################################


write.csv(def_file, file = paste(output_dir ,"/", project_name,"_STAR_check.csv", sep = ""), row.names=FALSE)











