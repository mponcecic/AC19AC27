################################################################################
#                           FASTQC TRIMMED FASTQ           
################################################################################

# Summary
# ---------

# Quality control of the trimmed fastq files

# Folder
# Input: Project_folder/01_TRIMMED
# Output: Project_folder/02_FASTQC


################################################################################
#                             FASTQC VERSION
################################################################################

# FastQC version is the latest in March 2023: 
# FastQC v0.12.1
#
# FastQC is located in: 
# /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/fastqc 

# Link: https://github.com/s-andrews/FastQC


################################################################################
#                             PIPELINE
################################################################################


### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
# Load R libraries 
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/DEG_Rocky")


#-------------------------------------------------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project_name <- "AC58"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "/vols/GPArkaitz_bigdata/mponce/"
# path <- "W:/mponce/"

### Process information ###
cluster <- "FAST"
walltime <- c("00:10:00") 
cpu <- 1
memory <- c("3")
#-------------------------------------------------------------------------------------------------------------------------------------------------------


# Input directories 
dir_infiles <- paste(path, project_name, "/02_TRIMMED", sep = "" )

# Create output directory
dir_outfiles <- paste(path, project_name, sep = "")
dir.create(file.path(dir_outfiles,"03_FASTQC"))
dir_outfiles <- paste(dir_outfiles,"/03_FASTQC",sep='')
setwd(dir_outfiles)

### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path=dir_infiles, pattern = "trmd.fastq.gz")
# Filter sample name 
samples_names <- gsub("_trmd.fastq.gz", "", samples)

### Generate PBS ###
for (i in 1:length(samples)) {
  # Job name
  name.job <- paste(samples_names[i],"_FastQC",sep='');
  # Command
  command <- paste("/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/FastQC/fastqc", samples[i], "-o ", dir_outfiles)
  # SBATCH File
  filename <- paste(name.job,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    paste("#SBATCH --job-name=",name.job,sep=''),
    paste("#SBATCH --partition=",cluster,sep=''),
    paste("#SBATCH --cpus-per-task=",cpu,sep=''),
    paste("#SBATCH --time=",walltime,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",name.job,".out",sep=''),
    paste("#SBATCH -e ",name.job,".err",sep=''),
    c(paste("cd ", dir_infiles)),
    c(command),
    file=filename,sep = "\n",append=F)
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
  }


q()

