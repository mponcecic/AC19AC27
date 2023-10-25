##################################################
#             INTERACTIVE JOB ROCKY              #
##################################################

# Resume 
# ------------------------------------------------
# PBS are generated and run in the cluster on this 
# file. The aim of each file is to performed the 
# Fastq Control Quality (fastqc) in the samples of
# AC35 project.
# 
# Notes:
# This files is based on the previous work of Saioa
# 0_Saioa_CellRanger_count_Rocky_15092022


# # In this case memory is FAST 
# # Ask for a node and wait until it allocates one to you
# salloc -N 1 -n 1 --mem=1G -t 00:30:00 --partition=FAST --job-name=interactive
# # Enter the node assigned to you (gno02 or other)
# ssh c02
# # Jobs queues
# squeue
# squeue -u mponce


### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
# Load R libraries 
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_R")


### General Project ###
# Project Name 
project_name <- "AC58"

# Output and input directories SSH 
dir_infiles <- paste("/vols/GPArkaitz_bigdata/mponce/", project_name, "/05_DEGs/01_TRIMMED", sep = "" )
dir_outfiles <- paste("/vols/GPArkaitz_bigdata/mponce/", project_name, "/05_DEGs", sep = "")

# Create output directory
dir.create(file.path(dir_outfiles,"02_FASTQC"))
dir_outfiles <- paste(dir_outfiles,"/02_FASTQC",sep='')
setwd(dir_outfiles)

### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path=dir_infiles, pattern = "trmd.fastq.gz")
# Filter sample name 
samples_names <- gsub("_trmd.fastq.gz", "", samples)


### Process information ###
cluster <- "FAST"
walltime <- c("00:10:00") 
cpu <- 1
memory <- c("3")


####################################################
#              FastQC info               
####################################################
# FastQC version is the latest in March 2023: 
# FastQC v0.12.1
#
# FastQC is located in: 
# /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/fastqc 



####################################################
# How to run in INDAR
# 
# 1. - input folder in DATA_shared
# 
# 2. - output folder inside my user folder in Bigdata 
# and the save in the correspondig FastQC folder for 
# the project
# 
# 3. - This script will generate a file listing all 
# commands (one per sample) to send it to the home-made 
# job scheduler for INDAR
####################################################


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

