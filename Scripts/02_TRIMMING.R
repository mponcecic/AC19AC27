################################################################################
#                             TRIMMING FASTQ  
################################################################################
 

# Summary
# ---------

# Trimmed the fastq file prior to mapping with STAR. This step is optional. 
# In this script, we can perform a standard trimming and a more specific trimming 
# from the projects sequenced with the SMARTer Stranded Total RNA-Seq kit v2 - 
# Pico input Mammalian when the insert is size smaller than 150 bps. Common steps 
# are: 
#
# The first step is to remove the adapters from both reads simultaneously which 
# is performed using -a and -A for the first and second read, followed by the 
# corresponding adapters. We use the Illumina universal adapter which can be seen 
# in https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm
#
# Adapter R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Adapter R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
#
# The second step correspond to the filtering based on the minimum read size 
# (-m=30) and the quality (-q 10).
# 
# This are the steps common from both libraries. The unique step in SMARTer-Pico 
# is the trimming the three random nucleotides present in reads is performed using 
# -u and -U 3 after removing the adapter. 
# 
# The differences in the steps are a consequence of how the reads looks like: 
# 
# Library - Truseq
# 
# Read 1
# ---------------INSERT--------------- 
# ---------------INSERT---------------Adapter  
#
# Read 2 
# ---------------INSERT--------------- 
# ---------------INSERT---------------Adapter
# 
# Library -  SMARTer Pico 
# 
# Read 1
# ---------------INSERT--------------- 
# ---------------INSERT---------------XXX
# ---------------INSERT---------------XXX-Adapter  INSERT < 150 nucleotides
#
# Read 2 
# XXX---------------INSERT--------------- 
# XXX---------------INSERT---------------Adapter
# 
# XXX 3 random nucleotides 


# Folder
# Input: W:/DATA_shared/Sequencing_name/
# Output: Project_folder/01_TRIMMED



################################################################################
#                         CUTADAPT VERSION
################################################################################

# Conda version is the latest in March 2022
# Conda 4.12.0
#
# Condaconda is located in
# /opt/ohpc/pub/apps/anaconda3/cic-env
# 
# Cutadapt version 4.3 with Python 3.10.10 in environment cutadaptenv
# 
# Link: https://cutadapt.readthedocs.io/en/stable/index.html



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
project_name <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "/vols/GPArkaitz_bigdata/user/"
# path <- "W:/user/"


# Date of the log file 0_Sample_info_XXX.log
logdate <- ""


### Process information ###
partition <- "FAST"
time <- c("00:40:00")
memory <- c("1")
cpu <- 4


### Cutadapt parameters ###

#   -a/A Adapter for the read 1 and 2, respectively. 
#   -q Trim low-quality bases from 5' and/or 3' ends of each read before adapter 
#     removal 
#   -m Discard reads shorter then the length 
#   --pair-filter=any Default option in cutadapt; Any of the reads in paired-end 
#     read have to much a criterion
#   -j Number of CPU cores. Set 0 to autodetect the number of cores available
#   -u/U Remove LEN bases from each read. Applied before adapter trimming. In the 
#     case of using SMARTer Pico, this has to be applied after the adapter trimming


# Library preparation kit
# Options: normal, pico (special trimming) and new (add parameter)
seq_library <- "normal"

## Adapter R1
# Illumina universal adapter
a <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
## Adapter R2
# Illumina universal adapter
A <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Quality threshold per base
q <- 10
# Length threshold
# The read must present the at least this minimum length
# m <- 20
m <- 30

# Remove 3 nucleotides in the reads 
u <- "3"
U <- "3"

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Output directories 
dir_out <- paste(path, project_name, sep = "")

# Load log file 
logfile <- read.table(paste(dir_outfiles, "/log/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory
dir_infiles <- paste(logfile$filedirRocky, "FASTQs", sep = "")

# Create output directory
dir.create(file.path(dir_out,"02_TRIMMED"))
dir_outfiles <- paste(dir_out,"/02_TRIMMED",sep='')
setwd(dir_outfiles)


### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path=dir_infiles, pattern = "_1.fastq.gz")
# Filter sample name
samples <- gsub("_1.fastq.gz", "", samples)

# Minnimum length parameter
min_length <- paste("-m=", m, sep = "")



###############################################################
#                      TRIMMING
############################################################### 


### Generate PBS ###
for (i in 1:length(samples)) {
  # Job name
  job_name <- paste(samples[i],"_Trimmed",sep='');

  # Input file
  input1 <- paste(dir_infiles, "/", samples[i],"_1.fastq.gz", sep="")
  input2 <- paste(dir_infiles, "/", samples[i],"_2.fastq.gz", sep="")
  
  # Output file: Adapter trimmed, q<10 and  m>30
  output1 <- paste(dir_outfiles, "/", samples[i],"_1_trmd.fastq.gz", sep="")
  output2 <- paste(dir_outfiles, "/", samples[i],"_2_trmd.fastq.gz", sep="")
  
  
  if(seq_library == "normal"){
    # Cutadapt command
    command <- paste("cutadapt -a", a,"-A", A,"-j 0 -q", q, min_length, " --pair-filter=any -o", output1,"-p", output2, input1, input2, sep=" ")
    
    # SBATCH File
    filename <- paste(job_name,".sh",sep='');
    cat(
      c("#!/bin/sh"),
      c("#SBATCH  --export=ALL"),
      paste("#SBATCH --job-name=",job_name,sep=''),
      paste("#SBATCH --partition=",partition,sep=''),
      paste("#SBATCH --cpus-per-task=",cpu,sep=''),
      paste("#SBATCH --time=",time,sep=''),
      paste("#SBATCH --mem=",memory,"GB",sep=''),
      paste("#SBATCH -o ",job_name,".out",sep=''),
      paste("#SBATCH -e ",job_name,".err",sep=''),
      c(paste("cd ", dir_infiles)),
      c("source /opt/ohpc/pub/apps/anaconda3/cic-env"),
      c("conda config --set channel_priority strict"),
      c("conda activate cutadaptenv"),
      c(command),
      file=filename,sep = "\n",append=F)
    
  } else {
    # Intermediate file to apply the u/U parameter after removing the adapter
    outputm1 <- paste(dir_outfiles, "/", samples[i],"_1_m1.fastq.gz", sep="")
    outputm2 <- paste(dir_outfiles, "/", samples[i],"_2_m1.fastq.gz", sep="")
  
    
    # Cutadapt command
    command <- paste("cutadapt -a", a,"-A", A,"--pair-filter=any -o", outputm1,"-p", outputm2, input1, input2, sep = " ")
    command2 <- paste("cutadapt -u", u,"-U", U,"-j 0 -q", q, min_length, " --pair-filter=any -o", output1,"-p", output2, outputm1, outputm2, sep = " ")
    
    # SBATCH File
    filename <- paste(job_name,".sh",sep='');
    cat(
      c("#!/bin/sh"),
      c("#SBATCH  --export=ALL"),
      paste("#SBATCH --job-name=",job_name,sep=''),
      paste("#SBATCH --partition=",partition,sep=''),
      paste("#SBATCH --cpus-per-task=",cpu,sep=''),
      paste("#SBATCH --time=",time,sep=''),
      paste("#SBATCH --mem=",memory,"GB",sep=''),
      paste("#SBATCH -o ",job_name,".out",sep=''),
      paste("#SBATCH -e ",job_name,".err",sep=''),
      c(paste("cd ", dir_infiles)),
      c("source /opt/ohpc/pub/apps/anaconda3/cic-env"),
      c("conda activate cutadaptenv"),
      c(command),
      c(command2),
      file=filename,sep = "\n",append=F)
  } 
  
  
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
  }


#######################################################################
#                            LOG FILE                        
#######################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project_name <- project_name
log_data$outputdir <- dir_outfiles
log_data$inputdir <- dir_infiles
log_data$Library <- seq_library
log_data$AdaptR1 <- a 
log_data$AdaptR2 <- A
log_data$q <- q
log_data$m <- m
log_data$nucR1 <- u
log_data$nucR2 <- U


write.table(as.data.frame(log_data), paste(dir_out, "/log/2_Trimming_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")


q()


