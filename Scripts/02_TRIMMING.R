###############################################################
#                   STEP 1: TRIMMING 
###############################################################
 
# Conda version is the latest in March 2022
# Conda 4.12.0
#
# Condaconda is located in
# /opt/ohpc/pub/apps/anaconda3/cic-env
# 
# Cutadapt version 4.3 with Python 3.10.10 in environment cutadaptenv


###############################################################
#               Experimental Information
###############################################################
#
# Project: AC35
#
# The raw data analyzed is contained in /vols/GPArkaitz_bigdata/DATA_shared/AC35  
# This data is the RNAseq from 24 mice
# The inserts size 40 - 300 bp
# Library average size 363 bp
# Paired-end reads
# Library kit TruSeq Stranded mRNA


###############################################################
#                           Problems
###############################################################
# 
# Problem 1
#-----------
#
# The sequencing data was produced with SMARTer Stranded Total RNA-Seq kit v2 - Pico input 
# Mammalian which used the SMARTer Pico v2 adapter. As a result, the first three nucleotides of 
# second read of paired-end reads are random nucleotides that must be trimmed to avoid multimapping.  
# 
#
# Problem 2
#-----------
#
# According to our data, the insert size can change from 40 to  300 pbs and the adapter size 
# is 139 pbs which means that when the  insert size is smaller than 150 pbs, the reads present 
# adapters content. In addition, the read length is 2x150 bp sequencing what cause an unusual 
# output when the inserts are smaller than 150 pbs. 
#
# Causing the transcription of the 3 random nucleotides on read 1, as well as, a partial 
# transcription of the read 2 adapter.


###############################################################
#                      How do the reads look
###############################################################
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
# XXX 3 nucleotides may be remove regardless of the quality


###############################################################
#                      Solution
############################################################### 
#
# 1. Adapters trimming
#   i. Cut the adapter corresponding to read 2 in the end (3') of the read 1***
#   ii. Cut the adapter corresponding to read 1 in the beginning (5') of the read 2
#
# 2. Three random nucleotides must be trimmed prior to mapping which correspond to the SMARTer 
# Stranded Pico v2 libraries. 
# 
# 3. Filter reads based on the quality (q > 10) and the read size (read > 30bp) to avoid
# poor alingment of the reads and false positive multimapping. The tendency in the quality
# of the reads tend to decrease at the end of the reads 
#
# *** Insert size smaller than 150 bps
#
#
# First step
#------------
#
# The first step is to remove the adapters from both reads simultaniously which is performed
# using -a and -A for the first and second read, followed by the corresponding adapters.
# The adapter can be seen in https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm
#
# Adapter R1: 
#           AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#
# Adapter R2:
#           AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
#
#
# Second step
#--------------
#
# The second step correspond to the trimming of the three nucleotides, as well as, the 
# filtering based on the minimum read size (-m=30) and the quality (-q 10).
# Trimming the three random nucleotides present in reads is performed using -u and -U 3. 


### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
# Load R libraries 
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_R")

### General Project ###

### General Project ###
# Output and input directories SSH 
dir_outfiles <- "/vols/GPArkaitz_bigdata/mponce/AC58" 
dir_infiles <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-58_TotalRNAseq/FASTQs"


# Create output directory
dir.create(file.path(dir_outfiles,"05_DEGs"))
dir_outfiles <- paste(dir_outfiles,"/05_DEGs",sep='')
dir.create(file.path(dir_outfiles,"01_TRIMMED"))
dir_outfiles <- paste(dir_outfiles,"/01_TRIMMED",sep='')
setwd(dir_outfiles)

### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path=dir_infiles, pattern = "_1.fastq.gz")
# Filter sample name
samples <- gsub("_1.fastq.gz", "", samples)


### Process information ###
partition <- "FAST"
time <- c("00:40:00")
memory <- c("1")
cpu <- 4


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

  # Command
  command <- paste("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0  -q 10 -m=30 --pair-filter=any -o", output1,"-p", output2, input1, input2, sep=" ")
  
  # SBATCH File
  filename <- paste(job_name,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    paste("#SBATCH --job-name=",job_name,sep=''),
    paste("#SBATCH --partition=",partition,sep=''),
    paste("#SBATCH --cpus-per-task=",cpu,sep=''),
    paste("#SBATCH --time=",time,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",job_name,".out",sep=''),
    paste("#SBATCH -e ",job_name,".err",sep=''),
    c("#SBATCH --mail-type=FAIL"),
    c("#SBATCH --mail-user=mponce@cicbiogune.es"),
    c(paste("cd ", dir_infiles)),
    c("source /opt/ohpc/pub/apps/anaconda3/cic-env"), 
    c("conda config --add channels defaults"),
    c("conda config --add channels bioconda"),
    c("conda config --add channels conda-forge"),
    c("conda config --set channel_priority strict"),
    #c("conda create -n cutadaptenv cutadapt"),       Download cutadapt environment
    c("conda activate cutadaptenv"),
    c(command),
    file=filename,sep = "\n",append=F)
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
  }

  # Pruebas
  # system("sbatch CAS23_Trimmed.sh");
  # Sys.sleep(3)

q()


