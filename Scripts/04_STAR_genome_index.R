#################################################################################
#                         STEP 2.1: STAR GENOME INDEX
#################################################################################

# Summary
#---------

# The reads from the AC58 project has been trimmed and followed by a quality 
# control check. The next step is to map sequence reads to a known transcriptome 
# or annotated genome, transforming the reads into genomic coordinates.
#
# The reads will be mapped using STAR which is a genome reference based aligner. 
# To run STAR two steps must be followed:
#  - Step 1. Build STAR genome index (00_STAR_genome_index.sh)
#  - Step 2. Mapping reads with STAR is performed in this script
#
# In this script, the first step is performed, building the genome index. The 
# reference genome for the read length 101 based on the genome assembly hg38 of 
# 2020 which is compatible with GRCh38, and the corresponding genome annotation 
# for the 2-pass mapping. The 2-pass mapping aims to map novel splice junctions.

# ISSUE 
#----------
# 
# The reference genome index in:
#               /vols/GPArkaitz_bigdata/DATA_shared/Genomes/Indexes/Human_101
# is not valid for the STAR version. So we run this script to generate the new 
# genome index. This step can be avoid for the next DEGs projects which fullfill 
# the characteristics of this project.



###############################################################
#               Experimental Information
###############################################################

# Project: AC58

# The raw data analyzed is contained in /vols/GPArkaitz_bigdata/DATA_shared/AC58
# Human cell line
# Total bulk RNAseq paired-end
# The inserts size 40 - 300 bp
# Library average size 327 bp
# Read length 101
# Library kit TruSeq Stranded Total RNA



###############################################################
#                       STAR Version
###############################################################

# R version is the latest in March 2023
# R version 4.3.0 (2023-04-21 ucrt)

# STAR version 2.7.10b


###############################################################
#                       STAR Manuals
###############################################################

## Manuals used 

# STAR manual 2.7.10b
# October 20, 2022
# Link: https://github.com/alexdobin/STAR


# STAR manual 2.7.0a
# January 23, 2019
# Link: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf




###############################################################
# How to run in ROCKY
# 
# 1. Input folder is inside my user folder in BigData
#
# 2. Output folder is inside my user folder in BigData
#
# 3. This script will generate all the "PBS" file, one per each 
# sample, which includes the commands to successfully run STAR. 
# The PBS are directly send to the job scheduler (SLURM) to run 
# the jobs in INDAR.
#
###############################################################



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
dir.create(file.path(dir_outfiles,"03_STAR"))
dir_outfiles <- paste(dir_outfiles,"/03_STAR",sep='')
setwd(dir_outfiles)
dir.create(file.path(dir_outfiles,"Human_101"))
dir_outfiles <- paste(dir_outfiles,"/Human_101",sep='')
# Set directory
setwd(dir_outfiles)

fasta <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes/refdata-gex-GRCh38-2020-A/fasta/genome.fa" 
gtf <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
read <- "100" 


### Process information ###
partition <- "FAST"
time <- c("03:00:00")
memory <- c("40")
cpu <- 8
ram <- 30000000000 



###############################################################
#                           STAR
###############################################################

### Generate PBS ###

job_name <-  paste("Genome_index", read, project_name, sep = "_")

# STAR path
path <- "source /opt/ohpc/pub/apps/star/STAR-2.7.10a/cic-env"

# Command
command <- paste("STAR --runThreadN", cpu,    
                 "--runMode genomeGenerate",
                 "--genomeDir", dir_outfiles,                     
                 "--genomeFastaFiles", fasta,              
                 "--sjdbGTFfile", gtf,
                 "--sjdbOverhang", read, 
                 sep=" ") 

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
  # c("#SBATCH --exclude=gn[10-13]"),
  c("#SBATCH --mail-type=FAIL"),
  c("#SBATCH --mail-user=mponce@cicbiogune.es"),
  "\n\n",
  c(path), 
  "\n",
  c(paste("cd ", dir_infiles)), 
  "\n",
  c(command),
  file=filename,sep = "\n",append=F)

# Run PBS
system(paste("sbatch",filename,sep=' '));
Sys.sleep(3)

q()
