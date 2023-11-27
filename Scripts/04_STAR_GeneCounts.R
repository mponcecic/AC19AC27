################################################################################
#                             GENE COUNTS                         
################################################################################

# Summary
# ---------
# 
# After aligning the reads with STAR, we used the output ReadsPerGene.out.tab 
# which is manipulated to obtain the raw counts file used in the differential 
# expressed analysis. To accomplish this aim we first need to decide the reads
# directionality and which one is a reflect of the coding strand.  
# 
# In addition, we perform the counts per million normalization method on the 
# raw counts data based on all the samples. 


# Folder
# Input: Project_folder/04_STAR
# Output: Project_folder/04_STAR



################################################################################
#                       How to choose the stranded
################################################################################


## Gene Counts (ReadsPerGene.out.tab) 

# Column 1. Gene ID
# Column 2. counts for unstraned RNA-seq
# Column 3. counts for 1st read strand aligned with RNA
# Column 4. counts for 2nd read strand aligned with RNA 


# First step: Stranded or unstranded?
# ------------------------------------
# 
# Based on the name of the protocol, Illumina Truseq mRNA stranded, reads are 
# stranded. The second column which corresponds to the unstranded DNA is not an 
# option. 


# Second step: Stranded or unstranded?
# ------------------------------------
#
# According to library preparation manual, we choose the fourth column 
# (counts for 2nd read strand aligned with RNA). Because the first cDNA was 
# synthesized using SMARTScribe Reverse Transcriptase. In our case, the second 
# reads align with the coding DNA strand were the gen is found. 
# 
# In case, you don't know the strandedness of the protocol, you must select the 
# columns with more counts. The ratio should be 1:50 or even 1:100.
# 
# For a deeper understanding keep reading
# 
# The DNA is double stranded which means the the strand can be template/antisense 
# or coding/sense when being classified based on the gen. The coding strand is 
# read in the 5' -> 3' direction, meaning this is the direction of the gene.
# 
# The RNA polymerase synthesizes the mRNA in the  5' -> 3' direction, which means 
# that the mRNA presents the same information as the coding strand. Bear this idea 
# in mind because we align the reads towards a genome and not a transcriptome.
# 
# At this point, the Illumina kit takes part. A reverse transcriptase is used to 
# obtain the first cDNA, after adding a primer in the 3' end. The cDNA is read in 
# the 3' -> 5' direction which refers to the template/antisense DNA strand. Remember, 
# the reverse transcriptase added 3 nucleotides at the end of the synthesis. 
# 
# Afterwards, the 5' and 3' primers are added and the reads are synthesized. As 
# a result, the second read can be read in the  5' -> 3' direction and the first 
# read in the 3' -> 5' direction. 
# 
# In this case, the second reads is the one presenting the reading direction of 
# the mRNA and the first strand is form the template strand. Being this the reason 
# why it's chosen which is a consequence of reverse transcriptase. 
# 
# 
# Schema
#                     
# DNA (coding strand) ----------------- GEN -----------------------
#   5' -> 3'
#
# mRNA                ----------------- GEN -----------------------
#   5' -> 3'
# 
# first cDNA        XXX---------------INSERT---------------AdapterN3
#   3' -> 5'
#
# Read 1
#   3' -> 5'
# ---------------INSERT--------------- 
# ---------------INSERT---------------XXX
# 
# Read 2 
#   5' -> 3' 
# XXX---------------INSERT--------------- 
# XXX---------------INSERT---------------Adapter3
# 
# 
# *The 5' end is represented on the left side and the 3' end on the right
# 


# References
# -------------------
# https://chipster.csc.fi/manual/library-type-summary.html -- Directionality
# https://github.com/alexdobin/STAR/issues/842 -- When to choose column 3/4
# https://www.biostars.org/p/3423/ -- Directionality
# https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html#counting-reads-per-genes



################################################################################
#                             LOAD DATA                           
################################################################################


#------------------------------------------------------------------------------------------------------------------------
### General Project ###
# Project name
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Select column to obtain the raw counts
strand <- 4
#------------------------------------------------------------------------------------------------------------------------


# Load libraries
source(paste(path, project, "/utils/Libraries.R", sep = ""))


# Input and output directory
input_dir <- paste(path, project, "/04_STAR/", sep = "") 

# Set working directory
setwd(input_dir)

# List fastq.zip files
samples <- gsub(".out", "", list.files(pattern = ".out"))
# Filter sample name 
samples_names <- gsub("_STAR", "", samples)


#######################################################################
#                            PROCESS                           
#######################################################################


for (i in 1:length(samples_names)){
  
  # Number of reads counted per gene
  sample_counts <- read.table(paste(input_dir,"/", samples[i],"/ReadsPerGene.out.tab", sep=""), sep="\t", header = FALSE, as.is = TRUE)
  
  # Remove first four rows
  sample_counts <- sample_counts[-c(1:4),]
  
  # Select the genes names and the column of interest 
  counts <- data.frame(Genes=sample_counts$V1, sample_counts[,(strand)])
  
  # Rename the column of interest with the sample name
  colnames(counts)[2] <- paste(samples_names[i])
  
  # Merge all the samples information
  if (exists("RawCounts")) {
    RawCounts <- merge(RawCounts, counts, by="Genes")
  } else {RawCounts <- counts}
  
  ## Output
  #             Genes    CAS1 CAS2 ... CAS24
  # ENSMUSG00000051951    8    1         n
  # ENSMUSG00000089699    0    9         n
  # ENSMUSG00000102331    7    0         n
}

# No duplicate genes
any(duplicated(counts$Genes))

# Genes as rownames and sorted
rownames(RawCounts) <- RawCounts$Genes
RawCounts <- RawCounts[,-c(1)]
RawCounts <- RawCounts[order(row.names(RawCounts)),]

# Sorting by colnames 
RawCounts <- RawCounts[order(colnames(RawCounts))]


# Create a CPM
# 
# Input data must be a matrix/data frame with numeric values.
# log = TRUE; prior.count should be selected 0 for TPM, 0.25 for CPM and FPKM
# geneLength is used for length-normalized units such as TPM, FPKM or FPK
raw <- apply(RawCounts, 2, as.numeric)
counts_cpm <- convertCounts(raw, unit = "CPM", normalize = "none", log = FALSE)
rownames(counts_cpm) <- rownames(RawCounts)



#######################################################################
#                         SAVED DATA                
#######################################################################


# Save the count data frame
write.table(RawCounts, paste(input_dir, "/RawCounts_", project,".txt", sep=""), quote=F, row.names=T, sep="\t") 
write.table(counts_cpm, paste(input_dir, "/CPM_Counts_", project,".txt", sep=""), quote=F, row.names=T, sep="\t")



################################################################################
#                                 LOG FILE 
################################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$Directory <- input_dir
log_data$project_name <- project
log_data$Strandness <- strand

write.table(as.data.frame(log_data), paste(path, project, "/log/4_STAR_GeneCounts_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")


