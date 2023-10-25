#######################################################################
#                             GENE COUNTS                         
#######################################################################

# Summary
#---------
# 
# The alignment output, ReadsPerGene.out.tab, is manipulated to obtain 
# the raw counts file. To accomplish this aim we first need to decide 
# the directionality of our reads and which one is a refflect of the 
# coding strand.  



###############################################################
#               Experimental Information
###############################################################

# Project: AC35

# The raw data analyzed is contained in /vols/GPArkaitz_bigdata/DATA_shared/AC35  
# Mouse (24 mice)
# Bulk RNAseq paired-end
# The inserts size 40 - 300 bp
# Library average size 363 bp
# Read length 101
# Library kit TruSeq Stranded mRNA


###############################################################
#                         R Version
###############################################################

# R version is the latest in March 2023
# R version 4.2.3 


###############################################################
#                How to choose the stranded
###############################################################


## Gene Counts (ReadsPerGene.out.tab)

# Column 1. Gene ID
# Column 2. counts for unstraned RNA-seq
# Column 3. counts for 1st read strand aligned with RNA
# Column 4. counts for 2nd read strand aligned with RNA 


# First step: Stranded or unstranded?
# ------------------------------------
# 
# Based on the name of the protocol, Illumina Truseq mRNA stranded, 
# reads are stranded. The second column which corresponds to the 
# unstranded DNA is not an option. 
#
#
# Second step: Stranded or unstranded?
# ------------------------------------
#
# According to SMARTer Stranded Total RNA-Seq Kit v2, we choose 
# the fourth column (counts for 2nd read strand aligned with RNA).  
# Because the first cDNA was synthetized using SMARTScribe Reverse 
# Transcriptase. In our case, the second reads align with the coding
# DNA strand were the gen is found. 
# 
# For a deeper understanding keep reading
# 
# The DNA is double stranded which means the the strand can be 
# template/antisense or coding/sense when being classified besed on 
# the gen. The coding strand is read in the 5' -> 3' direction, meaning
# this is the direction of the gene.
# 
# The RNA polymerase synthetizes the mRNA in the  5' -> 3' direction, 
# which means that the mRNA presents the same information as the coding 
# strand. Bear this idea in mind because we align the reads towards a 
# genome and not a transcriptome.
# 
# At this point, the Illumina kit takes part. A reverse transcriptase is 
# used to obtain the first cDNA, after adding a primer in the 3' end. 
# The cDNA is read in the 3' -> 5' direction which refers to the 
# template/antisense DNA strand. 
# Remember, the reverse transcriptase added 3 nucleotides at the end 
# of the synthesis. 
# 
# Afterwards, the 5' and 3' primers are added and the reads are synthetized.
# As a result, the second read can be read in the  5' -> 3' direction and the
# first read in the 3' -> 5' direction. 
# 
# In this case, the second reads is the one presenting the reading direction 
# of the mRNA and the first strand is form the templat strand. Being this the 
# reason why it's chosen which is a consequence of reverse transcriptase. 
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
# Following the manual SMARTer Stranded Total RNA-Seq Kit v2 that 
# you can find in the documentation folder and different online 
# resources such as https://github.com/alexdobin/STAR/issues/842 
# https://www.biostars.org/p/3423/
# https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html#counting-reads-per-genes
# https://chipster.csc.fi/manual/library-type-summary.html

#######################################################################
#                            LOAD LIBRARIES                           #
#######################################################################

library(dplyr)

#######################################################################
#                         LOAD FILES AND DATA                           
#######################################################################

# Project name
project <- "AC58"

# Input and output directory
input_dir <- paste("W:/mponce/", project, "/05_DEGs/03_STAR/", sep = "") 

# Set working directory
setwd(input_dir)

# List fastq.zip files
samples <- gsub(".out", "", list.files(pattern = ".out"))
# Filter sample name 
samples_names <- gsub("_STAR", "", samples)

# Select column to obtain the raw counts
strand <- 4


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


#######################################################################
#                         SAVED DATA                
#######################################################################

# Save the count dataframe
write.table(RawCounts, paste(input_dir, "/RawCounts_", project,".txt", sep=""), quote=F, row.names=T, sep="\t") 

