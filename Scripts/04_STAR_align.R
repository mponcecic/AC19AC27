#################################################################################
  #                            MAPPING WITH STAR 
#################################################################################

# Summary
#---------

# Mapped the reads to the genome of reference. 
# 
# The reads will be mapped using STAR which is a genome reference based aligner. 
# To run STAR two steps must be followed:
#  - Step 1. Build STAR genome index (00_STAR_genome_index.sh)
#  - Step 2. Mapping reads with STAR is performed in this script
#
# The output files are mapped and unmapped reads in SAM format,
# genes counts, transcripts in SAM format and the splice junctions 
# sites. 

# Folder
# Input: Project_folder/02_TRIMMED //  W:/DATA_shared/Sequencing_name/
# Output: Project_folder/04_STAR

# Requirement
# ---------------------
#
# Genome index appropriate for your data in the following folder: 
#         /vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Indexes
# 
# The genome index are available for human and mouse with read length of 51, 101
# and 151. 


#################################################################################
#                             STAR VERSION
#################################################################################


# STAR version 2.7.10a

## Manuals used 

# STAR manual 2.7.10a
# October 20, 2022
# Link: https://github.com/alexdobin/STAR

# STAR manual 2.7.0a
# January 23, 2019
# Link: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf



################################################################################
#                         PARAMETERS SELECTION
################################################################################


# Run parameters
#     --runThreadN
#         Number of threads used for the alignment, we use the 
#         number of cores available. 
# 
# Genome parameters
#     --genomeDir
#         Path to the genome indexes directory
#         PATH = /vols/GPArkaitz_bigdata/mponce/AC35/02_STAR/Mouse_101
# 
# Read parameters
#     --readFilesIn 
#         Files with the path. The files contain the reads to 
#         align. We used paired-end reads so two files are given 
#         with their specific path.
#         PATH = /vols/GPArkaitz_bigdata/mponce/AC35/01_TRIMMED/PolyG
# 
#     --readFilesCommand Uncompressed the fastq files. In this 
#         case, we use the command zcat.
# 
# Limits parameters
#     --limitBAMsortRAM 
#         Maximum available RAM (bytes) for sorting BAM. If =0, it 
#         will be set to the genome index size. The value is: 
#         10 x Genome size (GB).
# 
#         We use 30 GB of RAM. If not the following error will appear:
# 
#         EXITING: fatal error trying to allocate genome arrays, exception thrown: std::bad_alloc
#         Possible cause 1: not enough RAM. Check if you have enough RAM 28703634434 bytes
#         Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v 28703634434
# 
# Output parameters
#     --outSAMtype 
#         Specify the format of the output files if SAM/BAM, sorted 
#         by coordinates, ... 
#         Default option is SAM
#         We choose SAM
# 
#     --outReadsUnmapped 
#         Specifies that unmapped reads will be given in a separate 
#         file from th euniqueliy mapped reads. 
#         We use Fastx option.
# 
#     --outFileNamePrefix
#         Define the output path and name of the folder
#         We use STAR_samplename
# 
#     --outFilterMultimapNmax 
#         default: 10
#         Maximum number of loci the read is allowed to map to. 
#         Alignments (all of them) will be output only if the read maps 
#         to no more loci than this value. Otherwise no alignments will be output, 
#         and the read will be counted as "mapped to too many loci" in the Log.final.out 
# 
#         If --outFilterMultimapNmax  1, limits the output to just the uniquely-mapping 
#         reads (i.e. reads that confidently map to only one locus). If the genes of 
#         interest have high sequence similarity, this option will probably eliminate 
#         a large number of reads that map to multiple loci.
# 
#         We chose --outFilterMultimapNmax  1
# 
# Quantification of Annotations
#     --quantMode 
#         The type of quantification request. In our case, we want 
#         the gene counts and the transcripts, so we chose:
#         TranscriptomicSAM GeneCounts
# 
#         The result of the GeneCounts in Mapped only includes the 
#         Uniquely mapped reads.
# 
# 2/pass Mappings
#     --twopassMode 
#         We used Basic which perform a basic 2-pass mapping, with 
#         all 1st pass junctions inserted into the genome indices on 
#         the fly. 
#

# ------------------------------------------------------------------------------
# Considerations 
# ------------------------------------------------------------------------------
#
#
# Consideration 1
# -----------------
# 
# The parameter outFilterMultimapNmax only applies to genomic alignments 
# which means that when mapping the reads towards the genome of reference
# it limits the number of read that can be multimapped. In our case, we chose
# 1, so only uniquely mapped read are chosen.
# 
# This parameter does not affect GeneCounts, but it was used in previous scripts
# so we will maintain this parameter because it can be of interest in the future.
# 
# Situation of interest (i.e. paralog genes, pseudogenes, ...)
# 
# If a read maps perfectly to one location and to another with a couple of 
# mismatches, such as in the situation where the read is from a gene, but 
# has similarity to a pseudogene, is it possible to control alignments reported 
# in terms of their optimality ? From my reading of the user manual, STAR would 
# report both alignments as long as they meet the user's multiple mapping location 
# count thresholds.
#
# These alignments are not output to Aligned.out.sam if Nmap > outFilterMultimapNmax.
#
# If you wan to get just the multimappers, the best way is to allow them (i.e. 
# use --outFilterMultimapNmax 100), and then filter them from the sam file by 
# the MAPQ field: the unique mappers will have MAPQ=255, while multimappers will 
# have MAPQ<255. Also, unique mappers have NH:i:1 attribute, while multimappers 
# have NH:i:<Nmult>, where Nmult is the number of loci they map to.
# 
#
# Consideration 2
# -----------------
# 
# The SPLICING JUNCTIONS output files includes all reads which means 
# that includes the uniquely mapped and the multimapping reads.
# If you want to modify the output used the parameter --outSJfilterReads 
# with the option Unique. 
#
#
# Consideration 3
# -----------------
# 
# Reason why we performed SAM instead of sorted by coordinates BAM format. 
# The --limitBAMsortRAM parameter showed the following error:
# 
#         EXITING because of fatal ERROR: not enough memory for BAM sorting: 
#         SOLUTION: re-run STAR with at least --limitBAMsortRAM 1288754585
#         May 09 16:46:23 ...... FATAL ERROR, exiting
#         
# The solution found was to remove this parameter, which will turn into the 
# default value, 0, which would take the genome indexes size as reference. 
# This solution make sense to me because the alignment was not affected and 
# BAM file will not be used. 
# 
# Finally, we decided to maintain the parameter and change the output to SAM
# format.



################################################################################
#                       STAR Output
################################################################################

## Log files 

# Log.out Detail information about the run 

# Log.progress.out Progress job information updated each minute. 
#     Consult while running 

# Log.final.out Summary mapping statistics after mapping job is 
#     complete. The parameters of interest are: 
#       - Uniquely mapped reads % which should be always over 80%.
#         90% very good library.
#       - % of reads mapped to multiple loci which give information 
#         about the multimapping. 
#       - % of reads unmapped: which gives information about the 
#         cause of the mismatches.


## Mapped reads with BAM format (Aligned.sortedByCoord.out.bam)
# Binary BAM format, thus saving time on converting SAM files to BAM. 
# It can also sort BAM files by coordinates, which is required by 
# many downstream applications.


## Unmapped reads (Unmapped.out.mate1 and Unmapped.out.mate2)

# For paired-end reads, if a read maps as a whole, but one of the mates 
# does not map, both mates will also be output in their corresponding file. 
# To indicate the mapping status of the read mates,the following tags are 
# appended to the read name:
#   00: mates were not mapped;
#   10: 1st mate mapped, 2nd unmapped
#   01: 1st unmapped, 2nd mapped


## Gene Counts (ReadsPerGene.out.tab)

# This file is the file that will be used in the analysis of DEGs and will be 
# used in the following step

# Column 1. Gene ID
# Column 2. counts for unstraned RNA-seq
# Column 3. counts for 1st read strand aligned with RNA
# Column 4. counts for 2nd read strand aligned with RNA 


## Splice junctions (SJ.out.tab)
# High confidence collapsed splice junctions in tab-delimited format. 
# Note that STAR defines the junction start/end as intronic bases, 
# while many other software define them as exonic bases.

# Column 1. Chromosome
# Column 2. first base of the intron
# Column 3. last base of the intron
# Column 4. strand. 1:+ / 2:-
# Column 5. intron
# Column 6. Note that in 2-pass mode, junctions detected in the 1st pass 
#           are reported as annotated (1), in addition to annotated 
#           junctions from GTF. Unannotated: 0
# Column 7. Number of uniquely mapping reads crossing the junction
# Column 8. Number of multi-mapping reads crossing the junction
# Column 9. Maximum spliced alignment overhang. It should be 10 


## Transcript coordinates (Aligned.toTranscriptome.out.bam)

# These transcriptomic alignments can be used with various transcript 
# quantification software that require reads to be mapped to transcriptome, 
# such as RSEM or eXpress. 

# Note, that STAR first aligns reads to entire genome, and only then searches 
# for concordance between alignments and transcripts.This approach offers 
# certain advantages compared to the alignment to transcriptome only, by not 
# forcing the alignments to annotated transcripts. Note that --outFilterMultimapNmax
# only applies to genomic alignments. If an alignment passes this filter, it 
# is converted to all possible transcriptomic alignments and all of them are output.



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
path <- "/vols/GPArkaitz_bigdata/mponce/"
# path <- "W:/mponce/"

# Trimmed fastq/ raw fastqs
trmd <- TRUE

# Log file date from 0_Sample_info {YYYYMMDD} format
# Add when you want to work with raw fastqs
logdate <- ""

### Process information ###
partition <- "FAST"
time <- c("03:00:00")
memory <- c("40")
node <- 1
cpu <- 6
ram <- as.numeric(memory)*10^9

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Load log file 
logfile <- read.table(paste(path, project_name, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Specie
specie <- logfile$Organism

# Read insert
read <- logfile$read


# Input directories
if(trmd == TRUE){
  # Input directory
  dir_infiles <- paste(path, project_name, "/02_TRIMMED", sep = "" )
  pattern = "_1_trmd.fastq.gz"
  pattern2 = "_2_trmd.fastq.gz"
}else{
  # Load log file 
  logfile <- read.table(paste(path, project_name, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)
  
  # Input directory
  dir_infiles <- paste(logfile$filedirRocky, "/FASTQs", sep = "")
  pattern = "_1.fastq.gz" 
  pattern2 = "_2.fastq.gz"
  }

# Output directory
dir_out <- paste(path, project_name, sep = "")
# Create output directory
dir.create(file.path(dir_out,"04_STAR"))
dir_outfiles <- paste(dir_out,"/04_STAR",sep='')
# Set directory
setwd(dir_outfiles)


# Genome index
indexfolder <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Indexes"
genome_dir <- paste(indexfolder, "/", specie, "_", read, sep = "")


### FASTQs ###
# List  fastq.gz files
samples <- list.files(path=dir_infiles, pattern = pattern)
# Filter sample name 
samples_names <- gsub(pattern,"", samples)



###############################################################
#                           STAR
###############################################################

### Generate PBS ###
for (i in 1:length(samples_names)) {
  # Job name
  job_name <- paste(samples_names[i],"_STAR",sep='');
  
  # Input files
  input1 <- paste(samples_names[i], pattern, sep="")
  input2 <- paste(samples_names[i], pattern2, sep="")
  
  # Output files in this directory
  output <- paste(dir_outfiles, "/", job_name, "/", sep="")
  
  # Program Path 
  star_path <- c("source  /opt/ohpc/pub/apps/star/STAR-2.7.10a/cic-env")
  
  # Command
  command <- paste("STAR --genomeDir", genome_dir,             
                   "--runThreadN", cpu,                         
                   "--readFilesCommand zcat",                  
                   "--readFilesIn", input1, input2,   
                   "--limitBAMsortRAM", ram,
                   "--twopassMode Basic ",                     
                   "--quantMode TranscriptomeSAM GeneCounts",
                   "--outSAMtype SAM",
                   "--outReadsUnmapped Fastx",                 
                   "--outFileNamePrefix", output,              
                   "--outFilterMultimapNmax 1",                
                   sep=" ") 
   
  # SBATCH File
  filename <- paste(job_name,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    paste("#SBATCH --job-name=",job_name,sep=''),
    paste("#SBATCH --partition=",partition,sep=''),
    c("#SBATCH --ntasks=1"),
    paste("#SBATCH --nodes=", node,sep =''), 
    paste("#SBATCH --cpus-per-task=",cpu,sep=''),
    paste("#SBATCH --time=",time,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",job_name,".out",sep=''),
    paste("#SBATCH -e ",job_name,".err",sep=''),
    "\n\n",
    c(star_path), 
    "\n",
    c(paste("cd ", dir_infiles)), 
    "\n",
    c(command),
    file=filename,sep = "\n",append=F)
  
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
}



#######################################################################
#                            LOG FILE                        
#######################################################################

# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
logfile$Date <- Sys.time()
logfile$Trimming <- trmd

write.table(as.data.frame(logfile), paste(dir_outfiles, "/4_STAR_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



q()
