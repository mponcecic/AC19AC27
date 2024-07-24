################################################################################
#                           MERGE AC19 AND AC27
################################################################################
# 
# Resume
# -------------------
# Merge together the counts from AC19 (first run) and AC27 (second run AC19). 

# Path
path <- "W:/mponce/"

# Project 
project <- "AC19AC27"

# Set working directory
dir_out <- paste(path, project, sep = "")
dir.create(file.path(dir_out, "04_STAR"), showWarnings = FALSE)
dir_output <- paste(dir_out, "/04_STAR", sep = "")

# Load gene count
counts_AC19 <- read.table("W:/mponce/AC19/04_STAR/RawCounts_AC19.txt", header = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")
counts_AC27 <- read.table("W:/mponce/AC27/04_STAR/RawCounts_AC27.txt", header = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")
print(dim(counts_AC27) == dim(counts_AC19))

# # Load Metadata
# metadata_AC19 <- read.csv("W:/mponce/AC19/Sample_info.csv")
# metadata_AC27 <- read.csv("W:/mponce/AC27/Sample_info.csv") 

# Check columns and rows in same order
stopifnot(identical(colnames(counts_AC19), colnames(counts_AC27)))
stopifnot(identical(rownames(counts_AC19), rownames(counts_AC27)))

# Merge data
gene_count <- counts_AC19 + counts_AC27

# Save data
write.table(gene_count, paste(dir_output, "/RawCounts_", project, ".txt", sep = ""), quote = FALSE, row.names = TRUE, sep = "\t")
