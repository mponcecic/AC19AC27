################################################################################
#                               RUN DESEeq2
################################################################################


# Summary 
#------------

# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory. 
# 
# This analysis is divided in different scripts. First, perform the quality 
# control of the data (05_DEG_v1_qc). Second, perform the methods to estimate the 
# DEGs which are DESeq2 (05_DEG_v2_DESeq2), EdgeR (05_DEG_v2_EdgeR), limma-voom 
# (05_DEG_v2_limma-voom) and Wilcoxon rank-sum test (05_DEG_v2_wilcoxon). Third,
# compare the results among the different methods used using 05_DEG_v3_Comparison 
# script.
# 
# In this script, DESeq2 method is perform based on the gene count matrix resulting 
# from the previous script.


# Recommendations and warnings
# ---------------------------------
# 
# If you want all the comparison to be included in the same input for applying 
# DESeq2, EdgeR, limma-voom or Wilcoxon test, this is not possible in this script. 
# We generate a unique matrix per each contrast which only includes the samples 
# of that contrast
# 
# In the design formula, interactions are not considered.
# 
# You will find a chunk of code where you can adjust all the parameters, 
# filtering, thresholds and colors needed to perform in this analysis. This 
# variables are deeply detailed in the code. Following the code, it's time to 
# load the data sets, be careful everything is matched.
# 
# Code chunk contain between the following separators can be modify and must be 
# modify to perform the correct analysis. Two chunk will be found one referring
# to the parameter selection and the order for loading the data.


# Authors: María Ponce



################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "XXX"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/user/"
path <- "W:/user/"


# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID <- "20240208144005"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "/utils/libraries_degs.R", sep = ""))

# Load functions scripts
source(paste(path, project, "/utils/functions_degs.R", sep = ""))


# Load log file 
logfile <- read.table(paste(path, project, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)


# Output directory
dir_out <- paste(path, project, "/05_DEG_ANALYSIS/", analysis_ID, sep = "")

# Input directory. Raw gene counts  
dir_infiles <- paste(dir_out, "/Results", sep = "")
 

# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- logfile$condition

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))

# Variance sources to include in the model
# Can be set NULL
# Mouse: RIN
# Human: RIN, dv200, Age, ... 
#
# Options
# var_exp <- c("Age", "dv200")
# var_exp <- NULL
var_exp <-  unlist(strsplit(logfile$Varexp, split = ","))

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}


# Filtering parameters 

# Minimum number of counts per gene
min_total <- logfile$min_total
# Minimum number of counts per sample
min_count <- logfile$min_count
# Large sample size 
n_large <- logfile$n_large
# Proportion
min_prop <- logfile$min_prop



### Threshold criteria 

## Significance level
# Default. pvalue < 0.05 
fdr_cutoff <- logfile$fdr_cutoff

## Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$lfc_cutoff

## Threshold label
threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")

## Multiple test correction
# Options
#   - FWER: Bonferroni 
#   - FDR: Benjamini-Hochberg, and the q-value
correction <- logfile$correction

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")),
                   Direction = unlist(str_split(logfile$colordir, pattern = ",")),
                   Shared = unlist(str_split(logfile$colorsh, pattern = ",")))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")
names(color_list[[trt]]) <- lvl_ord
names(color_list[["Direction"]]) <- c("Downregulated", "Not significant", "Upregulated")

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)

# Specie 
# specie = "Human"
specie <- logfile$Organism 

# Method used to study DEGs
analysis <- "DESeq2"

# Variance transformation method used 
vsd_type <- logfile$Variance

# Summary of the different comparisons results for this analysis found in the 
# 05_DEG_ANALYSIS/Result folder
sum_res <- data.frame()


################################################################################
#                               LOAD DATA
################################################################################


## Load metadata file 
sample_info <- read.table(paste(dir_infiles, "/Metadata_", project, "_", analysis_ID, ".txt", sep = ""))

## Load annotation data
# Genome annotation path
anot_path <- sub(pattern = "/.*", replacement = "", path)
if (anot_path == ""){anot_path <- "/vols/GPArkaitz_bigdata"}
# Reference genome
ref_genome <- read.table(paste(anot_path, "/DATA_shared/Genomes_Rocky/DEG_Annotation/Annotated_Genes_20231121.txt", sep = ""), header = TRUE)
genome <- ref_genome[which(ref_genome$Specie == specie), -4]
genome$Name <- paste(genome$Symbol, genome$Ensembl, sep = "_")
print(dim(genome))



################################################################################
#                   DESeq with filtered or raw data
################################################################################

for (h in 1:2) {
  
  if(h == 1){
    # FILTERED COUNTS
    analysis <- "DESeq2"
    # Load data
    raw_counts <- read.table(paste(dir_infiles, "/GeneCount_", vsd_type , "_blindFALSE_", project, "_filtered_", analysis_ID, ".txt", sep = ""))
    
  } else {
    # RAW COUNTS
    analysis <- "DESeq2_NoFilter"
    # Load data
    raw_counts <- read.table(paste(dir_infiles, "/GeneCount_", vsd_type , "_blindFALSE_", project, "_nofilter_", analysis_ID, ".txt", sep = ""))
  }
  print(analysis)
  print(dim(raw_counts))
  
  
  ################################################################################
  #                             DATA PROCESSING
  ################################################################################
  
  
  # Annotate all the genes in the matrix with the Symbol annotation
  
  # CAUTION: Symbol id present more than one Ensembl identifier. 
  # I propose to perform the analysis using the Ensembl identifier, if not the 
  # mean value of the genes with the same Symbol id should be performed.
  annot <- genome[which(genome$Ensembl %in% rownames(raw_counts)),]
  gene_names <- annot[match(rownames(raw_counts), annot$Ensembl), ]
  print(dim(gene_names))
  
  
  
  ################################################################################
  #                             CREATE DATAFRAMES
  ################################################################################
  
  # Create Workbook for the researchers
  # Each sheet corresponds to one of the comparison results
  # Columns: Name, Symbol, Ensembl, DEG, Direction, log2FC, padj, variables, 
  # normalized counts, raw count
  exc <- createWorkbook()
  
  # PCA data
  exc_pca <- createWorkbook()
  
  # Final results in *.txt
  # Data frame with all the comparison results, same as exc
  final_data <- data.frame()
  
  
  ################################################################################
  #                               COMPARISONS
  ################################################################################
  
  for (i in 1:length(contrast)){
    
    # Contrast 
    name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
    print(name)
    
    # Contrast levels 
    control <- contrast[[i]][3]
    experimental <- contrast[[i]][2]
    
    # Comparison levels
    comp_lvl <- c(control, experimental)
    
    # Select the contrast levels
    color_l <- color_list
    color_l[[trt]] <- color_l[[trt]][which(names(color_l[[trt]]) %in% c(experimental, control))]
    
    
    ##############################################################################
    #                         Create working directories
    ##############################################################################
    
    # Create a folder for the analysis. This folder contains several folders 
    # classified in Results and Figures. 
    
    # Load output directory
    dir.create(file.path(dir_out, name), showWarnings = FALSE)
    dir_outfolder <- paste(dir_out, "/", name, sep='')
    setwd(dir_outfolder)
    
    # Save files with comparison results separate
    dir.create(file.path(dir_outfolder , "Results"), showWarnings = FALSE)
    dir_files <- paste(dir_outfolder, "/Results", sep='')
    
    # Figures folder
    dir.create(file.path(dir_outfolder , analysis), showWarnings = FALSE)
    dir_fig <- paste(dir_outfolder ,"/", analysis, sep='')
    
    
    
    ##############################################################################
    #                              Data processing
    ##############################################################################
    
    
    # Metadata
    metadata <- sample_info[which(sample_info[[trt]] %in% comp_lvl),]
    metadata[,trt] <- factor(metadata[,trt], levels = c(comp_lvl))
    
    
    # Gene count matrix per comparison
    gene_counts <- raw_counts[, metadata$Sample]
    
    # Transformed data per comparison
    res_vst <- raw_counts[,which(colnames(raw_counts) %in% paste(vsd_type, metadata$Sample, sep = "_"))]
    colnames(res_vst) <- gsub(pattern = "VST_", replacement = "", colnames(res_vst))
    
    # Normalized counts
    res_norm <- raw_counts[,which(colnames(raw_counts) %in% paste("Norm_", metadata$Sample, sep = ""))]
    
    # DESeq2 design formula
    # Used to estimate the variance stabilization (VST) method proposed in DESeq2
    # 
    # Steps
    # 1. Check there are covariates
    # 2. Check the values of the covariate are not equal to avoid colinearity in 
    #   the model. If their values are equal the variable is not included in the 
    #   model.
    # 3. Create the design formula
    design_cond <- design_condition("DESeq2", trt, var_exp, metadata)
    print(design_cond)
    
    # Annotation for all the results 
    ref <- paste(analysis, "_", name, "_", project, sep = "")
    
    
    
    ##############################################################################
    #                                 DESeq2
    ##############################################################################
    
    
    # Step 1: Create DESeq object
    # ----------------------------------------------------------------------------
    #
    # DESeqDataSetFromMatrix contains all the information necessary to run DESeq. 
    # - Count matrix
    # - Metadata with the following columns: Sample, Time and RIN. Data associated 
    #       to the samples valuable to perform the contrast
    # - Design formula with the experimental conditions (Time) and the covariates (RIN)
    
    # The options of the SummarizedExperiment can be applied in this object.
    dds <- DESeqDataSetFromMatrix(countData = gene_counts,                     # Input: Count matrix 
                                  colData = metadata,                          # Input: Metadata
                                  design = eval(parse(text = design_cond)),    # Design matrix
                                  tidy = FALSE,                                # Default; Option: TRUE = First column of count data is the row names for the count matrix
                                  ignoreRank = FALSE                           # Default: Reserved for DEXSeq developers
    )
    
    
    # Step 2: Run DESeq
    # ----------------------------------------------------------------------------
    # 
    # Filtering by default
    #   - Genes with zero counts are remove
    #   - Gene with an extreme count outlier
    #   - Genes with low mean normalized counts
    
    dds <- DESeq(dds,                                                  # DESeqDataSet
                 test = "Wald",                                        # default; Wald Test
                 fitType = "parametric",                               # default; Dispersion estimates
                 sfType = "ratio",                                     # default; Size factors estimates
                 betaPrior = FALSE,                                    # default; betaPrior = FALSE
                 minReplicatesForReplace = 7,                          # default; Outliers replaced
                 useT = FALSE,                                         # default; Std Normal distr for Wald statistics
                 # minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,   # default;  lower bound on the estimate count while fitting the GLM
                 quiet = FALSE,                                        # default; Print a message per step performed
                 
                 ## Not used
                 # full = design(object),                                # default; Not used
                 # reduced,                                              # defaultt; Not used
                 # Parallelization
                 # parallel = FALSE,                                     # default; Not used
                 # BPPARAM = bpparam()                                   # default; Not used
    )
    
    # Fit of the dispersion estimates
    pdf(file = paste(dir_fig, "/00_Dispersion_DESeq2_", project, "_", analysis_ID, ".pdf", sep = ""), width = 6, height = 5)
    plotDispEsts(dds, genecol = "black", fitcol = "red", finalcol = "dodgerblue", legend = TRUE, log = "xy", cex = 0.45)
    dev.off()
    
    
    # Step 3: Contrast genes with Wald test
    # ----------------------------------------------------------------------------
    # 
    # Why independent filtering set TRUE?
    # 
    # Manual information:
    # 
    # The goal of independent filtering is to filter out those tests from the 
    # procedure that have no, or little chance of showing significant evidence, 
    # without even looking at their test statistic. Typically, this results in 
    # increased detection power at the same experiment-wide type I error. Here, 
    # we measure experiment-wide type I error in terms of the false discovery rate.
    # 
    # The results function of the DESeq2 package performs independent filtering 
    # by default using the mean of normalized counts as a filter statistic. A 
    # threshold on the filter statistic is found which optimizes the number of 
    # adjusted p values lower than a significance level alpha (we use the standard 
    # variable name for significance level, though it is unrelated to the dispersion 
    # parameter α). The adjusted p values for the genes which do not pass the 
    # filter threshold are set to NA
    # 
    # The default independent filtering is performed using the filtered_p function 
    # of the genefilter package. The filter threshold value and the number of 
    # rejections at each quantile of the filter statistic are available as 
    # metadata of the object returned by results.
    # 
    # Personal interpretation:
    # Independent filtering is a step applied after the analysis which increases 
    # the statistical power of the results. It estimates a parameter named theta 
    # which is used as a threshold and based on the mean of normalized counts, count 
    # normalization used un DESeq2. This value will change from one data set to 
    # another. The filtering of the adjusted pvalues is based on the gene counts 
    # which mean of normalized counts overcome the theta statistics. The ones 
    # which do not overcome this threshold are set to NA.
    
    # Log2 fold change result table for an specific comparison
    resl <- results(object = dds,                      # DESeqDataSet  
                    contrast = contrast[[i]],          # Constrast 
                    test = "Wald",                     # default; Wald Test
                    minmu = 0.5,                       # default; lower bound on the estimate count while fitting the GLM
                    
                    ## Independent filtering  
                    alpha = lfc_cutoff,                       # default = 0.1
                    pAdjustMethod = correction,        # default = "BH"
                    lfcThreshold = 0,                  # default
                    independentFiltering = TRUE,      # Independent filtering 
                    format = "DataFrame"               # default, Result format
                    # theta,                             # Quantiles at which to assess the number of rejections
                    # filter,                            # default, mean of normalized counts
                    # filterFun,                         # optional; for performing independent filtering and p-value adjustment
                    # 
                    # altHypothesis = "greaterAbs",     #c("greaterAbs", "lessAbs", "greater", "less"),
                    # listValues = c(1, -1),
                    # cooksCutoff,                      # Threshold Cook's distance
                    # 
                    
                    # 
                    # ## Not used
                    # addMLE = FALSE,                   # default; Not used
                    # tidy = FALSE,                     # default; Not used
                    # saveCols = NULL,                  # default; Not used
                    # name,                             # default; Not used
                    # 
                    # ## Parallelization
                    # parallel = FALSE,                 # default; Not used
                    # BPPARAM = bpparam(),              # default; Not used
    )
    
    
    # Results as data frame 
    res <- as.data.frame(resl)
    res <- res[match(rownames(gene_counts), rownames(res)),]
    
    # Shrinkage log2 fold change 
    res_lfcS <- as.data.frame(lfcShrink(dds, contrast = contrast[[i]], type = "ashr"))
    shrklog2FC <- res_lfcS$log2FoldChange 
    res <- cbind(res, shrklog2FC)
    
    # Change columns names to plot data 
    colnames(res) <- c("MeanExp","log2FC", "lfcSE", "stat", "pvalue", "padj", "shrklog2FC")
    
    # MA plot 
    pdf(paste(dir_fig, "/00_MA_plot_", ref, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 5)
    MA_plot(res, analysis, fdr_cutoff)
    dev.off()
    
    # Plot with the number of adjusted pvalues rejected 
    pdf(paste(dir_fig, "/00_Rejections_", ref, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4)
    plot(metadata(resl)$filterNumRej, type="b", ylab="number of rejections", xlab="quantiles of filter")
    lines(metadata(resl)$lo.fit, col="red")
    abline(v=metadata(resl)$filterTheta)
    dev.off()
    
    
    ##############################################################################
    #                            Data Processing
    ##############################################################################
    
    
    ## Results as a data frame
    res_df <- res
    res_df$Ensembl <- rownames(res_df)
    print(dim(res_df))
    
    
    ## Significative genes
    # Events with p-val NA are saved too
    # Based on the alpha and log2 FC thresholds  
    res_df$DEG[res_df$padj <= fdr_cutoff & abs(res_df$log2FC) > lfc_cutoff] <- "YES"
    res_df$DEG[res_df$padj > fdr_cutoff & abs(res_df$log2FC) <= lfc_cutoff] <- "NO"
    res_df$DEG[is.na(res_df$DEG)] <- "NO"
    
    ## Gene direction
    # Contains if DEGs are up (Upregulated) or downregulated (Downregulated)
    # Based on the alpha and log2 FC 
    res_df$Direction[res_df$padj <= fdr_cutoff & res_df$log2FC > lfc_cutoff] <- "Upregulated"
    res_df$Direction[res_df$padj <= fdr_cutoff & res_df$log2FC < -lfc_cutoff] <- "Downregulated"
    res_df$Direction[res_df$padj > fdr_cutoff & res_df$log2FC <= lfc_cutoff] <- "Not significant"
    res_df$Direction[is.na(res_df$Direction)] <- "Not significant"
    
    
    ## Annotated gene names in Symbol
    res_df <- merge(gene_names, res_df, by = "Ensembl") 
    print(head(res_df))
    print(dim(res_df))
    
    ## MERGE WITH GENE COUNTS
    # Row names to a variable
    genes <- gene_counts
    genes$Ensembl <-  rownames(genes)
    # Merge gene_counts and comparison results
    result <- merge(x = res_df, y = genes, by = "Ensembl")
    
    ## Differential expressed genes
    # Select differentially expressed genes
    df <- result[which(result$DEG == "YES"),]
    dim(df)
    
    ## Transform matrix 
    # Select the differentially expressed genes that overcame the test
    # Used to plot the data 
    m <- res_vst[which(rownames(res_vst) %in% df$Ensembl),]
    
    if(identical(rownames(m), df$Ensembl) == FALSE){m <- m[match(rownames(m), df$Ensembl),]}
    
    
    ##############################################################################
    #                                 Plot
    ##############################################################################
    
    ## HISTOGRAMS
    # Representation of the adjusted p-value and log2 fold-change for all and 
    # significant genes
    plot_hist <- hist_verif(res_df, df)
    ggsave(filename = paste("01_Histogram_verif_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = plot_hist, path = dir_fig, height = 4, width = 4, bg = "white")
    
    
    ##  CORRELATION 
    #
    # Execute a Pearson correlation which accept possible NAs.
    # The acceptance of the NA should be valuable in the future in case, we accept
    # comparisons when a sample is missing a PSI value.
    #
    # The results of the correlation matrix must be aligned with the results
    # in the heatmap and PCA.
    m0 <- m[which(rowSums(m) != 0), ] 
    pem <- cor(m0, method = "pearson", use = "na.or.complete")
    plot_h <- pheatmap(pem, color = colorRampPalette(brewer.pal(9, "Blues"))(255),
                       cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                       fontsize_row = 6, fontsize_col = 6, border_color = NA, treeheight_row = 0, treeheight_col = 0)
    
    pdf(paste(dir_fig, "/Cor_pearson_", ref, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
    print(plot_h)
    dev.off()
    
    ## PCA PLOTS
    plot_pcas <- pca_plot(m0, trt, metadata, color_l)
    
    ggsave(filename = paste("PCA_params_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[1]], path = dir_fig, height = 4, width = 4, bg = "white")
    ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
    ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
    ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = plot_pcas[[4]], path = dir_fig, height = 5, width = 6, bg = "white")
    addWorksheet(exc_pca, paste("Rotation", name, sep = "_"))
    addWorksheet(exc_pca, paste("GenesPCs", name, sep = "_"))  
    writeData(exc_pca, as.data.frame(plot_pcas[[5]]), sheet = paste("Rotation", name, sep = "_"))
    writeData(exc_pca, as.data.frame(plot_pcas[[6]]), sheet = paste("GenesPCs", name, sep = "_"))
    
    
    ## HEATMAP
    plot_heatmap <- heatmap_plot(m0, metadata, trt, color_l)
    
    pdf(paste(dir_fig, "/Heatmap_zscore_", ref, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
    print(plot_heatmap)
    dev.off()
    
    
    ## VOLCANO
    
    volcano <- volcano_plot(res_df, color_list = color_l, lfc_cutoff, fdr_cutoff)
    
    ggsave(filename = paste("Volcano_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = volcano[[1]], path = dir_fig, height = 5, width = 6, bg = "white")
    ggsave(filename = paste("Volcano_color_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = volcano[[2]], path = dir_fig, height = 5, width = 6, bg = "white")
    ggsave(filename = paste("Volcano_lablels_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = volcano[[3]], path = dir_fig, height = 5, width = 6, bg = "white")
    
    
    ## WATERFALL
    
    waterfall_pl <- waterfall_plot(df, color_l)
    waterfall_plot_top <- waterfall_top(df, color_l)
    
    ggsave(filename = paste("Waterfall_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = waterfall_pl, path = dir_fig, height = 5, width = 6, bg = "white")
    ggsave(filename = paste("Waterfall_top_genes_", ref, "_", analysis_ID, ".pdf", sep = ""), plot = waterfall_plot_top, path = dir_fig, height = 5, width = 6, bg = "white")
    
    
    
    ##############################################################################
    #                               Save data 
    ##############################################################################
    
    
    # Summary table
    sum_line <- c(name, analysis, design_cond, vsd_type, dim(gene_counts)[1], dim(res_df$DEG), sum(res_df$DEG == "YES"), sum(res_df$Direction == "Upregulated"), sum(res_df$Direction == "Downregulated"))
    sum_res <- rbind(sum_res, sum_line)
    
    
    # All results
    colnames(res_vst) <- paste(vsd_type, colnames(res_vst), sep = "_")
    data <- cbind(result, res_vst)
    data <- cbind(data, res_norm)
    data <- data %>% select(Name, Symbol, Ensembl, Biotype, DEG, Direction, log2FC, padj, shrklog2FC, MeanExp, lfcSE, stat, pvalue, everything())
    write.table(data, paste(dir_files, "/", ref, ";All_", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE)
    
    # Save data in the workbook
    addWorksheet(exc, name)
    writeData(exc, data, sheet = name)
    
    # All comparisons results
    result2 <- merge(x = res_df, y = raw_counts, by = "Ensembl")
    result2$Comparison <- name
    result2 <- result2 %>% select(Comparison, Name, Symbol, Ensembl, Biotype, DEG, Direction, log2FC, padj, shrklog2FC, MeanExp, lfcSE, stat, pvalue, everything())
    final_data <- rbind(final_data, result2)
    
  }
  
  # Save Workbook
  saveWorkbook(exc, file =  paste(dir_infiles, "/", analysis, "_", project, ";All_", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)
  
  # Save all comparisons 
  write.table(final_data, file = paste(dir_infiles, "/", analysis, "_", project, ";All_", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".txt", sep = ""), sep = " ", row.names = FALSE, col.names = TRUE)
  
  # Save selected genes in all comparison
  final_data <- final_data[which(final_data$DEG == "YES"), ]
  write.table(final_data, file = paste(dir_infiles, "/", analysis, "_", project, ";Selected_", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".txt", sep = ""), sep = " ", row.names = FALSE, col.names = TRUE)
  
  # Save PCA data
  saveWorkbook(exc_pca, file =  paste(dir_infiles, "/PCA_data_", analysis, "_", project, ";", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)
  
  
  ################################################################################
  #                                 LOG FILE 
  ################################################################################
  
  
  # Save log file information
  log_data <- c()
  log_data$Date <- Sys.time()
  log_data$analysis_ID <- analysis_ID
  log_data$Directory <- dir_out
  log_data$Analysis <- analysis
  log_data$project_name <- project
  log_data$Organism <- specie
  log_data$condition <- trt
  log_data$condition_order <- paste0(lvl_ord, collapse =",")
  log_data$Outliers <- paste(outliers, collapse = ",") 
  log_data$Varexp <- paste(var_exp, collapse = ",") 
  log_data$min_count <- logfile$min_count
  log_data$min_prop <- logfile$min_prop
  log_data$n_large <- logfile$n_large
  log_data$min_total <- log_data$min_total
  log_data$fdr_cutoff <- fdr_cutoff
  log_data$lfc_cutoff <- lfc_cutoff
  log_data$correction <- correction
  log_data$Variance <- vsd_type
  log_data$contrast <- paste(unlist(contrast), collapse = ",")
  log_data$colortrt <- paste(color_list[[1]], collapse = ",")
  log_data$colorheat <- paste(color_list[[2]], collapse = ",")
  log_data$colordir<-  paste(color_list[[3]], collapse = ",")
  log_data$colorsh <- paste(color_list[[4]], collapse = ",")
  
  write.table(as.data.frame(log_data), paste(path, project, "/log/5_DEG_v2_", analysis, "_", analysis_ID, ".log", sep = ""), row.names = FALSE, eol = "\r")
  
}


# Save Summary table 
colnames(sum_res) <- c("Comparison", "Analysis", "Design", "Transformation", "Genes", "DEGs", "Upregulated", "Downregulated")
write.csv(sum_res, paste(dir_infiles, "/Summary_tab_DESeq2_", project, "_", threshold, "_", analysis_ID, ".csv", sep = ""), row.names = FALSE)



