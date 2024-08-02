################################################################################
#                   COMPARE RESULTS BETWEEN DIFFERENT PROJECTS
################################################################################


# Project name 
project <- "AC19AC27"
project2 <- "AC19" 

# Path
path <- "W:/mponce/"

# Date of the log file 4_AS_qc_XXXX.log
analysis_ID <- "20240801092404"
analysis_ID2 <- "20240801090811"

# Comparison
comparison <- "noDOXvsDOX24H"

# Set working directories
dir_out <- paste(path, project, "/", sep = "")
dir_log <- paste(dir_out, "log", sep = "")
dir_input <- paste(dir_out, "05_DEG_ANALYSIS/", analysis_ID, "/Results/", sep = "")
dir_input2 <- paste(path, project2, "/05_DEG_ANALYSIS/", analysis_ID2, "/Results/", sep = "")

# Output directory
dir.create(file.path(dir_out, paste("06_Comparison_", project, "_", project2, sep = "")), showWarnings = FALSE)
dir_output <- paste(dir_out, "06_Comparison_", project, "_", project2, sep = "")
dir.create(file.path(dir_output, paste(analysis_ID, analysis_ID, sep = "_")))
dir_files <- paste(dir_output, "/", paste(analysis_ID, analysis_ID, sep = "_"), sep = "")
dir.create(file.path(dir_files, comparison))
dir_files <- paste(dir_files, "/", comparison, sep = "")

# Load libraries
source(paste(dir_out, "utils/libraries_DEGs.R", sep = ""))

# Load log file 
logfile <- read.table(paste(dir_log, "/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)

# Condition
trt <- logfile$condition

# Contrast levels
lvl_ord <- unlist(str_split(comparison, pattern = "vs"))

# Significance level chosen to select the differential spliced 
# events
pval_cutoff <- logfile$fdr_cutoff

# Cutoff abs(DPSI)
lfc_cutoff <- logfile$lfc_cutoff


# Outlier selection
# If a sample is considered an outlier, you can remove it 
outliers <- logfile$Outliers
if(is.na(outliers)){outliers <- NULL}

# Threshold label
threshold <- paste("padj_", pval_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")),
                   Direction = unlist(str_split(logfile$colordir, pattern = ",")),
                   Shared = unlist(str_split(logfile$colorsh, pattern = ",")),
                   Congruence = c(Yes = "#44A57CFF", No = "red"),
                   Dataset = c(AC19AC27 = "#E75B64FF", AC27 = "#D8AF39FF"))
names(color_list) <- c(trt, "Heatmap", "Direction", "Shared")
names(color_list[[trt]]) <- lvl_ord
names(color_list[["Direction"]]) <- c("Downregulated", "Not significant", "Upregulated")
color_list[[trt]] <- color_list[[trt]][which(names(color_list[[trt]]) %in% lvl_ord)]

## Functions
check_congruence <- function(x){
  a <- unlist(strsplit(x, split = ", "))
  if (all(a ==a[1])){
    return("Yes")
  } else {return("No")}
}
callback <- function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}


#### LOAD DATA
# Project 1
metadata <- read.table(paste(dir_input, "Metadata_", project, "_", analysis_ID, ".txt", sep = ""))
raw_data <- readxl::read_xlsx(paste(dir_input, "DESeq2_", project, ";All_VSTblindFALSE_", threshold, "_", analysis_ID, ".xlsx", sep = ""), sheet = comparison)
raw1 <- raw_data
sel1 <- raw_data[which(raw_data$DEG == "YES"),]

# Project 2
raw_data2 <- readxl::read_xlsx(paste(dir_input2, "DESeq2_", project2, ";All_VSTblindFALSE_", threshold, "_", analysis_ID2, ".xlsx", sep = ""), sheet = comparison)
raw2 <- raw_data2
sel2 <- raw_data2[which(raw_data2$DEG == "YES"),]


metadata <- metadata[which(metadata$Condition %in% c(lvl_ord)),]
metadata$Condition <- factor(metadata$Condition, levels = lvl_ord)
metadata2 <- metadata




## VENN DIAGRAM

# Venn diagram Spliceosome
list_events <- list(project = raw1$Name, project2 = raw2$Name)
names(list_events) <- c(project, project2)
label <- "Transcriptome"
ggvenn::ggvenn(list_events, digits = 2,stroke_size = 0.25, set_name_size = 3, text_size = 2.5, fill_alpha = 0.6, fill_color = c("#E75B64FF", "#D8AF39FF"))+ 
  labs(title = label, subtitle = comparison)+ theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
ggsave(filename = paste("Venn_diagram_", comparison,"_", label , "_", threshold,".pdf" ), path = dir_files, height = 4, width = 4, bg = "white", plot = last_plot())

# Venn diagram
list_events <- list(project = sel1$Name, project2 = sel2$Name)
names(list_events) <- c(project, project2)
label <- "DEGs"
ggvenn::ggvenn(list_events, digits = 2,stroke_size = 0.25, set_name_size = 3, text_size = 2.5, fill_alpha = 0.6, fill_color = c("#E75B64FF", "#D8AF39FF"))+ 
  labs(title = label, subtitle = comparison)+ theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
ggsave(filename = paste("Venn_diagram_", comparison,"_",  label , "_", threshold,".pdf" ), path = dir_files, height = 4, width = 4, bg = "white", plot = last_plot())



### OVERLAP DATA
data_set <- sel1 %>% mutate(Dataset = project) %>% select(Name, Symbol, Ensembl, log2FC, padj, Dataset, Direction)
data_set <- rbind(data_set, sel2 %>% mutate(Dataset = project2) %>% select(Name, Symbol, Ensembl, log2FC, padj, Dataset, Direction))
data_set <- data_set %>% mutate(Direction = case_when(Direction == "Upregulated" ~ "+", Direction  == "Downregulated" ~ "-"))
tab_dup <- as.data.frame(table(data_set$Name))
dup_event <- tab_dup[tab_dup$Freq > 1,]

duplicates_table <- data_set[which(data_set$Name %in% tab_dup$Var1),] %>% 
  group_by(Name) %>%
  summarise(N = n_distinct(Dataset),
            Dataset = paste(unique(Dataset), collapse = ", "),
            Direction = paste(Direction, collapse = ", "))
duplicates_table$Congruence <- sapply(duplicates_table$Direction, check_congruence)
duplicates_table <- as.data.frame(duplicates_table)
duplicates_table$Congruence[-which(duplicates_table$Dataset == paste(project, project2, sep = ", "))] <- NA

head(duplicates_table)
dim(duplicates_table)


ggdata <- duplicates_table %>% select(Name, Congruence)
sig2 <- raw2[which(raw2$Name %in% ggdata$Name), c("Name", "log2FC", "padj")]
sig1 <- raw1[which(raw1$Name %in% ggdata$Name), c("Name", "log2FC", "padj")]
colnames(sig1)[-1] <- paste(c("log2FC", "padj"), "_", project2, sep = "")
colnames(sig2)[-1] <- paste(c("log2FC", "padj"), "_", project, sep = "")

gdt <- merge(ggdata, sig2, by = "Name", all.x = TRUE)
gdt <- merge(gdt, sig1, by = "Name", all.x = TRUE)
dim(gdt)

# Check there is no event duplicated in the table
table(gdt$Name %in% duplicates_table$Name)
tab_dup <- as.data.frame(table(gdt$Name))
dup_event <- tab_dup[tab_dup$Freq > 1,]


gdt$Shared <- "No"
gdt$Shared[which(gdt$Name %in% duplicates_table$Name[which(duplicates_table$Dataset == paste(project, project2, sep = ", "))])] <- "Yes"
gdt$Labeled <- NA
gdt$Labeled[which(gdt$Name %in% duplicates_table$Name[which(duplicates_table$Dataset == paste(project, project2, sep = ", "))])] <- gdt$Name[which(gdt$Name %in% duplicates_table$Name[which(duplicates_table$Dataset == paste(project, project2, sep = ", "))])]

table(gdt$Congruence)
# No Yes 
# 4  29 

# Correlation 
cor_all <- cor.test(x = gdt$log2FC_AC19[-which(is.na(gdt$Congruence))], gdt$log2FC_AC19AC27[-which(is.na(gdt$Congruence))])

ggplot(gdt, aes(x = log2FC_AC19, y = log2FC_AC19AC27, color = factor(Congruence), shape = Shared, label = Labeled))+
  geom_point(size = 2.5)+
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
  geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
  geom_text_repel(hjust = "outward", max.overlaps = 10, size = 2, na.rm = TRUE, min.segment.length = 0.6)+
  # annotate("text", x = 1, y = Inf, label = paste("R=", round(cor_all$estimate,2) , "; p=", sprintf("%.2e", cor_all$p.value), sep = ""), hjust = 0, vjust = 1, )+
  scale_color_manual(values =  c("red", "#44A57CFF"))+
  labs(color ="Congruence", x = paste("log2FC", project2, sep = " "), y = paste("log2FC", project, sep = " "), title = "DEGs", subtitle = paste(comparison , "; R=", round(cor_all$estimate,2) , "; p=", sprintf("%.2e", cor_all$p.value), sep = ""))+
  theme_bw()+ 
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank(), 
        text = element_text(size = 6), axis.title = element_text(size = rel(1.25)))
ggsave(filename = paste("Correlation_", comparison,"_", label , "_", threshold,".pdf" ), path = dir_files, height = 6, width = 6, bg = "white", plot = last_plot())




#### HEATMAP



### Prepare sample information for the heatmap
sample_info <- metadata
sample_info$Dataset <- project
metadata2$Dataset <- project2
met2 <- data.frame(Sample = paste(metadata2$Sample, project2, sep = "_"), Dataset = rep(project2), Condition = metadata2$Condition) 
sample_info <- rbind(sample_info, met2)

shared_events_list <- gdt$Name[which(gdt$Shared == "Yes")]

# Column information
# Sample columns information
# Should only contain the samples of the tissue followed by the treatment
# sample_col <-  sample_info %>% dplyr::select(all_of(trt))
# rownames(sample_col) <- sample_info$Sample
# sample_col[trt] <- factor(sample_info[[trt]])
sample_col <-  sample_info %>% dplyr::select(Condition, Dataset)
rownames(sample_col) <- sample_info$Sample
sample_col$Dataset <- factor(sample_info$Dataset)
sample_col$Condition <- factor(sample_col$Condition)

label <- "DEGs"

m2 <- raw2[which(raw2$Name %in% shared_events_list), c("Name", paste("VST_", metadata2$Sample[which(metadata2$Dataset == project2)], sep=""))]
m1 <- raw1[which(raw1$Name %in% shared_events_list), c("Name", paste("VST_", sample_info$Sample[which(sample_info$Dataset == project)], sep=""))]
colnames(m2) <- gsub("VST_", "", colnames(m2))
colnames(m1) <- gsub("VST_", "", colnames(m1))
colnames(m2)[-1] <- paste(colnames(m2)[-1], project2, sep = "_")
rownames(m2) <- m2$Name
rownames(m1) <- m1$Name
m2 <- m2[,-1]
m1 <- m1[,-1]

identical(rownames(m2), rownames(m1))
m <- cbind(m1, m2)

# Row information
# Should only contain the differential events of the tissue
# # Info: Name (Gen_Event) and event type
# sample_row <- raw2[which(raw2$Name %in% gdt$Name[which(gdt$Shared == "Yes")]),] %>% distinct(Name, .keep_all = TRUE) %>% select(Name)
# sel_congruence <- gdt[which(!is.na(gdt$Congruence)), c("Name", "Congruence")]
# sample_row <- merge(sample_row, sel_congruence, by = "Name", all.y = TRUE)

color_l <- color_list

m[is.na(m)] <- 0
n <- nrow(m)

# Choosing heatmap font size
if (dim(m)[2] <= 10){font_col = 7} else if (dim(m)[2] > 10 & dim(m)[2] < 20){font_col = 5} else {font_col = 2}
if (dim(m)[1]>=60 & dim(m)[1]<200){font_row = 2} else if (dim(m)[1]>=200 & dim(m)[1]<250){font_row = 1} else if (dim(m)[1]<=60 & dim(m)[1]>=20){font_row = 2} else {font_row = 5} 

if(dim(m)[1]>= 100){
  pheatmap(m,
           scale = "row",
           color = color_l$Heatmap,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = sample_col,
           # annotation_row = sample_row,
           annotation_colors = color_l,
           show_rownames = FALSE,
           fontsize_col = font_col,
           treeheight_row = 0,
           border_color = NA,
           na_col = "#199A11",
           # annotation_legend = FALSE,
           clustering_method = "ward.D2",
           clustering_callback = callback, height = 6, width = 6) 
  # filename = paste(dir_files, "/Heatmap_Vehvskiki_", label, "_shared_", n,"_", threshold, ".pdf", sep = ""))
  
  pheatmap(m,
           scale = "row",
           color = color_l$Heatmap,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = sample_col,
           # annotation_row = sample_row,
           annotation_colors = color_l,
           show_rownames = FALSE,
           fontsize_row = font_row,
           fontsize_col = font_col,
           treeheight_row = 0,
           border_color = NA,
           na_col = "#199A11",
           # annotation_legend = FALSE,
           clustering_method = "ward.D2",
           clustering_callback = callback, height = 6, width = 6, filename = paste(dir_files, "/Heatmap_", comparison, "_", label, "_shared_", n ,"_", threshold, "_nocolscluster.pdf", sep = ""))
  
}else{
  pheatmap(m,
           scale = "row",
           color = color_l$Heatmap,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = sample_col,
           # annotation_row = sample_row,
           annotation_colors = color_l,
           show_rownames = TRUE,
           fontsize_row = font_row,
           fontsize_col = font_col,
           treeheight_row = 0,
           border_color = NA,
           na_col = "#199A11",
           # annotation_legend = FALSE,
           clustering_method = "ward.D2",
           clustering_callback = callback, height = 6, width = 6, filename = paste(dir_files, "/Heatmap_", comparison, "_", label, "_shared_", n,"_", threshold, ".pdf", sep = ""))
  
  pheatmap(m,
           scale = "row",
           color = color_l$Heatmap,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = sample_col,
           # annotation_row = sample_row,
           annotation_colors = color_l,
           show_rownames = TRUE,
           fontsize_row = font_row,
           fontsize_col = font_col,
           treeheight_row = 0,
           border_color = NA,
           na_col = "#199A11",
           # annotation_legend = FALSE,
           clustering_method = "ward.D2",
           clustering_callback = callback, height = 6, width = 6, filename = paste(dir_files, "/Heatmap_", comparison, "_", label, "_shared_", n ,"_", threshold, "_nocolscluster.pdf", sep = ""))
  
}




m2 <- raw2[which(raw2$Name %in% shared_events_list), c("Name", metadata2$Sample[which(metadata2$Dataset == project2)])]
m1 <- raw1[which(raw1$Name %in% shared_events_list), c("Name", sample_info$Sample[which(sample_info$Dataset == project)])]
colnames(m2)[-1] <- paste(colnames(m2)[-1], project2, sep = "_")
m2 <- m2[,-1]
m <- cbind(m1, m2)

data <- gdt %>% select(Name, Shared, Congruence, log2FC_AC19, padj_AC19, log2FC_AC19AC27, padj_AC19AC27)
data <- merge(data, m, by = "Name", all.x = TRUE)

vec_names <- data$Name
data <- data %>% separate(Name, into = c("Symbol", "Ensembl"), sep = "_") %>% 
  mutate(Name = vec_names) %>% 
  select(Name, Symbol, Ensembl, everything())
write.csv(data, paste(dir_files, "/Shared_DEGs_", project, "_", project2, "_", comparison, ".csv", sep = ""), row.names = FALSE)



