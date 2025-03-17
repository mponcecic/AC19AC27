################################################################################
#                               TRAJECTORIES
################################################################################


# Summary
# --------

# After running the differentially expressed genes or splicing pipelines you can 
# run this pipeline in case it is time course data. This script will unveil the 
# patterns in time-course data. This patterns will be refer as trajectories or 
# clusters.
# 
# Different approaches can be used to accomplish the task. In our case, we offer 
# a soft and hard clustering algorithms. The hard clustering method is K-Means and 
# the soft clustering method is Fuzzy C-Means (FCM) clustering. 
# You can just run both or the one you believe is better for your data.


#-------------------------------------------------------------------------------
#                             REQUIREMENTS
#-------------------------------------------------------------------------------





################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC19AC27"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
# path <- "/vols/GPArkaitz_bigdata/mponce/"
path <- "W:/mponce/"

# Date of the log file 5_DEG_qc_/4_AS_qc_XXXX.log
analysis_ID <- "20240801092404"

# Clustering methods
clustering <- c("Mfuzz", "k-means")
# clustering <- "Mfuzz"
clustering <- "k-means"

# Options: "Splicing", "DESeq2", "DESeq2_Nofilter"
analysis <- "DESeq2"

# Cluster centroids
# This parameters must be selected carefully 
# If you are running Mfuzz, DO NOT START AT 1
sec_c <- seq(2, 17, 1)

# Select the number of clusters in which you want the results
# If set NULL, the results will be only performed in the optimal number of clusters
n_clusters <- 3:9

# Level order must be numeric
# In AC19AC27, the comparison were made using noDOX, DOX24H, DOX42H so we need to change 
# it to a numeric version manually
lvl_ord <- c(0, 24, 42)
metadata_condition <- c(rep(lvl_ord, each = 4))


# Set the seed
seed = "123456"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the seed
set.seed(seed)

# Remame analysis
analysis_ID2 <- analysis_ID

# Analysis ID 
# Created using a timestap to trace all de outputs belonging to a certain setup
analysis_ID <- format(Sys.time(), "%Y%m%d%H%M%S")

# Output directory
dir_wd <- paste(path, project, sep = "")
dir.create(file.path(dir_wd, "06_TRAJECTORIES"), showWarnings = FALSE)
dir_out <- paste(dir_wd, "/06_TRAJECTORIES", sep = "")
dir.create(file.path(dir_out, analysis_ID), showWarnings = FALSE)
dir_out <- paste(dir_out, "/", analysis_ID, sep = "")


# Load libraries and functions
source(paste(dir_wd, "/utils/libraries_trajectories.R", sep = ""))
source(paste(dir_wd, "/utils/functions_trajectories.R", sep = ""))

# Validation indexes data frame
val_df <- data.frame()


if(analysis == "Splicing"){source(paste(dir_wd, "/utils/trajectories_AS_call.R", sep = ""))} else{source(paste(dir_wd, "/utils/trajectories_DEGs_call.R", sep = ""))}

print("DEG/AS in all the comparison")
print(dim(sel_counts)[1])


# Remove duplicates genes/events
norm_counts <- sel_counts
rownames(norm_counts) <- df$Name


# Color used by Mfuzz package
col <- c( "#FF8F00", "#FFA700", "#FFBF00",
          "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
          "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00",
          "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28",
          "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F",
          "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF",
          "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF",
          "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF",
          "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
          "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF",
          "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F",
          "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")



################################################################################
#                             CLUSTERING METHOD 
################################################################################

for (j in 1:length(clustering)) {
  
  # Clustering methods
  method <- clustering[j]
  print(method)
  
  ref <- paste(method, analysis, project, sep = "_")
  
  
  ##############################################################################
  #                         Create working directories
  ##############################################################################
  
  # Figures folder
  dir.create(file.path(dir_output, method), showWarnings = FALSE)
  dir_fig <- paste(dir_output ,"/", method, sep='')
 
  # Save files with comparison results separate
  dir.create(file.path(dir_fig , "Results"), showWarnings = FALSE)
  dir_outfolder <- paste(dir_fig, "/Results", sep='')
  
  
  ##############################################################################
  #                      Preprocessing data        
  ##############################################################################
  
  
  # Estimate mean/average of normalized counts per each condition
  # --------------------------------------------------------------------------
  norm_counts_av <- exprs_mean(norm_counts, sample_info, trt)
  
  # THE ROW NAMES OF PHENOTYPIC DATA AND THE COLUMN NAMES OF EXPRESSION MATRIX 
  # MUST BE SAME. IN ORDERT TO BUILD THE EXPRESSION SET  
  if((all(rownames(pheno)==colnames(norm_counts_av)))==TRUE){print("ROW NAMES OF PHENOTYPIC DATA AND THE COLUMN NAMES OF EXPRESSION MATRIX ARE EQUAL.")
  } else {norm_counts_av <- norm_counts_av[,match(rownames(pheno), colnames(norm_counts_av))]}
  
  
  
  ################################################################################
  #                           FUZZY C-MEANS CLUSTERING        
  ################################################################################
  
  if(method =="Mfuzz"){ 
    
    # Step 1: Build the expression set
    # --------------------------------------------------------------------------
    
    # Phenotypic data transformation
    # This data must be class AnnotatedDataFrame in order to be part of the 
    # expressionSet.
    # 
    # Information
    # data: phenotypic data
    # varMetadata: Description of the variables in phenotyic (optional)
    
    p_data <- new("AnnotatedDataFrame", data = pheno)
    
    # The expression set must have:
    # assayData (matrix): The rows represent events, probe set, genes, ... mentioned 
    #                 as features in Expression Set and columns represent samples.
    # 
    # Other valuable information can be added to:
    # phenoData: Phenotypic Data (AnnotatedDataFrame)
    # featureData: Information about each event (AnnotatedDataFrame)
    # experimentData: Experimental information such as researchers, lab, contact 
    #                 information, abstract, ... 
    # annotation: Character string with the annotation data
    # protocolData: Equipment-generated information and protocols (AnnotatedDataFrame)
    # ...
    
    exprs_set <- ExpressionSet(assayData = norm_counts_av, phenoData = p_data) 

    
    # Step 2: Estimate fuzziness (m)
    # -----------------------------------
    # Method proposed by Schwaemmle and Jensen
    m <- mestimate(exprs_set)
    print(m)
    
    
    # Step 3: Select the minimum number of centroids (c)
    # -----------------------------------------------------
    # Method: minimum distance (Dmin) between cluster centroid can be used 
    c <- Dmin(exprs_set, m = m, crange = sec_c, repeats = 3)
    
    # # Step 4: Optimal number of clusters
    # # -----------------------------------------------------
    # c_opt <- estimate_c(c, method)
    
   # Elbow plot
    plot_elbow <- elbow_plot(sec_c, c, c_opt)
    ggsave(file = paste("Select_cluster_number_c_", c_opt, "_", analysis_ID, ".pdf", sep = ""), path = dir_fig, height = 3, width = 3, plot = plot_elbow[[1]], bg = "white")
    # Save Elbow plot data 
    write.table(plot_elbow[[2]], file = paste(dir_outfolder, "/Elbow_plot_", ref, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE)
    
    
    for(l in n_clusters){
      
      # Select the number of clusters
      c_min <- l
      print(c_min)
      
      # Condition labels
      labs <- paste(threshold, "_m_", round(m,2),"_c_", c_min, sep = "")
      lab_folder <- paste("m_", round(m,2),"_c_", c_min, sep = "")
      
      
      ##########################################################################
      #                       Create figures folders
      ##########################################################################
      
      # Figures folder
      dir.create(file.path(dir_fig, lab_folder), showWarnings = FALSE)
      dir_c <- paste(dir_fig ,"/", lab_folder, sep='')
      
      
      # Overlap figures
      dir.create(file.path(dir_fig, "Overlap"), showWarnings = FALSE)
      dir_over <- paste(dir_fig ,"/Overlap", sep='')
      
      
      # Step 5: Fuzzy c-means 
      # ------------------------------------------------------------------------
      cl <- mfuzz(exprs_set, c = c_min, m = m)
      
      # List of the events present in each cluster
      # The threshold used is 0.03
      list_mem <- acore(exprs_set, cl = cl, min.acore = 0.03)
      
      
      # Step 6: Empty clusters
      # ------------------------------------------------------------------------
      empty <- c()
      
      for (j in 1:c_min){
        k <- length(which(cl$membership[, j] > 0.5))
        empty <- append(empty, k)
      }
      
      
      # Step 7: Validation indexes
      # ------------------------------------------------------------------------
      val <- val_indx(cl, m, c_min, norm_counts_av)
      val_df <- rbind(val_df, val)
      
      
      
      #############################################################################
      #                                Plots 
      #############################################################################
      
      
      # Overlap plot 
      over_cl <- overlap(cl)
      pdf(file = paste(dir_over, "/Overlap_", ref, "_", labs, "_", analysis_ID, ".pdf", sep = ""), width = 6, height = 6, bg = "white")
      overlap.plot(cl, overlap = over_cl, thres = 0.05)
      dev.off()
      
      
      # File for membership values to each event save in a sheet per cluster
      exc1 <- createWorkbook()
      
      data_memb <- c()
      
      for (k in 1:c_min){
        
        # Select the genes/events in the first trajectory group
        mem_data <- as.data.frame(list_mem[k])
        colnames(mem_data) <- c("Name", "Membership")
        # Total number of genes/events per trajectory group
        n_sel <- dim(mem_data)[1]
        
        # Select the valid genes/events
        idx <- which(rownames(norm_counts) %in% mem_data$Name)
        
        # Normalized/raw gene/events for all samples
        sel_norm <- norm_counts[idx, sample_info$Sample]
        
        # Stack normalized and scaled gene/events per condition
        sel_norm_av <- norm_counts_av[idx, ]
        s_norm_av <- stack(as.data.frame(sel_norm_av))
        s_norm_av$Name <- rep(mem_data$Name)
        s_norm_av$ind <- as.numeric(as.character(s_norm_av$ind))
        s_norm_av$Membership <- rep(mem_data$Membership)

        
        ####### Plots #######
        
        # Boxplot
        p_boxplot <- boxplot_tm(s_norm_av, trt, analysis)
        ggsave(filename = paste("Boxplot_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), path = dir_c, plot = p_boxplot, height = 4, width = 4, bg = "white")
        
        # Mfuzz plots
        p_fuzz <- trajectory_plot(s_norm_av, trt, k, lvl_ord, n_sel, col, analysis)
        ggsave(paste(ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), width = 4, height = 3, plot = p_fuzz[[1]], path = dir_c, bg = "white")
        ggsave(paste(ref, "_lines_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), width = 4, height = 3, plot = p_fuzz[[2]], path = dir_c, bg = "white")
        
        # Heatmap with all samples
        pc_heatmap <- heatmap_scale_plot(sel_norm, sample_info, trt, color_list)
        pdf(paste(dir_c, "/Heatmap_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
        print(pc_heatmap[[1]])
        dev.off()
        
        pdf(paste(dir_c, "/Heatmap_", ref, "_", labs, "_trajectory_", k,"_cl_cols_", analysis_ID, ".pdf" , sep = ""), height = 4, width = 4, bg = "white")
        print(pc_heatmap[[2]])
        dev.off()
        
        # Grid of the figures
        grid_plots <- grid.arrange(grobs = list(p_fuzz[[1]], p_boxplot, pc_heatmap[[1]][[4]]), layout_matrix = rbind(c(1,1, 3,3,3),c(2,2,3,3,3)))
        ggsave(paste("00_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), path = dir_c, width = 8, height = 6, plot = grid_plots, bg = "white")
        
        # # Subclustering dendogram plot 
        # dendn <- dend_plot(pc_heatmap, mem_data, k)
        # pdf(paste(dir_c, "/Dendrogram_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), height = 5, width = 7, bg = "white")
        # print(dendn[[1]])
        # dev.off()
        
        # Researcher data 
        dm <- inner_join(df, mem_data, by = c("Name"))
        dm$Cluster <- k 
        
        # Sort and save data output 
        dm <- order_output(dm, analysis)
        addWorksheet(exc1, k)
        writeData(exc1, k, dm, rowNames = FALSE)
        
        }
    
      ##########################################################################
      #                               Results
      ##########################################################################
      
      # Vector cluster
      cluster <- cl$cluster
      # Vector membership
      member <- cl$membership
      
      if (identical(names(cluster), df$Name) == TRUE){
        
        ### Create the result table
        # Data frame result
        result <- df
        
        # Add columns cluster and membership
        result$Cluster <- cluster
        result$Membership <- c()
        
        # Select the corresponding membership to each events and cluster
        for (j in 1:dim(member)[1]){
          result$Membership[j] <- member[names(cluster[j]), cluster[j]]
        }
        
        # Sort based on the membership and columns
        result <- result[order(result$Membership, decreasing = TRUE), ] 
        result <- order_output(result, analysis)
        
      } else {print("THE INFORMATION ORDER IS NO CORRECT.")}
      
      
      ################################################################################
      #                               Save data 
      ################################################################################
      
      exc <- createWorkbook()
      
      
      addWorksheet(exc, "Center")
      addWorksheet(exc, "Cluster")
      addWorksheet(exc, "Size")
      addWorksheet(exc, "Membership")
      addWorksheet(exc, "Error")
      addWorksheet(exc, "Empty_cluster")
      
      writeData(exc, "Center", cl$centers, rowNames = TRUE)
      writeData(exc, "Cluster", cluster, rowNames = TRUE)
      writeData(exc, "Size", cl$size, rowNames = TRUE)
      writeData(exc, "Membership", member, rowNames = TRUE)
      writeData(exc, "Error", cl$withinerror, rowNames = TRUE)
      writeData(exc, "Empty_cluster", empty, rowNames = TRUE)
      
      if(analysis == "Splicing"){
        file_res_name <- paste(ref, ";", labs, sep = "")
      } else{
        file_res_name <- paste(ref, ";",  mt, labs, sep = "")
      }
      
      saveWorkbook(exc, file = paste(dir_outfolder, "/", file_res_name, "_RESULTS_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)
      # Researchers output
      saveWorkbook(exc1, file = paste(dir_outfolder, "/", file_res_name, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)
      # Save researcher output in one txt
      write.table(result, file = paste(dir_outfolder, "/", file_res_name, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE, sep = " ", eol = "\t")
      
      }
    
    # Validation indexes
    # Plot
    p_val <- val_plot(val_df, c_opt)
    ggsave(filename = paste("Validation_", project, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, path = dir_fig, plot = p_val, bg = "white")
    # Save
    write.table(val_df, file = paste(dir_outfolder, "/Validation_indexes_", ref, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE, )
    
    
    
    ################################################################################
    #                             K-MEANS CLUSTERING        
    ################################################################################
    
    
  } else {
    
    # Step 1: Optimal number of clusters
    # -----------------------------------------------------
    wss_res <- wss_plot(norm_counts_av, sec_c, seed, method) 
    c_opt <- wss_res[[1]]
    
    
    # Step 2: Validation plots
    # -----------------------------------------------------
    # Elbow plot
    plot_elbow <- wss_res[[2]]
    ggsave(file = paste("Select_cluster_number_c_", c_opt, "_", analysis_ID, ".pdf", sep = ""), path = dir_fig, height = 3, width = 3, plot = plot_elbow, bg = "white")
    # Elbow plot
    plot_wss <- fviz_nbclust(norm_counts_av, kmeans, method = "wss", k.max = 25)
    ggsave(filename = paste("00_Select_cluster_number_WSS_", analysis_ID, ".pdf" , sep = ""), height = 6, width = 6, plot = plot_wss, path = dir_fig, bg = "white")
    # Silhouette plot
    plot_sil <- fviz_nbclust(norm_counts_av, kmeans, method = "silhouette", k.max = 25)
    ggsave(filename = paste("00_Select_cluster_number_Silhouette_", analysis_ID, ".pdf" , sep = ""), height = 6, width = 6, plot = plot_sil, path = dir_fig, bg = "white")
    
    
    # Data from with the results of kmeans
    comp <- data.frame()
        
    for(k in n_clusters){
      
      # Step 3: Select the number of clusters
      # ------------------------------------------------------------------------
      c_min <- k
      print(c_min)
      
      # Condition labels
      labs <- paste(threshold, "_c_", c_min, sep = "")
      
      # Figures folder
      dir.create(file.path(dir_fig, labs), showWarnings = FALSE)
      dir_c <- paste(dir_fig ,"/", labs, sep='')
      
      
      # Step 4: K-means 
      # ------------------------------------------------------------------------
      km <- kmeans(norm_counts_av, centers = c_min, iter.max = 100, nstart = 100, algorithm = "Hartigan-Wong")
      
      
      # Step 5: Validation indexes 
      # ----------------------------------------------------------------------------
      kmeans_res <-  data.frame(Cluster = c_min,
                                Total_sum_sq = km$totss,
                                Total_w_cl_sum_sq = km$tot.withinss,
                                Betw_cl_su_sq = km$betweenss, 
                                Iteration = km$iter)
      
      comp <- rbind(comp, kmeans_res)
      
      
      # Step 6: Estimate the membership values
      # ----------------------------------------------------------------------------
      # Estimate the membership value of each gene/event to all the trajectories in 
      # the cluster
      center <- km$centers
      mem_matrix <- membership_kmeans(norm_counts_av, center)
      
      
      # Create a matrix with the following information
      # 
      #                           Cluster          Name         Membership
      # WDR54_ENSG00000005448         8   WDR54_ENSG00000005448  0.1860047
      # KDM7A_ENSG00000006459        10   KDM7A_ENSG00000006459  0.2104042
      # BAIAP3_ENSG00000007516       10  BAIAP3_ENSG00000007516  0.2460961
      # HSD17B6_ENSG00000025423       5 HSD17B6_ENSG00000025423  0.2912813
      # FAM13B_ENSG00000031003        3  FAM13B_ENSG00000031003  0.2376018
      # MYOM2_ENSG00000036448        10   MYOM2_ENSG00000036448  0.5110685
      mem_all <- merge_memb(mem_matrix, km)
      
      ##############################################################################
      #                                 Plots
      ##############################################################################
      
      # K-means cluster PCA plot
      plot_pca <- fviz_cluster(km, norm_counts_av, geom = "point")+theme_DEGs
      ggsave(filename = paste("PCA_", ref, "_", labs, "_", analysis_ID, ".pdf", sep = ""), height = 6, width = 6, plot = plot_pca, path = dir_c, bg = "white")
      
      
      # File for membership values to each event save in a sheet per cluster
      exc1 <- createWorkbook()
      
      data_memb <- c()
      
      for (k in 1:c_min){
        
        # Select the genes/events in the trajectory group
        mem_data <- mem_all[which(mem_all$Cluster == k), ]
        # Total number of genes/events per trajectory group
        n_sel <- dim(mem_data)[1]
        
        # Select the valid genes/events
        idx <- which(rownames(norm_counts) %in% mem_data$Name)
        
        # Normalized/raw gene/events for all samples
        sel_norm <- norm_counts[idx, sample_info$Sample]
        # Select normalized and scaled gene/events per condition
        sel_norm_av <- norm_counts_av[idx, ]
        
        
        # Stack normalized and scaled gene/events per condition
        s_norm_av <- stack(as.data.frame(sel_norm_av))
        s_norm_av$Name <- rep(mem_data$Name)
        s_norm_av$ind <- as.numeric(as.character(s_norm_av$ind))
        s_norm_av$Membership <- rep(mem_data$Membership)
        
        
        ####### Plots #######
        
        # Boxplot
        p_boxplot <- boxplot_tm(s_norm_av, trt, analysis)
        ggsave(filename = paste("Boxplot_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), plot = p_boxplot, path = dir_c, width = 4, height = 3, bg = "white")
        
        # K-means plot
        p_kmeans <- trajectory_plot(s_norm_av, trt, k, lvl_ord, n_sel, col, analysis)
        ggsave(paste(ref, "_", labs, "_trajectory_", k,analysis_ID, ".pdf", sep = ""), width = 4, height = 3, plot = p_kmeans[[1]], path = dir_c, bg = "white")
        ggsave(paste(ref, "_lines_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), width = 4, height = 3, plot = p_kmeans[[2]], path = dir_c, bg = "white")
        ggsave(paste(ref, "_black_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), width = 4, height = 3, plot = p_kmeans[[3]], path = dir_c, bg = "white")
        
        # Heatmap
        pc_heatmap <- heatmap_scale_plot(sel_norm, sample_info, trt, color_list)
        pdf(paste(dir_c, "/Heatmap_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
        print(pc_heatmap[[1]])
        dev.off()
        
        pdf(paste(dir_c, "/Heatmap_", ref, "_", labs, "_trajectory_", k, "_cl_cols_", analysis_ID, ".pdf", sep = ""), height = 4, width = 4, bg = "white")
        print(pc_heatmap[[2]])
        dev.off()
        
        # Grid of the figures
        grid_plots <- grid.arrange(grobs = list(p_kmeans[[3]], p_boxplot, pc_heatmap[[1]][[4]]), layout_matrix = rbind(c(1,1, 3,3,3),c(2,2,3,3,3)))
        ggsave(paste("00_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), path = dir_c, width = 8, height = 6, plot = grid_plots, bg = "white")
       
        # Subclustering dendogram plot 
        dendn <- dend_plot(pc_heatmap, mem_data, k)
        pdf(paste(dir_c, "/Dendrogram_", ref, "_", labs, "_trajectory_", k, "_", analysis_ID, ".pdf", sep = ""), height = 5, width = 7, bg = "white")
        print(dendn[[1]])
        dev.off()
        
        # Add column subclustering
        mem_data <- dendn[[2]]
        # Select the genes/events in the trajectory group
        data_memb <- rbind(data_memb, mem_data)
        
        # Researcher data 
        dm <- inner_join(df, mem_data, by = "Name")
        
        # Sort and save data output 
        dm <- order_output(dm, analysis)
        colnames(dm)[which(colnames(dm) == "Membership")] <- "Dist_Centroid" 
        addWorksheet(exc1, k)
        writeData(exc1, k, dm, rowNames = FALSE)

      }
      
      
      ### Create the result table
      # Data frame result
      result <- inner_join(df, data_memb, by = "Name")
      # Sort based on the membership and columns
      result <- result[order(result$Membership, decreasing = TRUE), ] 
      result <- order_output(result, analysis)
      colnames(result)[which(colnames(result) == "Membership")] <- "Dist_Centroid" 
      
      
      ################################################################################
      #                               SAVE DATA 
      ################################################################################
      
      exc <- createWorkbook()
      
      addWorksheet(exc, "Center")
      addWorksheet(exc, "Cluster")
      addWorksheet(exc, "Size")
      addWorksheet(exc, "Membership")
      addWorksheet(exc, "Other_results")
      
      writeData(exc, "Center", km$centers, rowNames = TRUE)
      writeData(exc, "Cluster", km$cluster, rowNames = TRUE)
      writeData(exc, "Size", km$size, rowNames = TRUE)
      writeData(exc, "Membership", mem_matrix, rowNames = TRUE)
      writeData(exc, "Other_results", kmeans_res, rowNames = TRUE)
      
      if(analysis == "Splicing"){
        file_res_name <- paste(ref, ";", labs, sep = "")
      } else{
        file_res_name <- paste(ref, ";",  mt, labs, sep = "")
      }
      
      saveWorkbook(exc, file = paste(dir_outfolder, "/", file_res_name, "_RESULTS_", analysis_ID,".xlsx", sep = ""), overwrite = TRUE)
      # Researchers output
      saveWorkbook(exc1, file = paste(dir_outfolder, "/", file_res_name, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)
      # Save researcher output in one txt
      write.table(result, file = paste(dir_outfolder, "/", file_res_name, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE, sep = " ", eol = "\t")
      
    }
    
    # Save Elbow plot data 
    write.table(wss_res[[3]], file = paste(dir_outfolder, "/Elbow_plot_", ref, "_", analysis_ID, ".txt", sep = ""), row.names = FALSE)
    # Save validation indexes
    write.csv(comp, file = paste(dir_outfolder, "/Validation_indexes_", ref, "_", analysis_ID, ".csv", sep = ""), row.names = FALSE)
  }
  
  
}

# Color bar
pdf(file = paste(dir_out, "/Colorbar_", analysis_ID, ".pdf" , sep = ""), width = 10, height = 2)
maColorBar(seq(0,1,0.01), col = col, horizontal = TRUE, k = 12, main = "Membership score")
dev.off()


################################################################################
#                                 LOG FILE 
################################################################################

# Save log file information
log_data <- c()
log_data$Date <- Sys.time()
log_data$analysis_ID <- analysis_ID
log_data$dir_out <- dir_out
log_data$dir_in <- dir_in
log_data$project_name <- project
log_data$condition <- trt
log_data$Analysis <- analysis
log_data$Clustering <- clustering
log_data$threshold <- threshold
log_data$Outliers <- log_data$Outliers 
log_data$seed <- seed
if(analysis != "Splicing"){
  log_data$Variance <- mt
  log_data$condition_order <- log_data$condition_order
  log_data$Varexp <- log_data$Varexp 
  log_data$min_count <- logfile$min_count
  log_data$min_prop <- logfile$min_prop
  log_data$n_large <- logfile$n_large
  log_data$min_total <- log_data$min_total
  log_data$fdr_cutoff <- fdr_cutoff
  log_data$lfc_cutoff <- lfc_cutoff
  log_data$correction <- log_data$correction
  log_data$colordir <- log_data$colordir
  log_data$colorsh <- log_data$colorsh
  }else{
    log_data$length_cutoff <- log_data$length_cutoff
    log_data$pval_cutoff <- log_data$pval_cutoff
    log_data$DPSI_cutoff <- DPSI_cutoff
    log_data$contrast <- log_data$contrast
    log_data$colorcomplex <- log_data$logfile$colorcomplex
    log_data$colorqual <- logfile$colorqual
}
log_data$colortrt <- paste(color_list[[trt]], collapse = ",")
log_data$colorheat <- paste(color_list[["Heatmap"]], collapse = ",")
log_data$colormem <- paste(col, sep = ",")

write.table(as.data.frame(log_data), paste(dir_wd, "/log/6_trajectories_", analysis, "_", analysis_ID, ".log", sep = ""), row.names = FALSE, eol = "\r")

