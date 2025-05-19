################################################################################
#                             LOAD FUNCTIONS 
################################################################################



# Estimates mean per condition and scaled the data afterwards
exprs_mean <- function(matrix, metadata, trt){
  
  # In studies with replicates, only one value corresonding to each time point must 
  # be given. For this, we compute the mean of the values among the replicates per 
  # each time point.
  #
  # First, replicates must present the same name which is the time point. For 
  # example:
  # H48_2_E1 -----> H48;  H48_2_E2 -----> H48  H48_2_E3 -----> H48
  # 
  # Second, the mean is compute and the column names of the expression set must 
  # present the same order as the metadata which is the order that we want the  
  # data to be presented in the cluster. We change the name to the one given in 
  # the metadata.
  # 
  # Finally, the data is centered and scaled by genes/events to avoid events with 
  # small values to affect the clustering.
  
  matrix <- matrix[, match(metadata$Sample, colnames(matrix))]
  colnames(matrix) <- as.vector(metadata[[trt]])
  exprs <- t(apply(matrix, 1, function(x) tapply(x, colnames(matrix), mean)))
  exprs_sc <- t(scale(t(exprs), center = TRUE, scale = TRUE))
  
  return(exprs_sc)
}



# Determine the optimal number of clusters
estimate_c <- function(data, method){
  # Determine percent of variation associated with cluster based on the minimum centroid distance
  data_sum <- data/sum(data) * 100
  # Calculate cumulative percents for each cluster based on the minimum centroid distance
  cumu <- cumsum(data_sum)
  # Determine which cluster exhibits cumulative percent greater than 90% and % variation associated with the cluster as less than 5
  co1 <- which(cumu > 90 & data_sum < 5)
  # Determine the difference between variation of PC and subsequent PC
  if(method == "FCM"){co2 <- sort(which((data[1:length(data) - 1] - data[2:length(data)]) > 0.1), decreasing = T)[1] + 1
  } else {co2 <- sort(which((data[1:length(data) - 1] - data[2:length(data)]) > 700), decreasing = T)[1] + 1}
  
  max_dim <- min(co1, co2)
  
  return(max_dim)
}



# Validation indexes estimation
val_indx <- function(cluster_df, m, c_min, exprs_norm){
  
  # Vector with membership
  member <- cluster_df$membership
  # Vector with cluster centers values
  centers <- cluster_df$centers
  
  # Number of events 
  N <- dim(member)[1]
  
  
  # Validation indexes 
  # 
  # The validation indexes estimated here are PC, MPC, PE, MPE, Fuzzy silhouette 
  # and Xie and Beni index which are specific for soft clustering algorithms
  
  # Partition coefficient
  # pc <- sum(member^m)/N
  pc <- PC(member)
  
  # Partition coefficient
  # mpc <- (c_min*pc - 1)/(c_min-1)
  mpc <- MPC(member)
  
  # Partition entropy 
  # pe <- -sum(member*log(member))/N
  pe <- PE(member)
  
  # Modified partition entropy 
  mpe <- (N*pe - 1)/(N-c_min)
  
  # Fuzzy silhouette index
  sil <- SIL.F(Xca = exprs_norm, U = member, alpha = 1)
  
  # Xie and Beni index
  xb <- XB(Xca = exprs_norm, U = member, m = m, H = centers)
  
  
  # Output data frames
  val <- data.frame(c = c_min, PC = pc, MPC = mpc, PE = pe, MPE = mpe, XB = xb, Silhouette = sil)
  return(val)
}



# Sort columns names based on the data type
order_output <- function(matrix, analysis){
  if(is.null(which(colnames(matrix) == "Subcluster")) == TRUE){
    if(analysis == "Splicing"){res <- matrix %>% select(Name, GENE, EVENT, COORD, LENGTH, FullCO, COMPLEX, Cluster, Subcluster, Membership, everything())
    } else if(analysis == "DESeq2" | analysis == "DESeq2_NoFilter" ){res <- matrix %>% select(Name, Symbol, Ensembl, Biotype, Cluster, Subcluster, Membership, everything())
    }else{res <- matrix %>% select(Name, Symbol, Ensembl, Biotype, Cluster, Subcluster, Membership, everything())}
  } else{
    if(analysis == "Splicing"){res <- matrix %>% select(Name, GENE, EVENT, COORD, LENGTH, FullCO, COMPLEX, Cluster, Membership, everything())
    } else if(analysis == "DESeq2" | analysis == "DESeq2_NoFilter" ){res <- matrix %>% select(Name, Symbol, Ensembl, Biotype, Cluster, Membership, everything())
    }else{res <- matrix %>% select(Name, Symbol, Ensembl, Biotype, Cluster, Membership, everything())}
  }
}

# Estimate membership for K-means
membership_kmeans <- function(data, center){
  
  
  # This function estimates the membership values of each event/gene for each 
  # trajectory based on the euclidean distance.
  # 
  # The euclidean distance must be estimated between the normalized/raw scaled 
  # genes/events and the trajectory center values. For this, we create two matrices. 
  # The first row in both matrices corresponds to normalized/raw scaled genes/events. 
  # The second row corresponds to i trajectory center value and k trajectory center 
  # value. 
  # 
  # First matrix 
  #           0         4         24         48
  # [1,] 0.7339740 0.9881515 -0.8295939 -0.8925316
  # [2,] 0.8020594 0.2639387 -1.3132046  0.2472065
  # 
  # Second matrix 
  #          0          4         24         48
  # [1,]  0.7339740  0.9881515 -0.8295939 -0.8925316
  # [2,] -0.8732561 -0.8257354  0.8106061  0.8883855
  # 
  # The euclidean distance in first matrix is divided by the value in the second 
  # matrix, and added to the previous distance values. 
  # 
  # tmp = tmp + 1.435967/3.425682
  # vector with the distances = 1/tmp
  # 
  # This distances were added to a vector and them divided by the total value of 
  # the distances to estimate the proportions. Finally, added to the final table 
  # with all the results:
  # 
  # Output matrix
  # 
  # Rows are the genes/events
  # Columns are the trajectories in a k-means
  # 
  #                               1          2          3          4          5         ...
  # WDR54_ENSG00000005448   0.03977859 0.09013318 0.12663437 0.03508046 0.07139745      ...
  # KDM7A_ENSG00000006459   0.03586115 0.09128077 0.16608190 0.03659859 0.11744288      ...
  # BAIAP3_ENSG00000007516  0.02835844 0.10444986 0.16155755 0.02728648 0.09775259      ...
  # HSD17B6_ENSG00000025423 0.02697266 0.11622219 0.07598818 0.02780980 0.29128127      ...
  # FAM13B_ENSG00000031003  0.03608619 0.08278078 0.23760178 0.03303345 0.07297994      ...
  # MYOM2_ENSG00000036448   0.01902411 0.07295759 0.08928430 0.01884933 0.08602846      ...
  # ...                         ...       ...         ...         ...       ...         ...
  
  
  n_clust <- dim(center)[1]
  n_row <- dim(data)[1]
  
  mb <- matrix(NA, ncol = n_clust,nrow = n_row) 
  
  for(i in 1:n_row){
    mb_i <- 0
    for (j in 1:n_clust) {
      tmp <- 0
      for (k in n_clust) {tmp <- tmp + (as.vector(dist(rbind(data[i,], center[j,])))/as.vector(dist(rbind(data[i,], center[k,]))))}
      mb_i[j] <- 1/tmp
    }
    mb[i,] <- mb_i/sum(mb_i)
  }
  
  rownames(mb) <- rownames(data)
  colnames(mb) <- rownames(center)
  
  return(u=mb)
}
  

# Create a data frame with Membership and Cluster
merge_memb <- function(membership_m, res){
  
  # Create a dataframe with the columns Name, Cluster and Membership
  
  data <- as.data.frame(res$cluster)
  colnames(data) <- "Cluster"
  data$Name <- rownames(data)

  cluster <- res$cluster
  
  u_i <- 0
  
  for (i in 1:dim(data)[1]) {
    col_cl <- cluster[i]
    u_i[i] <- membership_m[i,col_cl]
  }
  
  data$Membership <- u_i
  return(data)
}



################################################################################
#                             PLOTS
################################################################################


# Elbow plot
elbow_plot <- function(sec_c, c, c_min = 0){
  data <- data.frame(x = sec_c, y = c)
  
  if(c_min == 0){
    A <- ggplot(data, aes(x = x, y = y))+
      geom_line()+
      geom_point()+
      labs(x = "Cluster number", y = "Min. centroid distance")+
      theme_bw()
  } else{
    A <- ggplot(data, aes(x = x, y = y))+
      geom_line()+
      geom_point()+
      geom_vline(xintercept = c_min, linetype = "dashed", color = "red", size = 0.5)+
      labs(x = "Cluster number", y = "Min. centroid distance")+
      theme_bw()
  }
  
  return(list(A, data))
}



# Elbow plot for K-means
wss_plot <- function(data, sec_c, seed, method){
  # Set seed 
  set.seed(seed = seed)
  
  # Max number clusters
  nc <- max(sec_c)

  # Estimate Within the sum of squares for different K
  wss <- c(c(nrow(norm_counts_av)-1)*sum(apply(norm_counts_av,2,var)), 
           sapply(2:nc, function(i) {
             sum(kmeans(norm_counts_av, centers = i, iter.max = 100, nstart = 100, algorithm = "Hartigan-Wong")$withinss)
             }))
  
  # Estimate Calinski-Harabasz index values for different k
  ch_indices <- sapply(sec_c, function(k) {
    kmeans_results <- kmeans(norm_counts_av, centers = k)
    fpc::cluster.stats(norm_counts_av, kmeans_results$cluster)$ch
    })

  
  # Data frame
  df <- data.frame(x = 1:nc, WSS_index = wss, CH_index = c(NA,ch_indices))
  
  # Estimate the optimal number of clusters
  c_opt <- estimate_c(wss, method)
  
  
  # Elbow plot
  B <- ggplot(df, aes(x = x, y = WSS_index))+
    geom_line()+
    geom_point()+
    geom_vline(xintercept = c_opt, linetype = "dashed", color = "red", size = 0.5)+
    labs(x = "Cluster number", y = "Within groups sum of square")+
    theme_bw()
  C <- ggplot(df, aes(x = x, y = CH_index))+
    geom_line()+
    geom_point()+
    labs(x = "Cluster number", y = "Calinski-Harabasz index")+
    theme_bw()
 return(list(c_opt, B, df, C))
 }


# Plot the validation indexes
val_plot <- function(val_df, cluster){
  
  # Stack matrix
  val_m <- stack(val_df[,-1])
  # Create a new column corresponding to the cluster
  val_m$Cluster <- rep(val_df$c)
  
  ggplot(data = val_m, aes(x = Cluster, y = values, color = factor(ind)))+
    geom_point()+
    geom_line()+
    scale_x_continuous(limits = c(min(val_df$c), max(val_df$c)), breaks = seq(min(val_df$c), max(val_df$c), 1))+
    labs(y = "Values of indexes", color = "Index")+
    scale_color_manual(values = c(PE = "darkblue", PC = "darkgreen", MPE = "darkviolet", MPC = "darkorange", XB = "goldenrod1", Silhouette = "coral2"))+
    theme(text = element_text(color = "black", size = 8), panel.grid = element_blank(), axis.line = element_line(size = 0.03), axis.ticks = element_line(size = 0.03), legend.position = "top", legend.title = element_text())
  
  
}

# Heatmap Scale by row
heatmap_scale_plot <- function(m, metadata, trt, color_l, callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}){
  
  # Heatmap
  #
  # Clustering of all the differentially expressed genes with different layers
  # of information such as the treatment. First, a Z-score is performed based on 
  # gene expression among the samples. Second, it's clustered by columns and rows
  # but the rows (genes) tree is not shown.
  # 
  # The heatmap function needs a
  #   - Matrix data (m): Samples as columns and genes as rows
  #   - Annotation samples (sample_col): Treatment corresponding to each sample
  
  # Column information
  # Sample columns information
  # Should only contain the samples of the tissue followed by the treatment
  sample_col <-  metadata %>% dplyr::select(all_of(trt))
  rownames(sample_col) <- metadata$Sample
  sample_col[trt] <- factor(metadata[[trt]])
  
  # Choosing heatmap font size
  if (dim(m)[2] <= 10){font_col = 7} else if (dim(m)[2] > 10 & dim(m)[2] < 20){font_col = 5} else {font_col = 2}
  if (dim(m)[1]>=60 & dim(m)[1]<200){font_row = 2} else if (dim(m)[1]>=200 & dim(m)[1]<250){font_row = 1} else if (dim(m)[1]<=60 & dim(m)[1]>=20){font_row = 2} else {font_row = 5} 
  
  
  # Heatmap cluster by rows
  A <- pheatmap(m,
                   scale = "row",
                   color = color_l$Heatmap,
                   annotation_col = sample_col,
                   annotation_colors = color_l,
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   cluster_cols = FALSE,
                   cluster_rows = TRUE, 
                   fontsize_col = font_col,
                   border_color = NA,
                   treeheight_col = 0,
                   clustering_method = "ward.D2", 
                   clustering_callback = callback) 
  
  # Heatmap cluster by columns
  B <- pheatmap(m,
                   scale = "row",
                   color = color_l$Heatmap,
                   annotation_col = sample_col,
                   annotation_colors = color_l,
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   cluster_cols = TRUE,
                   cluster_rows = FALSE, 
                   fontsize_col = font_col,
                   border_color = NA,
                   treeheight_row = 0,
                   clustering_method = "ward.D2", 
                   clustering_callback = callback) 
  
  return(list(A,B))
}


# Boxplot
boxplot_tm <- function(s_data, trt, analysis){
  
  labs = ifelse(analysis=="Splicing", "Scaled PSI", "Scaled log2 Norm. counts")

  ggplot(data = s_data, aes(x = as.factor(ind), y = values, fill = as.factor(ind)))+
    geom_boxplot(width = 0.5, alpha = 0.6, lwd = 0.05)+
    geom_point(size = 0.8, alpha = 0.7)+
    labs(x = paste(trt, sep = ""), y = paste(labs, sep = ""))+
    scale_fill_manual(values = color_list[[trt]])+
    theme(text = element_text(color = "black", size = 8), panel.grid = element_blank(), axis.line = element_line(linewidth = 0.03), axis.ticks = element_line(size = 0.03), legend.position = "none")
}





# Trajectory plot
trajectory_plot <- function(s_data, trt, k, lvl_order, n_sel, col, analysis){
  
  # This function aims to generate the same plots as FCM but considering the 
  # time points as numeric. This will plot the data in real scale.
  # 
  # Requirements: Stack matrix with the columns ind, values, Name and Membership
  
  # Labels based on the data 
  labs = ifelse(analysis=="Splicing", "Scaled PSI", "Scaled log2 Norm. counts")
  
  # Trajectory plot adjusted to time points
  A <- ggplot(data = s_data, aes(x = ind, y = values, group = reorder(Name, Membership), color = Membership))+
    geom_line(lwd = 0.5)+
    labs(y = labs, x = "", title = paste("Cluster ", k, sep = ""), subtitle = paste("n = ", n_sel, sep = ""))+
    scale_x_continuous(limits = c(min(lvl_order), max(lvl_order)), 
                       breaks = lvl_order, 
                       expand = c(0.03,0.03))+
    scale_colour_gradientn(colors = col, limits = c(0, 1))+
    theme(text = element_text(color = "black", size = 8), plot.title = element_text(face = "bold", hjust = 0.5, size = 10), plot.subtitle = element_text(hjust = 0.5), 
          axis.line = element_line(linewidth = 0.03), axis.ticks = element_line(linewidth = 0.03), legend.position = "none")
  
  # Trajectory plot with time point lines
  B <- A+geom_vline(xintercept = lvl_order, linetype = 'dashed', col = 'grey', lwd = 0.05, alpha = 0.7)
  
  # Trajectory plot with black lines and center of the clusters in red
  C <- ggplot(data = s_data, aes(x = ind, y = values, group = reorder(Name, Membership)))+
    geom_line(lwd = 0.5)+
    stat_summary(fun.y=mean, colour = "#DC143C", geom = "line", aes(group = 1), size = 0.8)+
    labs(y = labs, x = "", title = paste("Cluster ", k, sep = ""), subtitle = paste("n = ", n_sel, sep = ""))+
    scale_x_continuous(limits = c(min(lvl_order), max(lvl_order)), 
                       breaks = lvl_order, 
                       expand = c(0.03,0.03))+
    theme(text = element_text(color = "black", size = 8), plot.title = element_text(face = "bold", hjust = 0.5, size = 10), plot.subtitle = element_text(hjust = 0.5), 
          axis.line = element_line(linewidth = 0.03), axis.ticks = element_line(linewidth = 0.03), legend.position = "none")
  
  return(list(A, B, C))
}





# Dendrogram
dend_plot <- function(pheat, data, k, sb = 3){
  
  # Hierarchical clustering
  hc <- pheat[[1]]$tree_row
  
  # Subclustering groups
  clust <- dendextend:::cutree(hc, k = sb, order_clusters_as_data = FALSE)
  branch_colors <- colorspace::rainbow_hcl(n = sb)
  legend_data <- data.frame(Cluster = 1:sb, Color = branch_colors)
  
  # Plot
  plot <- plot(as.dendrogram(hc) %>% dendextend::color_branches(k = sb), 
               lwd = 2, cex = 0.01, leaflab = "none", main = paste("Cluster ", k, sep = ""))
  rect.hclust(hc,k = sb, border = branch_colors)
  legend("topleft", legend = legend_data$Cluster, fill = legend_data$Color, title = "Clusters", xpd = TRUE)
  
  splot <- recordPlot()
  
  # Sort genes to add the column
  data["Subcluster"] <- as.vector(clust[match(rownames(data), names(clust))])
  
  return(list(splot, data))
}



pca_plot <- function(sel_counts, trt, metadata, color_l) {
  
  # PCA plots
  #
  # To perform the PCA analysis we need a matrix with the PSI values per each 
  # sample. The samples must be place in the rows and the events in the columns. 
  # This matrix is centered and scale by using prcomp().
  #
  # A data frame with the original data and the metadata must be used to plot.
  # This data is used to increase interpretation power.
  #
  # The principal components which explains the variability of the data should
  # at least sum 65% of the variance of the data.
  
  if(sum(is.na(sel_counts)) >= 1){
    
    # PCA {FactoMineR}	R Documentation
    # Principal Component Analysis (PCA)
    # Description
    # Performs Principal Component Analysis (PCA) with supplementary individuals, 
    # supplementary quantitative variables and supplementary categorical variables.
    # Missing values are replaced by the column mean.
    ## Example with missing data
    ## use package missMDA
    # nb <- estim_ncpPCA(t(m), ncp.min = 0, ncp.max = 5, method.cv = "Kfold", nbsim = 100)
    # imputed <- imputePCA(t(m), ncp = nb$ncp)
    # m_pca <- PCA(imputed$completeObs)
    # m_pca <- PCA(t(m))
    
    m_0 <- m[which(rowSums(is.na(sel_counts))==0),]
    mt <- t(m_0)
  } else {
    mt <- t(sel_counts)
  }
  
  m_pca <- prcomp(mt, scale. = TRUE)
  
  var_perc <- get_eigenvalue(m_pca)[,2]
  
  # Scree plot
  # Percentage of variances explained by each principal component
  pca_scree <- fviz_eig(m_pca, choice = "variance", addlabels = TRUE, ggtheme = theme_classic(), main = "Screeplot")
  
  # PC1 vs PC2
  pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1, 2),
                           geom.ind = c("point", "text"), 
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]],  fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC2 (", round(var_perc[2],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC3
  pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1, 3),
                           geom.ind = c("point", "text"),  
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC3 (", round(var_perc[3],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC4
  pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1, 4),
                           geom.ind = c("point", "text"),  
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = color_l[[trt]]) +
    scale_fill_manual(values = color_l[[trt]]) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC4 (", round(var_perc[4],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # How to get the contribution (%) of each parameter to a PC?
  # 
  # First: Select the "rotation" and absolute value
  # --------------------------------------------------
  # The rotation are the principal components (the eigenvectors of the covariance 
  # matrix), in the original coordinate system. Typically a square matrix (unless 
  # you truncate it by introducing tolerance) with the same number of dimensions 
  # your original data had.
  pca_abs <- abs(m_pca$rotation)
  
  # Second: Scale data between 0 and 1
  # -----------------------------------
  pca_vals <- as.data.frame(sweep(pca_abs, 2, colSums(pca_abs), "/"))
  pca_vals <- pca_vals %>% mutate(ID = rownames(pca_vals)) %>% select(ID, everything())
  
  pca_rot <- as.data.frame(m_pca$rotation)
  pca_rot <- pca_rot %>% mutate(ID = rownames(pca_rot)) %>% select(ID, everything())
  
  # Return the plots
  return(list(pca_scree, pca_1vs2, pca_1vs3, pca_1vs4, pca_rot, pca_vals))
}

