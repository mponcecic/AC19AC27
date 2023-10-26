
# Library
library(RColorBrewer)

# RColorBrewer palette options
# Run display.brewer.all(colorblindFriendly = TRUE) to choose the palette you 
# want, if not the default option is Dark2

# Color list
# Better to choose the manual 
color_list <- list(Heatmap = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),
                   Direction = c(Downregulated = "#4169E1", `Not significant` = "grey", Upregulated = "#DC143C"),
                   Shared = c("#87CEEB","#228B22" ,"#32CD32","#FFD700"))

# Color list function
color_palette <- function(color_list, trt, lvl_order, palette = "Dark2"){
  
  # Description
  #
  # This function aims to automatized the generation of the color palette for 
  # the different levels in the condition contrast
  
  # RColorBrewer palette options
  # Run display.brewer.all(colorblindFriendly = TRUE) to choose the palette you 
  # want, if not the default option is Dark2
  
  # Number lists in the list
  n <- length(names(color_list))
  
  # If the number list is than 4 it means that the condition color palette was 
  # not manually added and must be added now to the color_list to visualize data
  if(n < 4){
    x <- brewer.pal(length(lvl_order), palette)
    names(x) <- lvl_order
    color_list[[trt]] <- x
    
    # To avoid problems with the color palette, we make sure the levels from the 
    # present color_list the correct levels from the condition  
  } else(if(all.equal.list(names(color_list[[trt]]), lvl_order) == TRUE){
    x <- brewer.pal(length(lvl_order), palette)
    names(x) <- lvl_order
    color_list[[trt]] <- x
  })
  
  return(color_list)
}


# Try
color_list <- color_palette(color_list, trt, lvl_order, "PuOr")
