################################################################
#                       FUNCTIONS
################################################################


# Summary
# ----------
#
# Contain all the functions from the script 00_Sample_info.R


# From pdf to data frame
pdf_to_tab <- function(x, last_sample){
  
  # Description
  # 
  # Transform the table from the pdf into a data frame to select the variables of 
  # interest
  # This can not be completely automatized because the tables differ from projects
  
  # Create a character vector using \n as the row separator
  x <- map(x, ~ str_split(.x, "\\n") %>% unlist())
  x <- reduce(x, c)
  
  # Select the beginning of the table
  # Must be always the same
  if(is.na(str_which(x, "Kibrary ID")[1])){
    tab_start <- str_which(x, "GAP ID")[1]
  } else {tab_start <- str_which(x, "Library ID")[1]}
  
  # Select the row of the end
  tab_end <- str_which(x, last_sample)[1]
  
  # Select the rows of interest
  tab <- x[(tab_start):(tab_end)]
  
  # Change column separator form spaces to "|" and remove the first four rows
  tab <- str_replace_all(tab, "\\s{2,}", "|")
  tab <- tab[-(1:4)]
  
  # Create a data frame of the table using "|" as a separator 
  # Remove the first column beacuse it's empty
  tab_df <- as.data.frame(do.call(rbind, strsplit(tab, "|", fixed = TRUE)))[,-1] 
  
  return(tab_df)
}