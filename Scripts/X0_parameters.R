################################################################################
#                           PARAMETER SELECTION                       
################################################################################


# Project name
project <- "AC58"

# Pathway to the folders and files
# Can be your personal folder in BigData
path <- "W:/mponce/"
vols_path <- "/vols/GPArkaitz_bigdata/mponce"

# Experimental condition
# Choose only one condition per script
# The name should be the same as in the metadata file or sample information file
trt <- "Time"

# Contrast levels
# The order is important because it will be used to order the data, as well as, 
# to create the contrast matrix, reference of the order, plot data, ...
# 
# The first must be the reference level
lvl_ord <- c("Control", "4", "24", "48")



# Date
logdate <- format(Sys.time(), "%Y%m%d.%H%M%S")

log_data <- c()

log_data$project_name <- project
log_data$condition <- trt
log_data$condition_order <- lvl_ord
log_data$path <- path
log_data$pathRocky <- path_rocky


write.table(log_data,paste("Logfile", logdate, ".log", sep = ""), row.names = FALSE)