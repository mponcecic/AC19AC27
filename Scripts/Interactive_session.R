################################################################################
#                           INTERACTIVE JOB
################################################################################

# Resume
# ---------
# 
# Ask for an interactive job in Rocky and access to R paths
# 
# Notes:
# This files is based on the previous work of Saioa 
# 0_Saioa_CellRanger_count_Rocky_15092022
# 
# Date: 27/10/2023

### Interactive job ###
# Ask for a node and wait until it allocates one to you
salloc -N 1 -n 1 --mem=G -t 00:30:00 --partition=FAST --job-name=interactive 

# Enter the node you were asign to
ssh c01

# Job queues
squeue -u mponce

### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/DEG_Rocky")