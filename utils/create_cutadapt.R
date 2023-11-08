################################################################################
#                       CREATE CUTADAPT ENVIRONMENT
################################################################################

# Summary
# ----------

# Generate the cutadapt environment in anaconda in your session in Rocky
# Perform this in the login node. There is no need to ask for a session. 

source /opt/ohpc/pub/apps/anaconda3/cic-env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n cutadaptenv cutadapt       