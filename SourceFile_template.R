#############################
# this is a template source file
# please change all paths accordingly
#############################

#############################
# Working directory
#############################
projectpath = "/PATH/TO/PROJECTS/CKDGen_ChrX/"
projectpath_main = "/PATH/TO/PROJECTS/CKDGen_ChrX/scripts/"
path_data = "/PATH/TO/PROJECTS/CKDGen_ChrX/data/"

#############################
# Own data to use
#############################
dataSourcePath = paste0(projectpath, "/data/")


#############################
# R library and R packages
#############################
.libPaths("/PATH/TO/RLibrary/VERSION_4.x/") 
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(toolboxH))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coloc))

#############################
# Other Tools 
#############################
path_plink2 = "/PATH/TO/TOOLS/plink2.0/20210203/unix_64/plink2"
path_gcta = "/PATH/TO/TOOLS/gcta_v1.94.0/gcta"

#############################
# Downloaded data sets 
#############################
path_GTExv8 = "/PATH/TO/DATA/GTEx_v8/"
