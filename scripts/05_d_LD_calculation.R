#' ---
#' title: "Analysis of LD between independent SNPs per loci"
#' subtitle: ""
#' author: "Katrin Horn"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' In this script the LD between SNP pairs is calculated to find per phenotype independent SNPs per loci.
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(projectpath_main)

#' # Read in list of independent SNPs
#' ***
snpList = fread("../results/05_b_Cojo_Select_Results.txt")
snpList

#' # Sort data by region and remove regions with just one SNP
#' ***
setkey(snpList, region)

tab = table(snpList[, region])
tab = tab[tab != 1]
tab
toDo = names(tab)


#' # define UKBB data files and run analyis
#####
data.UKBB = "../results/05_a_ukb20272_imp_chrX_v3_s486631"

result = foreach(t = toDo) %do% {
  dat = snpList[region == t, ]
  
  #calculate pairs of SNPs
  allPairs = combn(x = c(1:nrow(dat)), m = 2)
  
  #calculate LD for every pair of SNPs
  res = foreach(c = c(1:ncol(allPairs))) %do% {
    pos1 = allPairs[1, c]
    pos2 = allPairs[2, c]
    mySNP1 = dat[pos1, SNP]
    mySNP2 = dat[pos2, SNP]
    
    outFile = paste0("../temp/LD_calculation/05_d_", gsub(":", "_", mySNP1), "_", gsub(":", "_", mySNP2))
    plinkCall = paste0(path_plink2, " --bfile ", data.UKBB, " --ld ", mySNP1, " ", mySNP2, " --out ", outFile)
    system(plinkCall)
    
    #files have different length
    linenumber = system(paste0("grep -n r^2 ", outFile, ".log"), intern = T)
    linenumber = as.numeric(substring(linenumber, 1, 2)) - 1
    output = fread(paste0(outFile, ".log"), skip = linenumber, nrows = 1)
    dummy = list(region = t, mySNP1, mySNP2, var1 = output[, V1], r2 = output[, V3], var2 = output[, V4], Dprime = output[, V6])
    dummy
  }
  res = rbindlist(res)
  res
}
result = rbindlist(result)
result

#' # Save results 
# ***
write.table(result, file = "../results/05_d_results_LD_calculation.txt", col.names = T, row.names = F, quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

