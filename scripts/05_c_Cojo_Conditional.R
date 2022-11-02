#' ---
#' title: "Analysis for independent SNPs in all loci"
#' subtitle: "Cojo Select analysis"
#' author: "Katrin Horn, Janne Pott"
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
#' 
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(projectpath_main)

#' # Preparations ####
#' ***
#' ##  ####
#' 
#' ## Load results of COJO Select
res.select = fread("../results/05_b_Cojo_Select_Results.txt")
res.select[, checkCol := paste0(phenotype, "_", region)]

#select settings_regions with multiple independent SNPs
indep = res.select[duplicated(checkCol), checkCol]
indep

toDo = res.select[is.element(checkCol, indep), ]
toDo

#' # Prepare input data for GCTA ####
#' ***
done = foreach(l = c(1:nrow(toDo))) %do% {
  outFile = paste0("../temp/05_c_Cojo_cond_input/CojoCond_", toDo[l, phenotype], "_Region_", toDo[l, region], "_", toDo[l, SNP],".snplist")
  outFile = gsub(":", "_", outFile)
  write.table(toDo[l, SNP], file = outFile, col.names = F, row.names = F, quote = F)
}

#' # Run COJO Conditional ####
#' ***
gctaCall.1 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_eGFR_all_Region_9.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_eGFR_all_Region_9_rs181497961.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_eGFR_all_Region9_rs111410539")
system(gctaCall.1)

gctaCall.2 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_eGFR_all_Region_9.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_eGFR_all_Region_9_rs111410539.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_eGFR_all_Region9_rs181497961")
system(gctaCall.2)

gctaCall.3 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_UA_all_Region_21.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_UA_all_Region_21_X_133799101_AGT_A.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_UA_all_Region_21_rs7056552")
system(gctaCall.3)

gctaCall.4 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_UA_all_Region_21.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_UA_all_Region_21_rs7056552.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_UA_all_Region_21_X_133799101_AGT_A")
system(gctaCall.4)

gctaCall.5 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_UA_all_Region_22.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_UA_all_Region_22_rs782581283.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_UA_all_Region_22_rs4328011")
system(gctaCall.5)

gctaCall.6 = paste0(path_gcta," --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631 --chr 23 --maf 0.02 --cojo-file ../temp/05_b_Cojo_select_input/CojoInput_UA_all_Region_22.ma --cojo-cond ../temp/05_c_Cojo_cond_input/CojoCond_UA_all_Region_22_rs4328011.snplist --out ../temp/05_c_Cojo_cond_results/CojoCond_UA_all_Region_22_rs782581283")
system(gctaCall.6)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

