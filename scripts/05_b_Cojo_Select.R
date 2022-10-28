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
#' ## Load data files
eGFR_all = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
eGFR_all = eGFR_all[invalid_assoc == F, ]
eGFR_male = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
eGFR_male = eGFR_male[invalid_assoc == F, ]
eGFR_female = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
eGFR_female = eGFR_female[invalid_assoc == F, ]
UA_all = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
UA_all = UA_all[invalid_assoc == F, ]
UA_male = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
UA_male = UA_male[invalid_assoc == F, ]
UA_female = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
UA_female = UA_female[invalid_assoc == F, ]

#' ## Load loci definition
loci = fread("../results/01_Locus_Definitions.txt")

#' # Prepare input data for GCTA and run COJO SELECT ####
#' ***
myColNames = c("ID_UKBB","effect_allele","other_allele","EAF","beta","SE","P","N")
settings = c("eGFR_all", "eGFR_male", "eGFR_female", "UA_all", "UA_male", "UA_female")

registerDoMC(cores=16)

#' loop through settings and regions (specific for phenotype)
cojoSelect = foreach(s = settings) %do% {
  # s = settings[1]
  dummy = unlist(strsplit(s, split = "_"))[1]
  myRegions = loci[grepl(dummy, phenotype), region]
  
  tempresults = foreach(r = myRegions) %dopar% {
    # r = myRegions[1]
    r.start = loci[r == region, region_start]
    r.ende = loci[r == region, region_end]
    
    #shrink to region data
    regionData = copy(get(s))
    regionData = regionData[position >= r.start & position <= r.ende, ]
    
    #remove SNPs with no ID in UKBB data
    regionData = regionData[!is.na(ID_UKBB),]
    #numberSNPs = nrow(regionData)
    #numberSignifSNPs = sum(regionData[,logP] > 7.3)
    
    #rename columns
    regionData0 = copy(regionData)
    colsOut = setdiff(colnames(regionData), myColNames)
    regionData[, get("colsOut"):=NULL]
    setcolorder(regionData, myColNames)
    setnames(regionData, myColNames, c("SNP", "A1", "A2", "freq", "b", "se", "p", "N"))
    
    #save input data
    input_fn = paste0("../temp/05_b_Cojo_select_input/CojoInput_", s,"_Region_", r, ".ma")
    fwrite(regionData, file=input_fn, sep = "\t")
    system(paste0("chmod 770 ",input_fn))
    
    #do cojo select run #without parameter to be 
    output_fn = paste0("../temp/05_b_Cojo_select_results/", s, "_Region_", r)
    myCall = paste0(path_gcta,
                    " --bfile ../results/05_a_ukb20272_imp_chrX_v3_s486631",
                    " --chr 23",
                    " --maf 0.02", 
                    " --cojo-file ",input_fn,
                    " --cojo-slct --cojo-p 5e-8",
                    " --out ", output_fn)
    system(myCall)
    
    resultFile = paste0(output_fn, ".jma.cojo")
    if (file.exists(resultFile)) { 
      res = fread(resultFile) 
      matched = match(res[, SNP], regionData0[, ID_UKBB])
      res[, rsID := regionData0[matched, rsID]]
    } else { res = NA}
    result = data.table(phenotype = s, region = r, size.region = r.ende - r.start, res)
    
    result
  }
  tempresults = rbindlist(tempresults, fill = T)
  tempresults
}
cojoSelect = rbindlist(cojoSelect, fill = T)
cojoSelect = cojoSelect[!is.na(Chr), ]
cojoSelect[, res := NULL]
cojoSelect

system("chmod -R 770 ../temp/05_b_Cojo_select_results/")

#' Save results in file
write.table(cojoSelect, file = "../results/05_b_Cojo_Select_Results.txt", col.names = T, row.names = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

