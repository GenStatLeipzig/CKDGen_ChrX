#' ---
#' title: "Analysis for credible sets of independent loci"
#' subtitle: "Credible Sets"
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
#' ## Load information about independent SNPs 
toDo = fread("../results/05_b_Cojo_Select_Results.txt")

#add flag for SNPs which have to use conditional statistics
toDo[, setting := paste0(phenotype, "_", region)]
toDo[, useCondStats := FALSE]
toDo[is.element(setting, toDo[duplicated(setting), setting]), useCondStats := TRUE]
toDo


#' ## Load data of phenotypes for SNPs which can be done with unconditional statistics
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

#' ## Load loci information
loci = fread("../results/01_Locus_Definitions.txt")

#' # Definition of credible set function
#' ***
CredibleSetFunction = function(data,filenm) {
  
  dat = copy(data)
  
  #' ### Get Prior
  quant_0.975 = quantile(x = dat$beta, probs = 0.975, na.rm = T)
  quant_0.975   
  quant_0.025 = quantile(x = dat$beta, probs = 0.025, na.rm = T)
  quant_0.025  
  range_quant = quant_0.975-quant_0.025
  prior = range_quant/(2*1.96)
  print(prior)
  print(qnorm(p = 0.975, sd = prior) - quant_0.975)
  
  #' Should be close to 0 ...
  #' 
  #' ### Calculate Bayes factors ####
  #' Bayes factors, sum of Bayes factors and posterior probability (PP)
  dat[,abf.Wakefield:=gtx::abf.Wakefield(beta = beta, se = beta_se, priorsd = prior)]
  
  sum.abf.r1 = sum(dat[,abf.Wakefield], na.rm=T)
  sum.abf.r1
  dat[, PostProb:=abf.Wakefield/sum.abf.r1]
  summary(dat[,PostProb])
  ordering = order(dat[,PostProb], decreasing = TRUE)
  dat = dat[ordering,]
  dat[, SumProb:=cumsum(PostProb)]
  dat[, prior := prior]
  dat
  
  #' ### Summary ####
  myfile = paste0("../temp/06_Credible_Sets/",filenm,".txt")
  write.table(dat,file=myfile,col.names = T, row.names = F, quote = F)
  
  return(dat)
}

#' # Prepare input data and calculate credible set (unconditional statistics) ####
#' ***
neededCol = c("Chr", "SNP", "bp", "refA", "freq", "beta", "beta_se", "p", "n")  
 
uncond = foreach(s = toDo[useCondStats == FALSE, setting]) %do% {
  #get variables needed for data preparation
  line = which(toDo[, setting] == s)  
  pheno = toDo[line, phenotype]
  myRegion = toDo[line, region]
  region.start = loci[region == myRegion, region_start]
  region.stop = loci[region == myRegion, region_end]
  
  locusData = get(pheno)
  locusData = locusData[chromosome == 23 & position >= region.start & position <= region.stop, ]
  setnames(locusData, c("chromosome", "rsID", "position", "effect_allele", "EAF", "beta", "SE", "P", "N"), neededCol)
  colsOut = setdiff(colnames(locusData), neededCol)
  locusData[, get("colsOut") := NULL]
  
  snp = toDo[line, rsID]
  result = CredibleSetFunction(locusData, paste0(s, "_", gsub(":", "_", snp)))
  
  report = list(phenotype = pheno, region = myRegion, SNP = snp, CredSet.95 = nrow(result[SumProb <= 0.95, ])+1, CredSet.99 = nrow(result[SumProb <= 0.99, ])+1, prior = unique(result[, prior]))
  report
}
uncond = rbindlist(uncond)


#' # Prepare input data and calculate credible set (conditional statistics) ####
#' ***
step25 = fread("../results/05_b_step25_Matching_Table_UKBB.txt.gz") 
step25[, SNPinRef := pmax(SNPinRef.x, SNPinRef.y, na.rm = T)]

cond = foreach(s = toDo[useCondStats == TRUE, SNP]) %do% {
  #get variables needed for data preparation
  line = which(toDo[, SNP] == s)
  if (length(line) > 1) {
    pos = which(toDo[line, useCondStats] == T)
    line = line[pos]
  }
  pheno = toDo[line, phenotype]
  myRegion = toDo[line, region]
  setting = toDo[line, setting]
  
  #load data onf conditional statistics
  inFile = paste0("../temp/05_c_Cojo_cond_results/CojoCond_",pheno, "_Region_", myRegion, "_",gsub(":", "_", s),".cma.cojo")
  locusData = fread(inFile)
  setnames(locusData, c("bC", "bC_se"), c("beta", "beta_se"))
  colsOut = setdiff(colnames(locusData), neededCol)
  locusData[, get("colsOut") := NULL]
  
  #change SNP IDs to match out IDs
  matched = match(locusData[, SNP], step25[, UKBB.ID])
  locusData[, SNPinRef := step25[matched, SNPinRef]]
  locusData[!is.na(SNPinRef), SNP := SNPinRef]
  
  snp = toDo[line, rsID]
  result = CredibleSetFunction(locusData, paste0(setting, "_", gsub(":", "_", snp)))
  
  report = list(phenotype = pheno, region = myRegion, SNP = snp, CredSet.95 = nrow(result[SumProb <= 0.95, ])+1, CredSet.99 = nrow(result[SumProb <= 0.99, ])+1, prior = unique(result[, prior]))
  report
}
cond = rbindlist(cond)


#' Save results in file
overview = rbindlist(list(uncond, cond))
myOrder = order(overview[, phenotype], overview[, region])
overview = overview[myOrder, ]

#get summary of prior used
summary(overview[, prior])

write.table(overview, file = "../results/06_SNP_numbers_Credible_Sets.txt", col.names = T, row.names = F)

#' # Generate files for annotation with GWAS Pipeline ####
#' ***
colsNeeded = c("SNP", "CredSet", "PostProb", "SumProb")
phenotypes = unique(overview[, phenotype])

lists = foreach(p = phenotypes) %do% {
  myFiles = list.files(path = "../temp/06_Credible_Sets/", pattern = p)    
  
  dat = foreach(f = myFiles) %do% {
    myFile = fread(paste0("../temp/06_Credible_Sets/",f))
    pos = min(which(myFile[,SumProb]>0.99))
    myFile = myFile[c(1:pos), ]
    dummy = unlist(strsplit(f, split = "_"))
    if (grepl("eGFR_all_9", f) | grepl("UA_all_21", f) | grepl("UA_all_22", f)) { 
      myFile[, CredSet := paste0("Region", dummy[3], "_", dummy[4])] 
    } else {
      myFile[, CredSet := paste0("Region", dummy[3])]  
    }
    
    colsOut = setdiff(colnames(myFile), colsNeeded)
    myFile[, get("colsOut"):=NULL]
    setnames(myFile, "SNP", "snp")
    setcolorder(myFile, c("snp", "CredSet", "PostProb", "SumProb"))
    myFile
  }
  dat = rbindlist(dat)
  outFile = paste0("../temp/06_Credible_Sets/", p, "_SNPs_to_annotate.txt")
  write.table(dat, file = outFile, col.names = T, row.names = F, quote = F, sep = "\t")
}
lists = rbindlist(lists)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

