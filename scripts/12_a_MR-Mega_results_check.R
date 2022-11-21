#' ---
#' title: "Check of MR-MEGA results"
#' subtitle: "Analyses run by Kamal and Alexander Teumer"
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
#' Load data of analysed phenotypes in MR-MEGA
eGFR_ALL = fread("../temp/12_MR_Mega/results_combined/MrMega_eGFR_Chr_X_PCs_AutosomalMRMEGAV2.result")
eGFR_FEMALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_eGFR_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_F_PCs_MRMEGAV2.result")
eGFR_MALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_eGFR_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_M_PCs_MRMEGAV2.result")
UA_ALL = fread("../temp/12_MR_Mega/results_combined/MrMega_Urate_overall_Chr_X_PCs_AutosomalMRMEGAV2.result")
UA_FEMALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_Urate_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_F_PCs_MRMEGAV2.result")
UA_MALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_Urate_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_M_PCs_MRMEGAV2.result")
BUN_ALL = fread("../temp/12_MR_Mega/results_combined/MrMega_BUN_overall_Chr_X_PCs_AutosomalMRMEGAV2.result")
BUN_FEMALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_BUN_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_F_PCs_MRMEGAV2.result")
BUN_MALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_BUN_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_M_PCs_MRMEGAV2.result")
UACR_ALL = fread("../temp/12_MR_Mega/results_combined/MrMega_UACR_overall_Chr_X_PCs_AutosomalMRMEGAV2.result")
UACR_FEMALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_UACR_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_F_PCs_MRMEGAV2.result")
UACR_MALE = fread("../temp/12_MR_Mega/results_sex_stratified/MrMega_UACR_chrX_stratified_3PCs_mac10_MAF0.0025_oevarimp0.6_Autosomal_M_PCs_MRMEGAV2.result")

#' #Load data of Meta-GWAS
meta.eGFR_ALL = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
meta.eGFR_FEMALE = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
meta.eGFR_MALE = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
meta.UA_ALL = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
meta.UA_FEMALE = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
meta.UA_MALE = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
meta.BUN_ALL = fread("../data/CKDGen_ChrX_sumStat_BUN_ALL.gz")
meta.BUN_FEMALE = fread("../data/CKDGen_ChrX_sumStat_BUN_FEMALE.gz")
meta.BUN_MALE = fread("../data/CKDGen_ChrX_sumStat_BUN_MALE.gz")
meta.UACR_ALL = fread("../data/CKDGen_ChrX_sumStat_UACR_ALL.gz")
meta.UACR_FEMALE = fread("../data/CKDGen_ChrX_sumStat_UACR_FEMALE.gz")
meta.UACR_MALE = fread("../data/CKDGen_ChrX_sumStat_UACR_MALE.gz")


#' # Add invalid associations flag to MR-MEGA results and filter results
matched = match(eGFR_ALL[, MarkerName], meta.eGFR_ALL[, ID_Meta])
eGFR_ALL[, invalid_assoc := meta.eGFR_ALL[matched, invalid_assoc]]
eGFR_ALL = eGFR_ALL[invalid_assoc == F, ]
matched = match(eGFR_FEMALE[, MarkerName], meta.eGFR_FEMALE[, ID_Meta])
eGFR_FEMALE[, invalid_assoc := meta.eGFR_FEMALE[matched, invalid_assoc]]
eGFR_FEMALE = eGFR_FEMALE[invalid_assoc == F, ]
matched = match(eGFR_MALE[, MarkerName], meta.eGFR_MALE[, ID_Meta])
eGFR_MALE[, invalid_assoc := meta.eGFR_MALE[matched, invalid_assoc]]
eGFR_MALE = eGFR_MALE[invalid_assoc == F, ]

matched = match(BUN_ALL[, MarkerName], meta.BUN_ALL[, position])
BUN_ALL[, invalid_assoc := meta.BUN_ALL[matched, invalid_assoc]]
BUN_ALL = BUN_ALL[invalid_assoc == F, ]
matched = match(BUN_FEMALE[, MarkerName], meta.BUN_FEMALE[, ID_Meta])
BUN_FEMALE[, invalid_assoc := meta.BUN_FEMALE[matched, invalid_assoc]]
BUN_FEMALE = BUN_FEMALE[invalid_assoc == F, ]
matched = match(BUN_MALE[, MarkerName], meta.BUN_MALE[, ID_Meta])
BUN_MALE[, invalid_assoc := meta.BUN_MALE[matched, invalid_assoc]]
BUN_MALE = BUN_MALE[invalid_assoc == F, ]

matched = match(UA_ALL[, MarkerName], meta.UA_ALL[, ID_Meta])
UA_ALL[, invalid_assoc := meta.UA_ALL[matched, invalid_assoc]]
UA_ALL = UA_ALL[invalid_assoc == F, ]
matched = match(UA_FEMALE[, MarkerName], meta.UA_FEMALE[, ID_Meta])
UA_FEMALE[, invalid_assoc := meta.UA_FEMALE[matched, invalid_assoc]]
UA_FEMALE = UA_FEMALE[invalid_assoc == F, ]
matched = match(UA_MALE[, MarkerName], meta.UA_MALE[, ID_Meta])
UA_MALE[, invalid_assoc := meta.UA_MALE[matched, invalid_assoc]]
UA_MALE = UA_MALE[invalid_assoc == F, ]

matched = match(UACR_ALL[, MarkerName], meta.UACR_ALL[, ID_Meta])
UACR_ALL[, invalid_assoc := meta.UACR_ALL[matched, invalid_assoc]]
UACR_ALL = UACR_ALL[invalid_assoc == F, ]
matched = match(UACR_FEMALE[, MarkerName], meta.UACR_FEMALE[, ID_Meta])
UACR_FEMALE[, invalid_assoc := meta.UACR_FEMALE[matched, invalid_assoc]]
UACR_FEMALE = UACR_FEMALE[invalid_assoc == F, ]
matched = match(UACR_MALE[, MarkerName], meta.UACR_MALE[, ID_Meta])
UACR_MALE[, invalid_assoc := meta.UACR_MALE[matched, invalid_assoc]]
UACR_MALE = UACR_MALE[invalid_assoc == F, ]

#' # Combine all MR-MEGA results and sort by position
eGFR_ALL[, setting := "eGFR_ALL"]
setkey(eGFR_ALL, Position)
eGFR_FEMALE[, setting := "eGFR_FEMALE"]
setkey(eGFR_FEMALE, Position)
eGFR_MALE[, setting := "eGFR_MALE"]
setkey(eGFR_MALE, Position)
BUN_ALL[, setting := "BUN_ALL"]
setkey(BUN_ALL, Position)
BUN_FEMALE[, setting := "BUN_FEMALE"]
setkey(BUN_FEMALE, Position)
BUN_MALE[, setting := "BUN_MALE"]
setkey(BUN_MALE, Position)
UA_ALL[, setting := "UA_ALL"]
setkey(UA_ALL, Position)
UA_FEMALE[, setting := "UA_FEMALE"]
setkey(UA_FEMALE, Position)
UA_MALE[, setting := "UA_MALE"]
setkey(UA_MALE, Position)
UACR_ALL[, setting := "UACR_ALL"]
setkey(UACR_ALL, Position)
UACR_FEMALE[, setting := "UACR_FEMALE"]
setkey(UACR_FEMALE, Position)
UACR_MALE[, setting := "UACR_MALE"]
setkey(UACR_MALE, Position)

MrMegaResults = rbindlist(list(eGFR_ALL, eGFR_FEMALE, eGFR_MALE, BUN_ALL, BUN_FEMALE, BUN_MALE, 
                              UA_ALL, UA_FEMALE, UA_MALE, UACR_ALL, UACR_FEMALE, UACR_MALE))


#' # Check SNPs of table 1 in MR-MEGA results
#' ## Load SNP table 1
load("../results/04_lookup_TopHits_allOtherTraits.RData")
rsID = WideTable[, rsID]
dummy = strsplit(rsID, split =":")
newID = sapply(dummy, function(x) return(paste0("chr23:", x[2], ":", x[4], ":", x[3]))) 
newID[16] = "chr23:152898260:A:C"

#' ## Match best phenotype and setting for SNPs
a = eGFR_ALL[is.element(MarkerName, newID[c(2, 3, 5:6, 8:16)]), ]
a[, region := c(2, 3, 5:6, 8:16)]
b = eGFR_FEMALE[is.element(MarkerName, newID[c(7)]), ]
b[, region := 7]
c = eGFR_MALE[is.element(MarkerName, newID[c(1, 4)]), ]
c[, region := c(1, 4)]
d = UA_ALL[is.element(MarkerName, newID[c(17:23)]), ]
d[, region := c(17:23)]
topHitsInMrMega = rbindlist(list(a, b, c, d))
setkey(topHitsInMrMega, region)  
topHitsInMrMega = cbind(rsID, topHitsInMrMega)
topHitsInMrMega[, region := NULL]
topHitsInMrMega
write.table(topHitsInMrMega, file ="../results/12_SNPs_table1_in_MR_MEGA.txt", col.names = T, row.names = F, quote = F)


#' # Check for genome wide significant results which are not in our trans-ethnic MetaGWAS
MR2 = MrMegaResults[`P-value_association` < 5*10^-8, ]
table(MR2[, setting])

#annotate this SNPs with our loci
loci = fread("../results/01_Locus_Definitions.txt")
done = foreach(l = c(1:nrow(MR2))) %do% {
  myPos = MR2[l, Position]
  mySet = MR2[l, setting]
  mySet = unlist(strsplit(mySet, "_"))[1]
  myRegion = loci[region_start < myPos & region_end > myPos & grepl(mySet, phenotype) , region]
  if (identical(integer(0), myRegion)) { myRegion = 0.0}
  MR2[l, region := as.numeric(myRegion)]
}
done = rbindlist(done)
table(MR2[, setting], MR2[, region], dnn = c("Setting", "Region"))

#take a look at SNPs not in our loci
MR2[region == 0.0, ]

#this SNPs needs its own locus, big distance to other SNPs of UA_ALL and without locus number of metaGWAS
MR2[setting == "UA_ALL" & MarkerName == "chr23:99890204:T:C" , region := 0.1]

#find SNP with lowest p value per setting and only in region 0 (not overlapping our loci)
MR3 = MR2 %>% filter(region < 1) %>% group_by(setting, region) %>% tally()
MR3 = as.data.table(MR3)
SNPsToPlot = foreach(l = c(1:nrow(MR3))) %do% {
  mySetting = MR3[l, setting]
  myRegion = MR3[l, region]
  dat = MR2[setting == mySetting & region == myRegion, ]
  SNP = dat[`P-value_ancestry_het`== min(`P-value_ancestry_het`), MarkerName]
  if(length(SNP) > 1) {SNP = SNP[1]}
  list(setting = mySetting, region = myRegion, ID_Meta = SNP)
}
SNPsToPlot = rbindlist(SNPsToPlot)
SNPsToPlot

#add SNP from table 1 with low p-value anc_het
SNPsToPlot = rbindlist(list(list(setting = "UA_ALL", region = "22", ID_Meta = "chr23:152898261:A:G"), SNPsToPlot))
SNPsToPlot

#' # Look at SNPsToPlot in MR-MEGA results
done = foreach(l = c(1:nrow(SNPsToPlot))) %do% {
  x = MrMegaResults[MarkerName == SNPsToPlot[l, ID_Meta],]
  x[setting == SNPsToPlot[l, setting],]
}
done = rbindlist(done)
done
write.table(done, file = "../results/12_MR_Mega_Hits_including_SNP_from_table1.txt", col.names = T, row.names = F, quote = F)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
