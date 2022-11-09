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

#' ## Check SNPs regarding heterogeneity due to different ancestry 
topHitsInMrMega[`P-value_ancestry_het` < 5*10^-8, ]
plotSNP = topHitsInMrMega[`P-value_ancestry_het` < 5*10^-8, MarkerName]

#' ## Generate Forest Plot of chr23:152898261:A:G
metaDetails = fread("../temp/12_MR_Mega/meta_results/GWASMA_uric_acid_overall_2021-06-23_11-41-05.gz")
matched = match(plotSNP, metaDetails[, markerID])
metaDetails = metaDetails[matched, ]

#get column positions of importanct columns
posBeta = colnames(metaDetails)[which(substring(colnames(metaDetails),1,5) == "beta.")]
posSE = colnames(metaDetails)[which(substring(colnames(metaDetails),1,3) == "se.")]

#calculate mean and SE of all studies
fp.mean = as.numeric(metaDetails[1, get("posBeta"), with = F])
myFilt = !is.na(fp.mean)
fp.mean = fp.mean[!is.na(fp.mean)]
fp.se = as.numeric(metaDetails[1, get("posSE"), with = F])
fp.se = fp.se[!is.na(fp.se)]
fp.lower = fp.mean - fp.se
fp.upper = fp.mean + fp.se

#extract study names
studyname = posBeta[myFilt]
studyname = substring(studyname, 6, nchar(studyname))
dummy = strsplit(studyname, split = "_")
studyname = sapply(dummy, function(x) return(paste0(x[[1]], "_", x[[2]])))
fp.labeltext = as.matrix(data.frame(study = studyname))

#add result from MetaGWAS
fp.mean = c(fp.mean, metaDetails[1, betaFEM])
fp.lower = c(fp.lower, metaDetails[1, betaFEM] - metaDetails[1, seFEM])
fp.upper = c(fp.upper, metaDetails[1, betaFEM] + metaDetails[1, seFEM])
fp.labeltext = rbind(fp.labeltext, "MetaGWAS") 

#forest plot
fp.summary = c(rep(FALSE, length(fp.mean)-1), TRUE)

marker = metaDetails[1, markerID]
marker = gsub(":", "_", marker)

forestplot(labeltext = fp.labeltext, mean = fp.mean, lower = fp.lower, upper = fp.upper, align = "r", is.summary = fp.summary, 
           xlab = metaDetails[1, markerID])

fn = paste0("../temp/12_MR_Mega/12_MR_MEGA_", marker, ".tiff")  
tiff(filename = fn, width = 800, height = 1200, res = 300, compression = 'lzw')
par(mar=c(0,0,0,0), cex.lab = 2)
forestplot(labeltext = fp.labeltext, mean = fp.mean, lower = fp.lower, upper = fp.upper, align = "r", is.summary = fp.summary, 
           xlab = metaDetails[1, markerID])
dev.off()


#' # Check for genome wide significant results of heterogeneity due to different ancestry (P-value_ancestry_het)
MR2 = MrMegaResults[`P-value_ancestry_het` < 5*10^-8, ]
table(MR2[, setting])

#annotate this SNPs with our loci and shrink to best SNP per region
loci = fread("../results/01_Locus_Definitions.txt")
foreach(l = c(1:nrow(MR2))) %do% {
  myPos = MR2[l, Position]
  mySet = MR2[l, setting]
  mySet = unlist(strsplit(mySet, "_"))[1]
  myRegion = loci[region_start < myPos & region_end > myPos & grepl(mySet, phenotype) , region]
  if (identical(integer(0), myRegion)) { myRegion = 0}
  MR2[l, region := myRegion]
}
table(MR2[, setting], MR2[, region], dnn = c("Setting", "Region"))
#one SNP for eGFR_MALE needs to be treated separatly
MR2[MarkerName == "chr23:146968633:A:C", region := -1]

#find SNP with lowest p value per setting and region
MR3 = MR2 %>% group_by(setting, region) %>% tally()
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

#' # Look at SNPsToPlot in MR-MEGA results
done = foreach(l = c(1:nrow(SNPsToPlot))) %do% {
  x = MrMegaResults[MarkerName == SNPsToPlot[l, ID_Meta],]
  x[setting == SNPsToPlot[l, setting],]
}
done = rbindlist(done)
done
write.table(done, file = "../results/12_MR_Mega_Hits.txt", col.names = T, row.names = F, quote = F)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
