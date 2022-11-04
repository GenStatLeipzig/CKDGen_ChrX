#' ---
#' title: "Lookup of our candidate SNPs in GWAS catalog"
#' subtitle: "Lookup in Graham, Sakaue and Kanai GWAS results"
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
#' Load data of three important GWAS
graham = fread("../temp/11_LookUp/data_from_GWAS_catalog/gwas-association-downloaded_2022-08-05-pubmedId 31015462.tsv")
sakaue = fread("../temp/11_LookUp/data_from_GWAS_catalog/gwas-association-downloaded_2022-08-05-pubmedId 34594039.tsv")
kanai = fread("../temp/11_LookUp/data_from_GWAS_catalog/gwas-association-downloaded_2022-08-05-pubmedId 29403010.tsv")

#Filter data for chromosome X
graham = graham[CHR_ID == "X", ]
sakaue = sakaue[CHR_ID == "X", ]
kanai = kanai[CHR_ID == "X", ]

#combine data sets
knownGWAS = rbindlist(list(graham, sakaue, kanai))
knownGWAS

#keep only kidney related phenotypes
phenotypes2keep = c("Creatinine levels", "Estimated glomerular filtration rate", "Albumin-globulin ratio", "Creatinine levels", 
                    "Glomerular filtration rate", "Non-albumin protein levels", "Nephrotic syndrome", "Serum albumin levels", 
                    "Blood urea nitrogen levels", "Serum uric acid levels", "Serum creatinine levels", "Serum total protein level")

knownGWAS = knownGWAS[is.element(`DISEASE/TRAIT`, phenotypes2keep), ]
table(knownGWAS[, `DISEASE/TRAIT`])


#' # Load Meta-GWAS results of all phenotypes and sex combined (ALL) setting
#adding rsID only for matching
eGFR = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
shortID = sapply(eGFR[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
eGFR[, shortID := shortID]

CKD = fread("../data/CKDGen_ChrX_sumStat_CKD_ALL.gz")
shortID = sapply(CKD[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
CKD[, shortID := shortID]

BUN = fread("../data/CKDGen_ChrX_sumStat_BUN_ALL.gz")
shortID = sapply(BUN[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
BUN[, shortID := shortID]

UA = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
shortID = sapply(UA[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
UA[, shortID := shortID]

Gout = fread("../data/CKDGen_ChrX_sumStat_Gout_ALL.gz")
shortID = sapply(Gout[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
Gout[, shortID := shortID]

UACR = fread("../data/CKDGen_ChrX_sumStat_UACR_ALL.gz")
shortID = sapply(UACR[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
UACR[, shortID := shortID]

MA = fread("../data/CKDGen_ChrX_sumStat_MA_ALL.gz")
shortID = sapply(MA[, rsID], function(x) return(unlist(strsplit(x, split = ":"))[1]))
MA[, shortID := shortID]


#' # Look up of known GWAS hits in our results
cols2keep = c("STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "REGION" , "CHR_ID", "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE", "STRONGEST SNP-RISK ALLELE",
              "SNPS", "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG", "OR or BETA")
colsOut = setdiff(colnames(knownGWAS), cols2keep)
knownGWAS[, get("colsOut") := NULL]

cols2keep = c("N", "effect_allele", "EAF", "infoScore", "I2", "beta", "SE", "P", "logP", "invalid_assoc", "reason_for_exclusion")

matched = match(knownGWAS[, SNPS], eGFR[, shortID])
eGFR = eGFR[matched, ]
colsOut = setdiff(colnames(eGFR), cols2keep)
eGFR[, get("colsOut") := NULL]
setcolorder(eGFR, cols2keep)
setnames(eGFR, colnames(eGFR), paste0("eGFR.", colnames(eGFR)))
knownGWAS = cbind(knownGWAS, eGFR)

matched = match(knownGWAS[, SNPS], CKD[, shortID])
CKD = CKD[matched, ]
colsOut = setdiff(colnames(CKD), cols2keep)
CKD[, get("colsOut") := NULL]
setcolorder(CKD, cols2keep)
setnames(CKD, colnames(CKD), paste0("CKD.", colnames(CKD)))
knownGWAS = cbind(knownGWAS, CKD)

matched = match(knownGWAS[, SNPS], BUN[, shortID])
BUN = BUN[matched, ]
colsOut = setdiff(colnames(BUN), cols2keep)
BUN[, get("colsOut") := NULL]
setcolorder(BUN, cols2keep)
setnames(BUN, colnames(BUN), paste0("BUN.", colnames(BUN)))
knownGWAS = cbind(knownGWAS, BUN)

matched = match(knownGWAS[, SNPS], UA[, shortID])
UA = UA[matched, ]
colsOut = setdiff(colnames(UA), cols2keep)
UA[, get("colsOut") := NULL]
setcolorder(UA, cols2keep)
setnames(UA, colnames(UA), paste0("UA.", colnames(UA)))
knownGWAS = cbind(knownGWAS, UA)

matched = match(knownGWAS[, SNPS], Gout[, shortID])
Gout = Gout[matched, ]
colsOut = setdiff(colnames(Gout), cols2keep)
Gout[, get("colsOut") := NULL]
setcolorder(Gout, cols2keep)
setnames(Gout, colnames(Gout), paste0("Gout.", colnames(Gout)))
knownGWAS = cbind(knownGWAS, Gout)

matched = match(knownGWAS[, SNPS], UACR[, shortID])
UACR = UACR[matched, ]
colsOut = setdiff(colnames(UACR), cols2keep)
UACR[, get("colsOut") := NULL]
setcolorder(UACR, cols2keep)
setnames(UACR, colnames(UACR), paste0("UACR.", colnames(UACR)))
knownGWAS = cbind(knownGWAS, UACR)

matched = match(knownGWAS[, SNPS], MA[, shortID])
MA = MA[matched, ]
colsOut = setdiff(colnames(MA), cols2keep)
MA[, get("colsOut") := NULL]
setcolorder(MA, cols2keep)
setnames(MA, colnames(MA), paste0("MA.", colnames(MA)))
knownGWAS = cbind(knownGWAS, MA)

#' # Save results 
#(big table)
write.table(knownGWAS, file = "../results/11_Look_Up_GWAS_hits_results.txt", col.names = T, row.names = F)

#smaller table
cols2keep = c("STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "STRONGEST SNP-RISK ALLELE", "SNPS", "CHR_ID", "CHR_POS", "RISK ALLELE FREQUENCY", "P-VALUE",
              "eGFR.N", "eGFR.effect_allele", "eGFR.EAF", "eGFR.infoScore", "eGFR.I2", "eGFR.beta", "eGFR.SE", "eGFR.P", "eGFR.logP", "eGFR.invalid_assoc", 
              "eGFR.reason_for_exclusion", "UA.N", "UA.effect_allele", "UA.EAF", "UA.infoScore", "UA.I2", "UA.beta", "UA.SE", "UA.P", "UA.logP", 
              "UA.invalid_assoc", "UA.reason_for_exclusion")
colsOut = setdiff(colnames(knownGWAS), cols2keep)
knownGWAS[, get("colsOut") := NULL]
write.table(knownGWAS, file = "../results/11_Look_Up_GWAS_hits_eGFR_UA_only.txt", col.names = T, row.names = F)

#input for annotation pipeline
step25 = fread("../results/05_b_step25_Matching_Table_UKBB.txt.gz")

cols2keep = c("STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "STRONGEST SNP-RISK ALLELE", "SNPS", "CHR_ID", "CHR_POS", "RISK ALLELE FREQUENCY", "P-VALUE")
colsOut = setdiff(colnames(knownGWAS), cols2keep)
knownGWAS[, get("colsOut") := NULL]
matched = match(knownGWAS[, SNPS], step25[, UKBB.ID])
knownGWAS[, snp := step25[matched, SNPinRef.x]]
knownGWAS[, position := CHR_POS]
knownGWAS[, chr := 23]
knownGWAS[, CHR_ID := NULL]
knownGWAS[, CHR_POS := NULL]
setcolorder(knownGWAS, c("snp", "chr", "position", "STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "STRONGEST SNP-RISK ALLELE", "SNPS", 
                         "RISK ALLELE FREQUENCY", "P-VALUE"))
write.table(knownGWAS, file = "../../../10_metaGWAS/_results_combined/Lookup_GWAS_Hits/SNPs_to_annotate.txt", col.names = T, row.names = F)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
