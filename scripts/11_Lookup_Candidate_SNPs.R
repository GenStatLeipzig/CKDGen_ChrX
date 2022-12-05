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
phenotypes2keep = c("Creatinine levels", "Estimated glomerular filtration rate", "Glomerular filtration rate", 
                    "Blood urea nitrogen levels", "Serum uric acid levels", "Serum creatinine levels")

knownGWAS = knownGWAS[is.element(`DISEASE/TRAIT`, phenotypes2keep), ]
table(knownGWAS[, `DISEASE/TRAIT`])

#How many entries are there?
dim(knownGWAS)

#How many unique SNPs are there?
length(unique(knownGWAS[, SNPS]))


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

#calculate one-sided p-values for replication check
eGFR[, P := pnorm(-abs(beta/SE))]
eGFR[, logP := -log10(P)]
CKD[, P := pnorm(-abs(beta/SE))]
CKD[, logP := -log10(P)]
BUN[, P := pnorm(-abs(beta/SE))]
BUN[, logP := -log10(P)]
UA[, P := pnorm(-abs(beta/SE))]
UA[, logP := -log10(P)]
Gout[, P := pnorm(-abs(beta/SE))]
Gout[, logP := -log10(P)]
UACR[, P := pnorm(-abs(beta/SE))]
UACR[, logP := -log10(P)]
MA[, P := pnorm(-abs(beta/SE))]
MA[, logP := -log10(P)]


#' # Look up of known GWAS hits in our results
cols2keep = c("STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "REGION" , "CHR_ID", "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE", 
              "STRONGEST SNP-RISK ALLELE", "SNPS", "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG", "OR or BETA")
colsOut = setdiff(colnames(knownGWAS), cols2keep)
knownGWAS[, get("colsOut") := NULL]

cols2keep = c("position", "N", "effect_allele", "EAF", "infoScore", "I2", "beta", "SE", "P", "logP", "invalid_assoc", "reason_for_exclusion")

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

#smaller table with only eGFR, BUN and UA
#plus checks for reproducing hits
cols2keep = c("STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "STRONGEST SNP-RISK ALLELE", "SNPS", "CHR_ID", "CHR_POS", "RISK ALLELE FREQUENCY", "P-VALUE", "OR or BETA", 
              "eGFR.position", "eGFR.N", "eGFR.effect_allele", "eGFR.EAF", "eGFR.infoScore", "eGFR.I2", "eGFR.beta", "eGFR.SE", "eGFR.P", "eGFR.logP", "eGFR.invalid_assoc", 
              "eGFR.reason_for_exclusion", "BUN.N", "BUN.effect_allele", "BUN.EAF", "BUN.infoScore", "BUN.I2", "BUN.beta", "BUN.SE", "BUN.P", "BUN.logP",
              "BUN.invalid_assoc", "BUN.reason_for_exclusion", "UA.N", "UA.effect_allele", "UA.EAF", "UA.infoScore", "UA.I2", "UA.beta", "UA.SE", "UA.P", 
              "UA.logP", "UA.invalid_assoc", "UA.reason_for_exclusion")
colsOut = setdiff(colnames(knownGWAS), cols2keep)
knownGWAS[, get("colsOut") := NULL]

#check alleles
dummy = strsplit(knownGWAS[, `STRONGEST SNP-RISK ALLELE`], split = "-")
knownGWAS[, RISK_ALLELE := sapply(dummy, function(x) return(x[[2]]))]
knownGWAS[, same_allele_eGFR := RISK_ALLELE == eGFR.effect_allele]
knownGWAS[, same_allele_BUN := RISK_ALLELE == BUN.effect_allele]
knownGWAS[, same_allele_UA := RISK_ALLELE == UA.effect_allele]
#-> no allele given where check is false

#do check by hand
knownGWAS[, SNP_replicated := NA]
knownGWAS[, `OR or BETA` := as.character(`OR or BETA`)]
#eGFR, same allele, p < 0.05, 5.72 z score decrease -> same beta direction -> replicated
knownGWAS[1, SNP_replicated := 1]
knownGWAS[1, `OR or BETA` := "5.72 z score decrease"]
#eGFR, same allele, p < 0.05, 9.266 z score decrease -> same beta direction -> replicated
knownGWAS[2, SNP_replicated := 1]
knownGWAS[2, `OR or BETA` := "9.266 z score decrease"]
#creatinine, same allele, p < 0.05, 0.188414 unit increase -> different beta direction -> replicated with eGFR
knownGWAS[4, SNP_replicated := 1]
knownGWAS[4, `OR or BETA` := "0.188414 unit increase"]
#eGFR, same allele, p < 0.05, 5.454 z score increase -> same beta direction -> replicated
knownGWAS[7, SNP_replicated := 1]
knownGWAS[7, `OR or BETA` := "5.454 z score increase"]
#BUN, same allele, p < 0.05, 0.0132 unit increase -> same beta direction -> replicated
knownGWAS[9, SNP_replicated := 1]
knownGWAS[9, `OR or BETA` := "0.0132 unit increase"]
#BUN, same allele, p < 0.05, 0.014 unit increase -> same beta direction -> replicated
knownGWAS[11, SNP_replicated := 1]
knownGWAS[11, `OR or BETA` := "0.014 unit increase"]
#BUN, same allele, p < 0.05, 0.0118 unit decrease -> same beta direction -> replicated
knownGWAS[12, SNP_replicated := 1]
knownGWAS[12, `OR or BETA` := "0.0118 unit decrease"]
#BUN, same allele, p < 0.05, 0.011 unit increase -> same beta direction -> replicated
knownGWAS[13, SNP_replicated := 1]
knownGWAS[13, `OR or BETA` := "0.011 unit increase"]
#BUN, same allele, p < 0.05, 0.013 unit increase -> same beta direction -> replicated
knownGWAS[14, SNP_replicated := 1]
knownGWAS[14, `OR or BETA` := "0.013 unit increase"]
#UA, same allele, p < 0.05, 0.0184248 unit decrease -> same beta direction -> replicated
knownGWAS[15, SNP_replicated := 1]
knownGWAS[15, `OR or BETA` := "0.0184248 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0157643 unit decrease -> different beta direction -> replicated with eGFR
knownGWAS[16, SNP_replicated := 1]
knownGWAS[16, `OR or BETA` := "0.0157643 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0157191 unit decrease -> different beta direction -> replicated with eGFR
knownGWAS[17, SNP_replicated := 1]
knownGWAS[17, `OR or BETA` := "0.0157191 unit decrease"]
#UA, same allele, p < 0.05, 0.0138 unit increase -> same beta direction -> replicated
knownGWAS[18, SNP_replicated := 1]
knownGWAS[18, `OR or BETA` := "0.0138 unit increase"]
#UA, same allele, p < 0.05, 0.0094 unit decrease -> same beta direction -> replicated
knownGWAS[19, SNP_replicated := 1]
knownGWAS[19, `OR or BETA` := "0.0094 unit decrease"]
#UA, same allele, p < 0.05, 0.013 unit decrease -> same beta direction -> replicated
knownGWAS[20, SNP_replicated := 1]
knownGWAS[20, `OR or BETA` := "0.013 unit decrease"]
#UA, same allele, p < 0.05, 0.0094 unit decrease -> same beta direction -> replicated
knownGWAS[21, SNP_replicated := 1]
knownGWAS[21, `OR or BETA` := "0.0094 unit decrease"]
#creatinine, same alleles, p < 0.05, 0.0138 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[23, SNP_replicated := 1]
knownGWAS[23, `OR or BETA` := "0.0138 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0125 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[25, SNP_replicated := 1]
knownGWAS[25, `OR or BETA` := "0.0125 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0128 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[26, SNP_replicated := 1]
knownGWAS[26, `OR or BETA` := "0.0128 unit increase"]
#creatinine, same allele, p < 0.05, 0.0125 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[27, SNP_replicated := 1]
knownGWAS[27, `OR or BETA` := "0.0125 unit increase"]
#creatinine, same allele, p < 0.05, 0.1699 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[28, SNP_replicated := 1]
knownGWAS[28, `OR or BETA` := "0.1699 unit increase"]
#creatinine, same allele, p < 0.05, 0.0111 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[30, SNP_replicated := 1]
knownGWAS[30, `OR or BETA` := "0.0111 unit increase"]
#creatinine, same allele, p < 0.05, 0.014 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[32, SNP_replicated := 1]
knownGWAS[32, `OR or BETA` := "0.014 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0167 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[33, SNP_replicated := 1]
knownGWAS[33, `OR or BETA` := "0.0167 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0138 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[34, SNP_replicated := 1]
knownGWAS[34, `OR or BETA` := "0.0138 unit increase"]
#creatinine, same allele, p < 0.05, 0.0178 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[36, SNP_replicated := 1]
knownGWAS[36, `OR or BETA` := "0.0178 unit decrease"]
#creatinine, same allele, p < 0.05, 0.0102 unit increase -> different beta direction -> replicated by eGFR
knownGWAS[37, SNP_replicated := 1]
knownGWAS[37, `OR or BETA` := "0.0102 unit increase"]
#creatinine, same allele, p < 0.05, 0.0104 unit decrease -> different beta direction -> replicated by eGFR
knownGWAS[38, SNP_replicated := 1]
knownGWAS[38, `OR or BETA` := "0.0104 unit decrease"]
#UA, same allele, p < 0.05, 0.0167 unit increase -> same beta direction -> replicated
knownGWAS[39, SNP_replicated := 1]
knownGWAS[39, `OR or BETA` := "0.0167 unit increase"]
#UA, same alleles, p 0.05, 0.0145 unit decrease -> same beta direction -> replicated
knownGWAS[40, SNP_replicated := 1]
knownGWAS[40, `OR or BETA` := "0.0145 unit decrease"]
#UA, same allele, p < 0.05, 0.0114 unit increase -> same beta direction -> replicated
knownGWAS[41, SNP_replicated := 1]
knownGWAS[41, `OR or BETA` := "0.0114 unit increase"]
#UA, same alleles, p < 0.05, 0.0118 unit increase -> same beta direction -> replicated
knownGWAS[42, SNP_replicated := 1]
knownGWAS[42, `OR or BETA` := "0.0118 unit increase"]
#creatinine, no risk allele given, p < 0.05, 0.02069 unit decrease -> replicated
knownGWAS[44, SNP_replicated := 1]
knownGWAS[44, `OR or BETA` := "0.02069 unit decrease"]
#(e)GFR, no risk allele given, p < 0.05, 	0.02271 unit increase -> replicated
knownGWAS[46, SNP_replicated := 1]
knownGWAS[46, `OR or BETA` := "0.02271 unit increase"]

table(knownGWAS[, SNP_replicated], useNA = "ifany") 

knownGWAS[, c("same_allele_eGFR", "same_allele_BUN", "same_allele_UA", "RISK_ALLELE") := NULL]

write.table(knownGWAS, file = "../results/11_Look_Up_GWAS_hits_eGFR_BUN_UA_only.txt", col.names = T, row.names = F, sep = "|")


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
