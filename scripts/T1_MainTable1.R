#' ---
#' title: "Get Table 1"
#' subtitle: "CKDGen - Chr X"
#' author: "Janne Pott"
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
#' **Table 1: Genome-wide significant regions with their respective index SNPs of eGFR and UA**
#' 
#' We want to automatically generate relevant columns for Table 1. These are: 
#' 
#' * region, cytoband, index SNP, best setting, position, EA, OA, EAF, region range
#' * beta, SE, p-value, expl. var, p-val sex IA, flag fdr sig IA, gw_sig by setting, 
#' * ckd assoc (only for eGFR)
#' * number of independent variants
#' * colocalization result
#' 
#' It will not include:
#' 
#' * proposed candidate gene
#' * novelty
#' * kidney function
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(paste0(projectpath,"scripts/"))

load("../results/04_lookup_TopHits_allOtherTraits.RData")
tab0 = fread("../results/01_Locus_Definitions.txt")
sexIA_eGFR = fread("../results/02_sex_ia_egfr.txt")
sexIA_UA = fread("../results/02_sex_ia_uric_acid.txt")
IndepSignals = fread("../results/05_b_Cojo_Select_Results.txt")
coloc_sexIA = fread("../results/03_coloc_sex_ia.txt")

#' # Get general information ####
#' ***
#' * region, cytoband, index SNP, best setting, position, EA, OA, EAF, region range
#' 

tab1 = copy(WideTable)
tab1 = tab1[,c(1:4)]
setnames(tab1,"Cytoband","cytoband")
setnames(tab1,"rsID","indexSNP")
setnames(tab1,"topPheno","bestSetting")

dummy = unlist(strsplit(tab1$indexSNP,":"))
tab1[,position := dummy[seq(2,length(dummy),4)]]
tab1[,EA := dummy[seq(4,length(dummy),4)]]
tab1[,OA := dummy[seq(3,length(dummy),4)]]
tab1[16,EA:="A"]
tab1[16,OA:="C"]

regranges = tab0$region_end - tab0$region_start
regranges = c(regranges[1:7],regranges[7:22])
tab1[,regionRanges:= regranges]

#' # Get association statistics ####
#' ***
#' * EAF, beta, SE, p-value, expl. var
#' 
myPhenos = unique(tab1$bestSetting)

tab2 = foreach(i=1:length(myPhenos))%do%{
  #i=1
  myPheno = myPhenos[i]
  tab3= copy(WideTable)
  myNames = names(tab3)[grepl(myPheno,names(tab3))]
  colsOut<-setdiff(colnames(tab3),myNames)
  tab3[,get("colsOut"):=NULL]
  names(tab3) = c("numberOfStudies","N","I2","EAF",
                  "infoScore","beta","SE","P",
                  "P_oneSided","invalid_assoc")  
  
  filt = tab1$bestSetting == myPheno
  tab1[filt,EAF := tab3[filt,EAF]]
  tab1[filt,beta := tab3[filt,beta]]
  tab1[filt,SE := tab3[filt,SE]]
  tab1[filt,pvalue := tab3[filt,P]]
  tab1[filt,N := tab3[filt,N]]
  tab1[filt,explVar := beta^2/(beta^2 + N*SE^2)]
  tab1
}
tab1[,N:=NULL]

#' # Get sexIA and gw by setting ####
#' ***
#' * p-val sex IA, flag fdr sig IA, 
#' * gw_sig by setting, 
#' * ckd assoc (only for eGFR)
#' 
stopifnot(tab1$indexSNP == sexIA_eGFR$SNP)
tab1[,pvalue_sexIA:=sexIA_eGFR$p_diff]
filt = grepl("UA",tab1$bestSetting)
tab1[filt,pvalue_sexIA:=sexIA_UA[filt,p_diff]]
tab1[,pvalue_sexIA_adj := p.adjust(pvalue_sexIA, method = "fdr")]
tab1[pvalue_sexIA_adj<0.05,]  
tab1[pvalue_sexIA<0.05,]
tab1[,sexIA_FDRsig :=F]
tab1[pvalue_sexIA_adj<0.05,sexIA_FDRsig :=T]
tab1[,pvalue_sexIA_adj:=NULL]

stopifnot(tab1$indexSNP == WideTable$rsID)
tab1[,ALL_pval := WideTable[,P_eGFR_ALL]]
tab1[filt,ALL_pval := WideTable[filt,P_UA_ALL]]
tab1[,MALE_pval := WideTable[,P_eGFR_MALE]]
tab1[filt,MALE_pval := WideTable[filt,P_UA_MALE]]
tab1[,FEMALE_pval := WideTable[,P_eGFR_FEMALE]]
tab1[filt,FEMALE_pval := WideTable[filt,P_UA_FEMALE]]

tab1[,ALL_gwsig :=F]
tab1[ALL_pval<5e-8,ALL_gwsig :=T]
tab1[,MALE_gwsig :=F]
tab1[MALE_pval<5e-8,MALE_gwsig :=T]
tab1[,FEMALE_gwsig :=F]
tab1[FEMALE_pval<5e-8,FEMALE_gwsig :=T]

tab1[,ALL_pval :=NULL]
tab1[,MALE_pval :=NULL]
tab1[,FEMALE_pval :=NULL]

stopifnot(tab1$indexSNP == WideTable$rsID)
tab1[,CKD_pval := WideTable[,P_oneSided_CKD_ALL]]
tab1[,CKD_beta := WideTable[,beta_CKD_ALL]]
tab1[filt,CKD_pval := NA]
tab1[filt,CKD_beta := NA]
filt1 = grepl("_MALE",tab1$bestSetting)
filt2 = grepl("_FEMALE",tab1$bestSetting)
tab1[filt1,CKD_pval := WideTable[filt1,P_oneSided_CKD_MALE]]
tab1[filt1,CKD_beta := WideTable[filt1,beta_CKD_MALE]]
tab1[filt2,CKD_pval := WideTable[filt2,P_oneSided_CKD_FEMALE]]
tab1[filt2,CKD_beta := WideTable[filt2,beta_CKD_FEMALE]]
tab1[CKD_pval<0.05 & sign(beta)!=sign(CKD_beta) & region<16,CKD_assoc := T]
tab1[is.na(CKD_assoc) & region<16,CKD_assoc:=F]

tab1[,CKD_pval:=NULL]
tab1[,CKD_beta:=NULL]

#' # Get sexIA and gw by setting ####
#' ***
#' * number of independent variants
#' * colocalization result
IndepSignals[,dumID := paste(region,phenotype,sep="__")]
indep = IndepSignals[,.N,dumID]
indep[,dumID := gsub("_all","_ALL",dumID)]
indep[,dumID := gsub("_male","_MALE",dumID)]
indep[,dumID := gsub("_female","_FEMALE",dumID)]
tab1[,dumID := paste(region,bestSetting,sep="__")]
matched = match(tab1$dumID,indep$dumID)
table(is.na(matched))
tab1[,NR_indepSignals := indep[matched,N]]
tab1[,dumID := NULL]

coloc_sexIA = rbind(coloc_sexIA[1:7,],coloc_sexIA[7:22,])
tab1[,coloc_sexIA:="inconclusive"]
filt = coloc_sexIA$PP.H4.abf>0.75
tab1[filt,coloc_sexIA:="both sexes - shared"]
filt = coloc_sexIA$PP.H3.abf>0.75
tab1[filt,coloc_sexIA:="both sexes - independent"]
filt = coloc_sexIA$PP.H2.abf>0.75
tab1[filt,coloc_sexIA:="female driven"]
filt = coloc_sexIA$PP.H1.abf>0.75
tab1[filt,coloc_sexIA:="male driven"]

tab1

#' # Save ####
#' ***
description = data.table(column = names(tab1),
                         description = c("Number of associated region",
                                         "Genomic cytoband of index SNP",
                                         "SNP with lowest p-value in this region",
                                         "Best phenotype and setting",
                                         "Genomic position according to hg19",
                                         "Effect allele",
                                         "Other allele",
                                         "Range of genomic region (min. +/-500 kb) according to locus collapsing",
                                         "Effect allele frequency in best setting",
                                         "Beta estimate in best setting",
                                         "Standard error in best setting",
                                         "P-value in best setting",
                                         "Explained variance in best setting",
                                         "Unadjusted p-value of sex-interaction test, correcting for correlation with Spearmans rho (0.1604 in eGFR and 0.1285 in UA)",
                                         "TRUE/FALSE flag indicating significant sex interaction after FDR correction",
                                         "TRUE/FALSE flag indicating genome-wide significance in setting ALL for given phenotype",
                                         "TRUE/FALSE flag indicating genome-wide significance in setting MALE for given phenotype",
                                         "TRUE/FALSE flag indicating genome-wide significance in setting FEMALE for given phenotype",
                                         "TRUE/FALSE flag indicating nominal (one-sided) association with CKD in our data with discordant effect direction (only for eGFR)",
                                         "Number independent signals at this loci in best setting",
                                         "Co-localization result for each loci (PP4>0.75, both sexes-shared; PP3>0.75, both sexes-independent; PP2>0.75, female-driven; PP1>0.75, male-driven; no PP>0.75, inconclusive"))

tosave4 = data.table(data = c("tab1","description"), 
                     SheetNames = c("Table1","Description"))
if(dir.exists("../tables/")==F) dir.create("../tables/") 
excel_fn = "../tables/MainTable1.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
