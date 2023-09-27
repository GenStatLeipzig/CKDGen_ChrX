#' ---
#' title: "Lookup of eQTL associations"
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
#' I want to check for each lead SNP the association with eQTLs in any tissue. 
#' 
#' Data taken from GTEx v8 and Neptun
#' 
#' # Initialize ####
#' ***
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(projectpath_main)
tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

#' Load data ####
#' ***
loci = fread("../results/01_Locus_Definitions.txt")
myTab = copy(loci)
myTab = myTab[,c(1,5,6,8,17,12,15,16,10,19)]

dummy = data.table(region = 7,
                   rsID = "rs149995096:100479327:C:T",
                   phenotype = "eGFR_FEMALE",
                   position = 100479327,
                   beta = -0.00312952,
                   EAF = 0.18097165571591,
                   effect_allele = "T",
                   other_allele = "C",
                   N = 39188, 
                   P = 2.20758e-09)
dumm2 = data.table(region = 8,
                   rsID = "rs12851072:102925123:T:A",
                   phenotype = "eGFR_ALL",
                   position = 102925123,
                   beta = 0.00142628,
                   EAF = 0.449420598771622,
                   effect_allele = "A",
                   other_allele = "T",
                   N = 773185,
                   P = 5.30351e-08)
dumm4 = data.table(region = 19,
                   rsID = "rs4557889:102359307:A:G",
                   phenotype = "UA_ALL",
                   position = 102359307,
                   beta = 0.0131624,
                   EAF = 0.6640011,
                   effect_allele = "G",
                   other_allele = "A",
                   N = 511971, 
                   P = 1.2735e-07)

myTab = rbind(myTab,dummy,dumm2, dumm4)
setorder(myTab,region,position)
myTab

myeQTLs2<-dir(path = "../temp/07_coloc/",pattern = ".RData")
myeQTLs2 = myeQTLs2[grepl("GTEx_v8_filtered_",myeQTLs2) | grepl("NEPTUNE_filtered_",myeQTLs2)]

dumTab = foreach(i=1:length(myeQTLs2))%do%{
  #i=49
  myeQTL = myeQTLs2[i]
  load(paste0("../temp/07_coloc/",myeQTL))
  data0 = data0[pos_b37 %in% myTab$position,]
  
  matched = match(data0$pos_b37,myTab$position)
  data0[,CKD_EA := myTab[matched,effect_allele]]
  data0[,CKD_OA := myTab[matched,other_allele]]
  data0[,CKD_EAF := myTab[matched,EAF]]
  data0[,CKD_beta := myTab[matched,beta]]
  data0[, CKD_EA2 := CKD_EA]
  data0[effect_allele != CKD_EA, CKD_EA2 := CKD_OA]
  data0[, CKD_beta2 := CKD_beta]
  data0[effect_allele != CKD_EA, CKD_beta2 := (-1)*CKD_beta]
  data0[,sameEffectDir_CKD_eQTL := T]
  data0[sign(beta) != sign(CKD_beta2),sameEffectDir_CKD_eQTL := F]
  #data0[pval<0.05,table(sameEffectDir_CKD_eQTL)]
  data0[,CKD_pheno := myTab[matched,phenotype]]
  data0[,CKD_N := myTab[matched,N]]
  data0[,CKD_pval := myTab[matched,P]]
  
  data0
}
eQTL_lookup = rbindlist(dumTab, fill=TRUE)
eQTL_lookup[pval<0.05, table(sameEffectDir_CKD_eQTL)]
save(eQTL_lookup,file=paste0("../temp/13_lookup_eQTL_leadSNPs_",tag,".RData"))

#' # Filter data ####
#' *** 
#' I only want this for the genes which were colocalized
#' 
load("../results/07_b_coloc_eQTLs.RData")
ColocTable = ColocTable[PP.H4.abf>=0.75]
genes = ColocTable[,unique(gene)]
genes

ColocTable[,tissue := gsub("GE in ","",trait2)]
ColocTable[,dumID := paste(gene,tissue,sep="_")]

eQTL_lookup2 = copy(eQTL_lookup)
eQTL_lookup2 = eQTL_lookup2[gene %in% genes]
eQTL_lookup2[,dumID := paste(gene,tissue,sep="_")]
eQTL_lookup2 = eQTL_lookup2[dumID %in% ColocTable$dumID]

eQTL_lookup2[gene == "CDKL5" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "CDK16" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "USP11" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "ARMCX2" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "TCEAL3" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "MORF4L2" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "ACSL4" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]
eQTL_lookup2[gene == "SLC25A5" & pval<0.05 ,c(1:6,22,7,27,29,30,31)]

#' Summary: 
#' 
#' - CDKL5 (Kidney_Cortex_Tubulointerstitial): FALSE, eGFR_ALL
#' - CDK16 (Stomach): FALSE, eGFR_ALL
#' - USP11 (Thyroid): FALSE, eGFR_ALL
#' - ARMCX2 (Kidney_Cortex_Tubulointerstitial): 
#'    - TRUE, UA_ALL
#'    - FALSE, eGFR_ALL
#' - TCEAL3 (Muscle_Skeletal): 
#'    - TRUE, UA_ALL
#'    - FALSE, eGFR_ALL
#' - MORF4L2 (Whole_Blood):
#'    - FALSE, UA_ALL
#'    - TRUE, eGFR_ALL 
#' - ACSL4 (Whole_Blood): TRUE, eGFR_ALL
#' - SLC25A5 (Kidney_Cortex_Tubulointerstitial): TRUE, eGFR_ALL
#' 
save(eQTL_lookup2,file=paste0("../results/13_lookup_eQTL_leadSNPs_",tag,".RData"))
names(eQTL_lookup2[,c(1:8,12:14,20:22,24:27,29,30,31)])
eQTL_lookup2 = eQTL_lookup2[,c(1:8,12:14,20:22,24:27,29,30,31)]
description = data.table(column = names(eQTL_lookup2),
                         description = c("ENSG ID",
                                         "Distance SNP-gene start",
                                         "eQTL Number of samples with minor allele",
                                         "eQTL Number of minor alleles",
                                         "eQTL Minor allele frequency",
                                         "eQTL pvalue","eQTL beta","eQTL SE","eQTL OA","eQTL EA","eQTL number of samples",
                                         "Gene","Cytoband","Tissue",
                                         "CKD EA", "CKD OA", "CKD EAF","CKD beta","CKD beta (allele switch to map to eQTLS)",
                                         "Same effect direction (after correcting the CKD beta for same effect allele)?",
                                         "CKD Phenotype"))

tosave4 = data.table(data = c("eQTL_lookup2", "description"), 
                     SheetNames = c("eQTL_lookup", "Description"))
excel_fn = paste0("../results/13_lookup_eQTL_effect_directions_",tag,".xlsx")

WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
