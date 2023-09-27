#' ---
#' title: "Lookup of testosterone associations"
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
#' I want to check for each lead SNP the association with testosterone in males.
#' 
#' Data taken from [Ruth et al.](https://www.nature.com/articles/s41591-020-0751-5)
#' 
#' # Initialize ####
#' ***
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(projectpath_main)

#' Load data ####
#' ***
testo = fread("../../../../2204_meta_steroidhormones/_literatur/UKBB/32042192-GCST90012113-EFO_0004908-Build37.f.tsv.gz")
testo = testo[chromosome == 23,]

data1 = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")

loci = fread("../results/01_Locus_Definitions.txt")
myTab = copy(loci)
myTab = myTab[,c(1,5,6,8,17)]

dummy = data.table(region = 7,
                   rsID = "rs149995096:100479327:C:T",
                   phenotype = "eGFR_FEMALE",
                   position = 100479327,
                   beta = -0.00312952)

myTab = rbind(myTab,dummy)
setorder(myTab,region,position)
myTab

data1 = data1[rsID %in% myTab$rsID,]
testo = testo[variant_id %in% data1$ID_UKBB]

matched = match(testo$variant_id,data1$ID_UKBB)
table(testo$variant_id == data1$ID_UKBB[matched])
testo[,matchingID := data1[matched,rsID]]

matched = match(myTab$rsID,testo$matchingID)
myTab[,beta_testo := testo[matched,beta]]
myTab[,SE_testo := testo[matched,standard_error]]
myTab[,pval_testo := testo[matched,p_value]]
myTab[,EAF_testo := testo[matched,effect_allele_frequency]]
myTab[,EA_testo := testo[matched,effect_allele]]

matched = match(myTab$rsID,data1$rsID)
myTab[,EAF_eGFR := data1[matched,EAF]]
myTab[,EA_eGFR := data1[matched,effect_allele]]


myTab[region %in% c(1,4,16),]

fwrite(myTab, file = "../results/13_Lookup_testo_leadSNPs.txt")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
