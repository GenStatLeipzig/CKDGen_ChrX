#' ---
#' title: "Get Table 2"
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
#' **Table 2: Overlapping regions of eGFR and UA**
#' 
#' We want to automatically generate relevant columns for Table 2. These are: 
#' 
#' * region NR for eGFR and UA
#' * cytoband
#' * index SNP for eGFR and UA 
#' * best setting for eGFR and UA 
#' * pairwise LD of index SNPs
#' * colocalization result (PP3 & PP4)
#' * interpretation
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(paste0(projectpath,"scripts/"))

load("../results/04_lookup_TopHits_allOtherTraits.RData")
coloc_overlap = fread("../results/08_coloc_overlap.txt")
LD_overlap = fread("../results//08_a_LD_overlapping_regions.txt")

coloc_overlap
coloc_overlap2 = coloc_overlap[c(1:2,4,3,5:8)]
coloc_overlap2

LD_overlap
LD_overlap2 = LD_overlap[c(1:2,4,3,5:8)]
LD_overlap2

#' # Get general information ####
#' ***
#' * region NR for eGFR and UA
#' * cytoband
#' * index SNP for eGFR and UA 
#' * best setting for eGFR and UA 
#' 

tab1 = copy(WideTable)
tab1 = tab1[,c(1:4)]
setnames(tab1,"Cytoband","cytoband")
setnames(tab1,"rsID","indexSNP")
setnames(tab1,"topPheno","bestSetting")

filt1 = is.element(tab1$indexSNP,LD_overlap2$eGFR.SNP)
tab2 = data.table(region_eGFR = tab1[filt1,region])
tab2[,region_UA := LD_overlap2[,UA.region]]
tab2[,regions := paste(region_eGFR, region_UA,sep=", ")]
tab2[,region_eGFR :=NULL]
tab2[,region_UA :=NULL]

tab2[,cytoband := tab1[filt1,cytoband]]
tab2[5,cytoband := "Xq22.1-2"]

tab2[,indexSNP_eGFR := LD_overlap2$eGFR.SNP]
tab2[,indexSNP_UA := LD_overlap2$UA.SNP]

tab2[,bestSetting_eGFR := tab1[filt1,bestSetting]]
matched = match(tab2$indexSNP_UA,tab1$indexSNP)
tab2[,bestSetting_UA := tab1[matched,bestSetting]]
tab2[,bestSetting_eGFR := gsub("eGFR_","",bestSetting_eGFR)]
tab2[,bestSetting_UA := gsub("UA_","",bestSetting_UA)]
tab2[,bestSettings := paste(bestSetting_eGFR, bestSetting_UA,sep=", ")]
tab2[,bestSetting_eGFR :=NULL]
tab2[,bestSetting_UA :=NULL]

tab2[,indexSNP_eGFR := gsub(":.*","",indexSNP_eGFR)]
tab2[,indexSNP_UA := gsub(":.*","",indexSNP_UA)]
tab2[8,indexSNP_eGFR := "chr23:152898260"]
tab2

#' # Get LD and coloc information ####
#' ***
tab2[,LD_r2 := LD_overlap2$r2]

tab2[,PP_H3 := coloc_overlap2$PP.H3.abf]
tab2[,PP_H4 := coloc_overlap2$PP.H4.abf]

tab2

#' # Get interpretation ####
#' ***
#' * overlapping == yes: LD_r2>0.8 or PP_H4>0.75
#' * overlapping == no: LD_r2<0.8 and PP_H3>0.75
#' * overlapping == inconclusive: LD<0.8 and PP_H3<0.75 and PP_H4<0.75
tab2[,overlapping:="inconclusive"]
tab2[LD_r2>0.8 | PP_H4>0.75,overlapping:="yes"]
tab2[LD_r2<0.8 & PP_H3>0.75,overlapping:="no"]
tab2


#' # Save ####
#' ***
description = data.table(column = names(tab2),
                         description = c("Index numbers of overlapping regions",
                                         "Genomic cytoband of regions",
                                         "SNP with lowest p-value in this region for eGFR",
                                         "SNP with lowest p-value in this region for UA",
                                         "Best setting for each index SNP (eGFR, UA)",
                                         "Pairwise LD r^2 between the index SNPs (reference data: UKBB)",
                                         "Posterior Probability of co-localization hypothesis 3 (independent signals)",
                                         "Posterior Probability of co-localization hypothesis 4 (shared signal)",
                                         "Interpretation: yes: index SNPs in LD (r2>0.8) or shared signal (H4>0.75); no: index SNPs not in LD (r2<0.8) and independent signals (H3>0.75); inconclusive: index SNPs not in LD (r2<0.8), but PP inconclusive regarding shared or independent signals"))

tosave4 = data.table(data = c("tab2","description"), 
                     SheetNames = c("Table2","Description"))
if(dir.exists("../tables/")==F) dir.create("../tables/") 
excel_fn = "../tables/MainTable2.xlsx"
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
