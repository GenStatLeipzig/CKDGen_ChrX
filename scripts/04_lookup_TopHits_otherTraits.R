#' ---
#' title: "Comparison of eGFR/UA hits with CKD/BUN/Gout"
#' subtitle: "CKDGen - Chr X"
#' author: "Janne Pott, Katrin Horn"
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
#' This is a script to look up the top hits of eGFR and UA in the other kidney related traits (CKD, BUN, Gout).
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Loading loci ####
#' ***
loci = fread("../results/01_Locus_Definitions.txt")
loci

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

#' # Create ToDo List ####
#' ***
statistics = list.files(path = "../data/",pattern = ".gz", recursive = TRUE)
pheno = gsub("CKDGen_ChrX_sumStat_","",statistics)
pheno = gsub(".gz","",pheno)
dummy = unlist(strsplit(pheno,"_"))

ToDoList = data.table(phenotype = pheno,
                      trait = dummy[seq(1,length(dummy),2)],
                      setting = dummy[seq(2,length(dummy),2)],
                      statistics = statistics,
                      statistic_path = paste0("../data/",statistics))
ToDoList

#' # Loop to extract data ####
#' ***
dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Woriking on ",myRow$phenotype)
  
  GWAMA = fread(myRow$statistic_path)
  matched = match(myTab[, rsID], GWAMA[, rsID])
  sum(is.na(matched))   #should be 0
  GWAMA = GWAMA[matched,]
  myNames = c("rsID","phenotype","Cytoband","numberOfStudies","N","EAF","beta","SE","P","I2","infoScore","invalid_assoc")
  colsOut<-setdiff(colnames(GWAMA),myNames)
  GWAMA[,get("colsOut"):=NULL]
  GWAMA[,region := myTab[,region]]
  GWAMA[,topPheno := myTab[,phenotype]]
  if(myRow$trait %nin% c("eGFR","UA")){
    GWAMA[,P_oneSided := pnorm(-abs(beta/SE))]
    GWAMA[,beta_toComp := myTab[,beta]]
    GWAMA[region<16 & sign(beta)==sign(beta_toComp),P_oneSided := 1-P_oneSided]
    GWAMA[region>=16 & sign(beta)!=sign(beta_toComp),P_oneSided := 1-P_oneSided]
    GWAMA[,beta_toComp := NULL]
    
  }
  GWAMA
  
}
dumTab = rbindlist(dumTab,fill = T,use.names = T)
dumTab
dummy = copy(dumTab)[1:23,]
dummy = dummy[,c(13,2,1,14)]

#' Get wide table with all results
dumTab2<-dcast(dumTab,
               formula = rsID ~ phenotype,
               value.var = c("numberOfStudies","N","I2","EAF","infoScore","beta","SE","P","P_oneSided","invalid_assoc"),
               sep = "_")

matched = match(dummy$rsID,dumTab2$rsID)
WideTable = cbind(dummy,dumTab2[matched,-1])

#' Get short table with matchin top settings 
dumTab_ALL = copy(dumTab)
dumTab_ALL = dumTab_ALL[grepl("_ALL",topPheno),]
dumTab_ALL = dumTab_ALL[grepl("_ALL",phenotype),]
dumTab_ALL[,phenotype := gsub("_ALL","",phenotype)]
dumTab_FEM = copy(dumTab)
dumTab_FEM = dumTab_FEM[grepl("_FEM",topPheno),]
dumTab_FEM = dumTab_FEM[grepl("_FEM",phenotype),]
dumTab_FEM[,phenotype := gsub("_FEMALE","",phenotype)]
dumTab_MAL = copy(dumTab)
dumTab_MAL = dumTab_MAL[grepl("_MAL",topPheno),]
dumTab_MAL = dumTab_MAL[grepl("_MAL",phenotype),]
dumTab_MAL[,phenotype := gsub("_MALE","",phenotype)]

dumTab3_ALL<-dcast(dumTab_ALL,
               formula = rsID ~ phenotype,
               value.var = c("N","EAF","beta","SE","P","P_oneSided"),
               sep = "_")
dumTab3_FEM<-dcast(dumTab_FEM,
                   formula = rsID ~ phenotype,
                   value.var = c("N","EAF","beta","SE","P","P_oneSided"),
                   sep = "_")
dumTab3_MAL<-dcast(dumTab_MAL,
                   formula = rsID ~ phenotype,
                   value.var = c("N","EAF","beta","SE","P","P_oneSided"),
                   sep = "_")

dumTab3 = rbind(dumTab3_ALL,dumTab3_FEM,dumTab3_MAL,fill=TRUE)
matched = match(dummy$rsID,dumTab3$rsID)
ShorterTable = cbind(dummy,dumTab3[matched,-1])

#' check significance and sign
ShorterTable[,check_BUN_eGFR := F]
ShorterTable[P_oneSided_BUN<0.05 & sign(beta_eGFR)!=sign(beta_BUN),check_BUN_eGFR := T]
ShorterTable[,table(check_BUN_eGFR)]
ShorterTable[,check_BUN_UA := F]
ShorterTable[P_oneSided_BUN<0.05 & sign(beta_UA)==sign(beta_BUN),check_BUN_UA := T]
ShorterTable[,table(check_BUN_UA)]

ShorterTable[,check_CKD_eGFR := F]
ShorterTable[P_oneSided_CKD<0.05 & sign(beta_eGFR)!=sign(beta_CKD),check_CKD_eGFR := T]
ShorterTable[,table(check_CKD_eGFR)]
ShorterTable[,check_CKD_UA := F]
ShorterTable[P_oneSided_CKD<0.05 & sign(beta_UA)==sign(beta_CKD),check_CKD_UA := T]
ShorterTable[,table(check_CKD_UA)]

ShorterTable[,check_UACR_eGFR := F]
ShorterTable[P_oneSided_UACR<0.05 & sign(beta_eGFR)!=sign(beta_UACR),check_UACR_eGFR := T]
ShorterTable[,table(check_UACR_eGFR)]
ShorterTable[,check_UACR_UA := F]
ShorterTable[P_oneSided_UACR<0.05 & sign(beta_UA)==sign(beta_UACR),check_UACR_UA := T]
ShorterTable[,table(check_UACR_UA)]

ShorterTable[,check_MA_eGFR := F]
ShorterTable[P_oneSided_MA<0.05 & sign(beta_eGFR)!=sign(beta_MA),check_MA_eGFR := T]
ShorterTable[,table(check_MA_eGFR)]
ShorterTable[,check_MA_UA := F]
ShorterTable[P_oneSided_MA<0.05 & sign(beta_UA)==sign(beta_MA),check_MA_UA := T]
ShorterTable[,table(check_MA_UA)]

ShorterTable[,check_Gout_eGFR := F]
ShorterTable[P_oneSided_Gout<0.05 & sign(beta_eGFR)!=sign(beta_Gout),check_Gout_eGFR := T]
ShorterTable[,table(check_Gout_eGFR)]
ShorterTable[,check_Gout_UA := F]
ShorterTable[P_oneSided_Gout<0.05 & sign(beta_UA)==sign(beta_Gout),check_Gout_UA := T]
ShorterTable[,table(check_Gout_UA)]

names(ShorterTable)
#' eGFR - UA - BUN - CKD - UACR - MA - Gout
x = c(5,12,19,26,33,40)

myNames = names(ShorterTable)[c(1:4,
                                x+6,       # eGFR
                                x+4,       # UA
                                x,47,48,   # BUN
                                x+1,49,50, # CKD
                                x+5,51,52, # UACR
                                x+3,53,54, # MA
                                x+2,55,56)]# Gout
colsOut<-setdiff(colnames(ShorterTable),myNames)
colsOut
setcolorder(ShorterTable,myNames)

ShorterTable[,P_oneSided_eGFR:=NULL]
ShorterTable[,P_oneSided_UA:=NULL]

#' Get Markus short table with all settings
#'
dumTab_ALL = copy(dumTab)
dumTab_ALL = dumTab_ALL[grepl("_ALL",phenotype),]
dumTab_FEM = copy(dumTab)
dumTab_FEM = dumTab_FEM[grepl("_FEM",phenotype),]
dumTab_MAL = copy(dumTab)
dumTab_MAL = dumTab_MAL[grepl("_MAL",phenotype),]

dumTab3_ALL<-dcast(dumTab_ALL,
                   formula = rsID ~ phenotype,
                   value.var = c("N","EAF","beta","SE","P","P_oneSided"),
                   sep = "_")
dumTab3_FEM<-dcast(dumTab_FEM,
                   formula = rsID ~ phenotype,
                   value.var = c("N","EAF","beta","SE","P","P_oneSided"),
                   sep = "_")
dumTab3_MAL<-dcast(dumTab_MAL,
                   formula = rsID ~ phenotype,
                   value.var = c("N","EAF","beta","SE","P","P_oneSided"),
                   sep = "_")

matched_all = match(dummy$rsID,dumTab3_ALL$rsID)
matched_mal = match(dummy$rsID,dumTab3_MAL$rsID)
matched_fem = match(dummy$rsID,dumTab3_FEM$rsID)

ShorterTable_v2 = cbind(dummy,
                        dumTab3_ALL[matched_all,c(36,35,37,38,41,40,39)],
                        dumTab3_MAL[matched_mal,c(36,35,37,38,41,40,39)],
                        dumTab3_FEM[matched_fem,c(36,35,37,38,41,40,39)])
x = c(5,12,19)

myNames = names(ShorterTable_v2)[c(1:4,
                                x,       # eGFR
                                x+1,       # UA
                                x+2,   # BUN
                                x+3, # CKD
                                x+4, # UACR
                                x+5, # MA
                                x+6)]# Gout
colsOut<-setdiff(colnames(ShorterTable_v2),myNames)
colsOut
setcolorder(ShorterTable_v2,myNames)

names(ShorterTable_v2) = gsub("P_oneSided_","",names(ShorterTable_v2))
names(ShorterTable_v2) = gsub("P_","",names(ShorterTable_v2))

ShorterTable_v2[,5:25] =  -log10(ShorterTable_v2[,5:25])

#' # Save data ####
#' ***
#' Both as text files and one excel file with two sheets

WriteXLS(x = c("WideTable","ShorterTable"), 
         ExcelFileName = "../results/04_lookup_TopHits_otherTraits.xlsx",
         SheetNames = c("AllSettings","MatchingSettings"),
         AutoFilter=T, 
         BoldHeaderRow=T)
save(WideTable, file = "../results/04_lookup_TopHits_allOtherTraits.RData")
save(ShorterTable, file = "../results/04_lookup_TopHits_matchingSettings.RData")
save(ShorterTable_v2, file = "../results/04_lookup_TopHits_allSettings_logP.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

