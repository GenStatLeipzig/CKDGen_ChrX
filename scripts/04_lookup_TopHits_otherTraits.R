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
myTab = myTab[,c(1,5,6,8)]

dummy = data.table(region = 7,
                   rsID = "rs149995096:100479327:C:T",
                   phenotype = "eGFR_FEMALE",
                   position = 100479327)

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
ShorterTable[,check_BUN := F]
ShorterTable[P_oneSided_BUN<0.05 & sign(beta_eGFR)!=sign(beta_BUN) & region<16,check_BUN := T]
ShorterTable[region>=16,check_BUN := NA]
ShorterTable$check_BUN

ShorterTable[,check_CKD := F]
ShorterTable[P_oneSided_CKD<0.05 & sign(beta_eGFR)!=sign(beta_CKD) & region<16,check_CKD := T]
ShorterTable[region>=16,check_CKD := NA]
ShorterTable$check_CKD

ShorterTable[,check_Gout := F]
ShorterTable[P_oneSided_Gout<0.05 & sign(beta_UA)==sign(beta_Gout) & region>=16,check_Gout := T]
ShorterTable[region<16,check_Gout := NA]
ShorterTable$check_Gout

names(ShorterTable)
x = c(5,12,19,26,33,40)
myNames = names(ShorterTable)[c(1:4,x+6,x,47,x+1,48,x+4,x+2,49)]
colsOut<-setdiff(colnames(ShorterTable),myNames)
ShorterTable[,get("colsOut"):=NULL]
setcolorder(ShorterTable,myNames)

ShorterTable[,P_oneSided_eGFR:=NULL]
ShorterTable[,P_oneSided_UA:=NULL]


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

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

