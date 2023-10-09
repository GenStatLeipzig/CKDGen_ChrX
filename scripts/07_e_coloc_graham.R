#' ---
#' title: "Co-localization Part 5: Coloc with Graham data"
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
#' *Co-localization analysis of Graham loci*:
#' 
#' One feedback comment from Anna Koetgen was that the locus number 9 was already reported and is not novel as described. 
#' 
#' Here, we test if it shares indeed a signal with the SumStats provided by [Graham et al.](https://csg.sph.umich.edu/willer/public/eGFR2018/). First we test eGFR_ALL, eGFR_MALE, and eGFR_FEMALE, and then we test the two conditional statistics for eGFR_ALL. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Load data ####
#' ***
#' ## Graham data (downloaded 31.01.2023)  ####
erg_Graham = fread("../../../06_lit/SumStats_Graham_2018/2018_Graham_et_al_eGFR_meta.tbl")
erg_Graham = erg_Graham[grepl("X:",MarkerName)]

setnames(erg_Graham,"Alt","effect_allele")
setnames(erg_Graham,"Ref","other_allele")
setnames(erg_Graham,"Freq","EAF")
setnames(erg_Graham,"Weight","nSamples")

erg_Graham[,other_allele := toupper(other_allele)]
erg_Graham[,effect_allele := toupper(effect_allele)]

erg_Graham[,pos := gsub("_.*","",MarkerName)]
erg_Graham[,pos := gsub("X:","",pos)]
erg_Graham[,pos := as.numeric(pos)]

erg_Graham[,dumID1 := paste(pos,effect_allele,other_allele,sep=":")]
erg_Graham[,dumID2 := paste(pos,other_allele,effect_allele,sep=":")]
erg_Graham[,dumID := paste0("X:",pos)]

erg_Graham[,MAF := EAF]
erg_Graham[EAF>0.5,MAF := 1-EAF]

erg_Graham

#' ## CLDGen data (our GWAMA)  ####
erg1_A = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
erg1_F = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
erg1_M = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")

erg1_A = erg1_A[invalid_assoc == F,]
erg1_M = erg1_M[invalid_assoc == F,]
erg1_F = erg1_F[invalid_assoc == F,]

erg1_A = erg1_A[position %in% erg_Graham$pos]
erg1_M = erg1_M[position %in% erg_Graham$pos]
erg1_F = erg1_F[position %in% erg_Graham$pos]

erg1_A[,dumID1 := paste(position,effect_allele,other_allele,sep=":")]
erg1_M[,dumID1 := paste(position,effect_allele,other_allele,sep=":")]
erg1_F[,dumID1 := paste(position,effect_allele,other_allele,sep=":")]
erg1_A[,dumID := paste("X",position,sep=":")]
erg1_M[,dumID := paste("X",position,sep=":")]
erg1_F[,dumID := paste("X",position,sep=":")]

filt1 = is.element(erg1_A$dumID1,erg_Graham$dumID1) | is.element(erg1_A$dumID1,erg_Graham$dumID2)
erg1_A = erg1_A[filt1,]
filt2 = is.element(erg1_M$dumID1,erg_Graham$dumID1) | is.element(erg1_M$dumID1,erg_Graham$dumID2)
erg1_M = erg1_M[filt2,]
filt3 = is.element(erg1_F$dumID1,erg_Graham$dumID1) | is.element(erg1_F$dumID1,erg_Graham$dumID2)
erg1_F = erg1_F[filt3,]

dupPos = erg1_A[duplicated(position),position]
erg1_A = erg1_A[!is.element(position,dupPos),]
dupPos = erg1_M[duplicated(position),position]
erg1_M = erg1_M[!is.element(position,dupPos),]
dupPos = erg1_F[duplicated(position),position]
erg1_F = erg1_F[!is.element(position,dupPos),]

allPos = c(erg1_A$position,erg1_F$position,erg1_M$position)
allPos = unique(allPos)
erg_Graham = erg_Graham[pos %in% allPos]

#' # Coloc for all loci ####
#' ***
myTab = fread("../results/01_Locus_Definitions.txt")

dumTab = foreach(i=1:dim(myTab)[1])%do%{
  #i=1
  myRow = myTab[i,]
  
  graham = copy(erg_Graham)
  graham = graham[pos < myRow$region_end & pos >myRow$region_start]
  
  CKDGenA = copy(erg1_A)
  CKDGenA = CKDGenA[position < myRow$region_end & position >myRow$region_start]
  CKDGenM = copy(erg1_M)
  CKDGenM = CKDGenM[position < myRow$region_end & position >myRow$region_start]
  CKDGenF = copy(erg1_F)
  CKDGenF = CKDGenF[position < myRow$region_end & position >myRow$region_start]
  
  setorder(CKDGenA,position)
  setorder(CKDGenM,position)
  setorder(CKDGenF,position)
  setorder(graham,pos)
  
  grahamA = copy(graham)
  grahamA = grahamA[dumID1 %in% CKDGenA$dumID1 | dumID2 %in% CKDGenA$dumID1]
  grahamM = copy(graham)
  grahamM = grahamM[dumID1 %in% CKDGenM$dumID1 | dumID2 %in% CKDGenM$dumID1]
  grahamF = copy(graham)
  grahamF = grahamF[dumID1 %in% CKDGenF$dumID1 | dumID2 %in% CKDGenF$dumID1]
  
  stopifnot(CKDGenA$dumID == grahamA$dumID)
  stopifnot(CKDGenM$dumID == grahamM$dumID)
  stopifnot(CKDGenF$dumID == grahamF$dumID)
  
  
  my_res_A<- coloc::coloc.abf(dataset1=list(beta=CKDGenA$beta,
                                          varbeta=(CKDGenA$SE)^2,
                                          N=CKDGenA$N,
                                          snp=CKDGenA$dumID,
                                          MAF=CKDGenA$MAF,
                                          type="quant"),
                            dataset2=list(pvalues=grahamA$`P-value`,
                                          N=grahamA$nSamples ,
                                          snp=grahamA$dumID,
                                          MAF=grahamA$MAF,
                                          type="quant"))
 
  my_res_M<- coloc::coloc.abf(dataset1=list(beta=CKDGenM$beta,
                                            varbeta=(CKDGenM$SE)^2,
                                            N=CKDGenM$N,
                                            snp=CKDGenM$dumID,
                                            MAF=CKDGenM$MAF,
                                            type="quant"),
                              dataset2=list(pvalues=grahamM$`P-value`,
                                            N=grahamM$nSamples ,
                                            snp=grahamM$dumID,
                                            MAF=grahamM$MAF,
                                            type="quant"))
  
  my_res_F<- coloc::coloc.abf(dataset1=list(beta=CKDGenF$beta,
                                            varbeta=(CKDGenF$SE)^2,
                                            N=CKDGenF$N,
                                            snp=CKDGenF$dumID,
                                            MAF=CKDGenF$MAF,
                                            type="quant"),
                              dataset2=list(pvalues=grahamF$`P-value`,
                                            N=grahamF$nSamples ,
                                            snp=grahamF$dumID,
                                            MAF=grahamF$MAF,
                                            type="quant"))
  
  my_res2A<-my_res_A[[1]]
  my_res2M<-my_res_M[[1]]
  my_res2F<-my_res_F[[1]]
  
  my_res2 = rbind(my_res2A,my_res2M,my_res2F)
  x2<-as.data.table(my_res2)
  x2[,locusNR:=myRow$region]
  x2[,trait1:=c("CKDGen_eGFR_ALL","CKDGen_eGFR_MALE","CKDGen_eGFR_FEMALE")]
  x2[,trait2:="Graham_eGFR_ALL"]
  x2
}
dumTab = rbindlist(dumTab)
dumTab[PP.H4.abf>0.75]
dumTab[PP.H3.abf>0.75]
dumTab[PP.H1.abf>0.75]

#' # Coloc with conditional values ####
#' ***
cond1 = fread("../temp/05_c_Cojo_cond_results/CojoCond_eGFR_all_Region_9_rs111410539.cma.cojo")
cond2 = fread("../temp/05_c_Cojo_cond_results/CojoCond_eGFR_all_Region_9_rs181497961.cma.cojo")

cond1[,dumID := paste0("X:",bp)]
cond2[,dumID := paste0("X:",bp)]
cond1[,dumID1 := paste0(bp,":",refA)]
cond2[,dumID1 := paste0(bp,":",refA)]

cond1[,MAF := freq]
cond1[freq>0.5,MAF := 1-freq]
cond2[,MAF := freq]
cond2[freq>0.5,MAF := 1-freq]

cond1 = cond1[!is.na(bC),]
cond2 = cond2[!is.na(bC),]

graham1 = copy(erg_Graham)
graham1[,dumID1 := paste0(pos,":",effect_allele)]
graham1[,dumID2 := paste0(pos,":",other_allele)]
graham1 = graham1[dumID1 %in% cond1$dumID1 | dumID2 %in% cond1$dumID1]

graham2 = copy(erg_Graham)
graham2[,dumID1 := paste0(pos,":",effect_allele)]
graham2[,dumID2 := paste0(pos,":",other_allele)]
graham2 = graham2[dumID1 %in% cond2$dumID1 | dumID2 %in% cond2$dumID1]

my_res_cond1<- coloc::coloc.abf(dataset1=list(beta=cond1$bC,
                                        varbeta=(cond1$bC_se)^2,
                                        N=cond1$n,
                                        snp=cond1$dumID,
                                        MAF=cond1$MAF,
                                        type="quant"),
                          dataset2=list(pvalues=graham1$P,
                                        N=graham1$nSamples ,
                                        snp=graham1$dumID,
                                        MAF=graham1$MAF,
                                        type="quant"))

my_res_cond2<- coloc::coloc.abf(dataset1=list(beta=cond2$bC,
                                              varbeta=(cond2$bC_se)^2,
                                              N=cond2$n,
                                              snp=cond2$dumID,
                                              MAF=cond2$MAF,
                                              type="quant"),
                                dataset2=list(pvalues=graham2$P,
                                              N=graham2$nSamples ,
                                              snp=graham2$dumID,
                                              MAF=graham2$MAF,
                                              type="quant"))

my_res21<-my_res_cond1[[1]]
my_res22<-my_res_cond2[[1]]

my_res2 = rbind(my_res21,my_res22)
x2<-as.data.table(my_res2)
x2[,locusNR:=9]
x2[,trait1:=c("CKDGen_eGFR_ALL_rs111410539","CKDGen_eGFR_ALL_rs181497961")]
x2[,trait2:="Graham_eGFR_ALL"]
x2

dumTab = rbind(dumTab,x2)
dumTab[locusNR == 9,]

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


