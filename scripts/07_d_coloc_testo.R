#' ---
#' title: "Co-localization Part 4: Run coloc against Testosterone"
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
#' *Co-localization analysis of testosterone associations*: 
#' 
#' We tested for overlapping causal variants between kidney trait associations and testosterone loci (in males and females).
#' 
#' We used data from Ruth et al (UKBB)
#'
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath,"scripts/"))
source("../helperFunctions/colocFunction_jp.R")

#' # Get Testosterone data ####
#' ***
#' PMID: 32042192
#' 
#' * GCST90012112 (Testosterone in females)
#' * GCST90012113 (Testosterone in males)
#' * GCST90012114 (Testosterone sex-combined)
#' 
TT = fread(paste0(path_Testosterone,"32042192-GCST90012112-EFO_0004908-Build37.f.tsv.gz")) 
TT = TT[chromosome==23,]
TT = TT[effect_allele_frequency>=0.01,]
TT = TT[effect_allele_frequency<=0.99,]
TT[,pval2 := as.numeric(p_value)]
TT[,N := 230454 ]
TT[,phenotype := "testosterone_females" ]
TT_fem = copy(TT)

TT = fread(paste0(path_Testosterone,"32042192-GCST90012113-EFO_0004908-Build37.f.tsv.gz")) 
TT = TT[chromosome==23,]
TT = TT[effect_allele_frequency>=0.01,]
TT = TT[effect_allele_frequency<=0.99,]
TT[,pval2 := as.numeric(p_value)]
TT[,N := 194453 ]
TT[,phenotype := "testosterone_males" ]
TT_mal = copy(TT)

TT = fread(paste0(path_Testosterone,"32042192-GCST90012114-EFO_0004908-Build37.f.tsv.gz")) 
TT = TT[chromosome==23,]
TT = TT[effect_allele_frequency>=0.01,]
TT = TT[effect_allele_frequency<=0.99,]
TT[,pval2 := as.numeric(p_value)]
TT[,N := 425097 ]
TT[,phenotype := "testosterone_all" ]
TT_all = copy(TT)

TT_fem[,dumID := paste(chromosome,base_pair_location,effect_allele,other_allele,sep=":")]
TT_mal[,dumID := paste(chromosome,base_pair_location,effect_allele,other_allele,sep=":")]
TT_all[,dumID := paste(chromosome,base_pair_location,effect_allele,other_allele,sep=":")]

table(duplicated(TT_fem$dumID))
table(duplicated(TT_mal$dumID))
table(duplicated(TT_all$dumID))

table(is.element(TT_mal$dumID,TT_fem$dumID))
table(is.element(TT_all$dumID,TT_fem$dumID))

TT_fem = TT_fem[dumID %in% TT_mal$dumID & dumID %in% TT_all$dumID,]
TT_mal = TT_mal[dumID %in% TT_fem$dumID & dumID %in% TT_all$dumID,]
TT_all = TT_all[dumID %in% TT_fem$dumID & dumID %in% TT_mal$dumID,]

table(TT_fem$dumID == TT_mal$dumID)
table(TT_fem$dumID == TT_all$dumID)

plot(TT_fem$effect_allele_frequency, TT_mal$effect_allele_frequency)
plot(TT_fem$beta, TT_mal$beta)
abline(0,1)
plot(-log10(TT_fem$pval2), -log10(TT_mal$pval2))
abline(0,1)

TT = rbind(TT_fem,TT_mal,TT_all)
save(TT, file = "../temp/07_TestoAssoc.RData")

#' # Match to CKDGen data ####
#' ***
CKDGen_loci = fread("../results/01_Locus_Definitions.txt")

dumTab = foreach(i=1:22)%do%{
  #i=1
  dat = copy(TT)
  myRow = CKDGen_loci[i,]
 
  dat = dat[base_pair_location<myRow$region_end] 
  dat = dat[base_pair_location>myRow$region_start] 
  dat[,region:=myRow$region]
  dat
}
dumTab = rbindlist(dumTab)
TT2 = copy(dumTab)
TT2[,chrPos_b37:= paste("chrX",base_pair_location,sep=":")]
setnames(TT2,"base_pair_location","pos_b37")
setnames(TT2,"standard_error","se")
setnames(TT2,"pval2","pval")
setnames(TT2,"N","n_samples")
TT2[,maf:= effect_allele_frequency]
TT2[effect_allele_frequency>0.5,maf:= 1-effect_allele_frequency]

TT2[pval<5e-8,.N,by=c("region","phenotype")]
plot(TT2[grepl("female",phenotype),effect_allele_frequency], TT2[grepl("_male",phenotype),effect_allele_frequency],
     xlab = "EAF in females",ylab = "EAF in males")
plot(TT2[grepl("female",phenotype),beta], TT2[grepl("_male",phenotype),beta],
     xlab = "beta in females",ylab = "beta in males")
abline(0,1)
plot(TT2[grepl("female",phenotype),-log10(pval)], TT2[grepl("_male",phenotype),-log10(pval)],
     xlab = "-log10(pval) in females",ylab = "-log10(pval) in males")
abline(0,1)

#' # Run coloc ####
#' ***
myPhenos = c("eGFR_ALL","eGFR_MALE","eGFR_FEMALE","UA_ALL","UA_MALE","UA_FEMALE")
source("../helperFunctions/colocFunction_jp.R")

dumTab = foreach(k=1:length(myPhenos))%do%{
  #k=1
  myPheno = myPhenos[k]
  myPheno2 = gsub(myPheno,pattern="\\_.*",replacement = "")
  message("Working on phenotype ",myPheno)
  myFileName = paste0("../data/CKDGen_ChrX_sumStat_",myPheno,".gz")
  
  data_GWAS = fread(myFileName)
  data_GWAS = data_GWAS[invalid_assoc==F,]
  data_GWAS[,chrPos_b37 := paste("chrX",position,sep=":")]
  setnames(data_GWAS,"position","pos_b37")
  setnames(data_GWAS,"SE","se")
  setnames(data_GWAS,"P","pval")
  setnames(data_GWAS,"N","n_samples")
  setnames(data_GWAS,"MAF","maf")
  
  myeQTLs2 = unique(TT2$phenotype)
  
  dumTab2 = foreach(i=c(1:length(myeQTLs2)))%dopar%{
    #i=1
    data_eQTLs = copy(TT2)
    message("Working on trait ",myeQTLs2[i])
    
    data_eQTLs = data_eQTLs[phenotype==myeQTLs2[i],]
    
    dumTab3 = foreach(j=c(1:dim(CKDGen_loci)[1]))%do%{
      #j=1
      myRow = CKDGen_loci[j,]
      
      data_eQTL1 = copy(data_eQTLs)
      data_eQTL1 = data_eQTL1[region == myRow$region]
      dup2 = data_eQTL1[duplicated(pos_b37),]
      data_eQTL1 = data_eQTL1[!is.element(pos_b37,dup2$pos_b37),]
      
      data_GWAS1 = copy(data_GWAS)
      data_GWAS1 = data_GWAS1[pos_b37 %in% data_eQTL1$pos_b37,]
      dups = data_GWAS1[duplicated(pos_b37),]
      data_GWAS1 = data_GWAS1[!is.element(pos_b37,dups$pos_b37),]
      
      data_eQTL1 = data_eQTL1[pos_b37 %in% data_GWAS1$pos_b37,]
      
      setorder(data_eQTL1,pos_b37)
      setorder(data_GWAS1,pos_b37)
      stopifnot(data_eQTL1$chrPos_b37 == data_GWAS1$chrPos_b37)
      
      res = colocFunction_jp(tab1 = data_GWAS1,tab2 = data_eQTL1,trait1 = myPheno,trait2 = myTissue2,
                             locus = moreInfo$cytoband,locus_name = myRow$gene,plotting = F,
                             col_SNPID = "chrPos_b37", col_pos = "pos_b37",
                             col_beta = "beta",col_se = "se",col_P = "pval",
                             col_N = "n_samples",col_MAF = "maf",col_effect = "effect_allele",col_other="other_allele")
      x2<-as.data.table(res)
      x3<-t(x2)
      x4<-as.data.table(x3)
      names(x4)<-names(res)
      x4[,cytoband:=myRow$region]
      x4[,trait1:=myPheno]
      x4[,trait2:=myeQTLs2[i]]
      x4
      
    }
    dumTab3 = rbindlist(dumTab3)
    dumTab3
  }
  ColocTable = rbindlist(dumTab2)
  ColocTable
}

ColocTable = rbindlist(dumTab)

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]

save(ColocTable,file="../results/07_d_coloc_testo.RData")

ColocTable = ColocTable[,c(7,8,9,1:6)]
names(ColocTable)[1] = "locus_NR"

description = data.table(column = names(ColocTable),
                         description = c("number of locus as used throughout the project","CKD trait and sex-setting","testosterone sex-setting",
                                         "Number of SNPs included in co-localization analysis per test",
                                         "Posterior probability for hypothesis 0: neither trait associated",
                                         "Posterior probability for hypothesis 1: only trait 1 associated (CKD trait)",
                                         "Posterior probability for hypothesis 2: only trait 2 associated (Testosterone)",
                                         "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                         "Posterior probability for hypothesis 4: both trait associated, shared signal"))

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/07_d_coloc_testo.xlsx"

WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

