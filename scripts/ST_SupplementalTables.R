#' ---
#' title: "Get Supplemental Tables"
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
#' **Supplemental Tables**
#' 
#' 1) Description of Studies (as received from participating studies) **--> not included here**
#' 2) Genotyping & Imputation of Studies (as received from participating studies) **--> not included here**
#' 3) Sample Sizes & SNP Numbers, and inflation factor $\lambda$ per phenotype
#' 4) Comparisons between the sexes (interaction and co-localization --> see script 02 and 03)
#' 5) Cross-phenotype comparision (--> see script 04)
#' 6) Genome-wide significant & independent hit per region (--> see script 05)
#' 7) Annotation of credible sets (a-e for eGFR and UA in their respective settings)
#' 8) Co-localization with eQTLs (--> see script 07)
#' 9) Replication in *HUNT* (--> see script 09)
#' 10) Replication of *Graham et al*, *Kanai et al* and *Sakaue et al* results (--> see script 11)
#' 11) Summary of MR-Mega results (--> see scripts 12)**--> not yet included**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

#' # Get content table (tab0) ####
#' ***
tab0 = data.table(Table = paste0("S",c(1:11)),
                  Title = c("Study descriptions",
                            "Genotyping & imputation information per study",
                            "Sample sizes, SNP Numbers and inflation factor per phenotype and setting",
                            "Within-phenotype comparisons of eGFR and UA associations between the sexes",
                            "Cross-phenotype comparisions",
                            "List of all independent signals per locus",
                            "Annotation of credible sets (sub table for each phenotype and setting)",
                            "Co-localization with eQTLs",
                            "Replication of eGFR ALL hits in HUNT",
                            "Replication of known associations provided in GWAS Catalog",
                            "Summary of MR-Mega results"),
                  Source = c("done manually",
                             "done manually",
                             "done within script ST_SupplementalTables.R",
                             "../results/02_ and ../results/03_",
                             "../results/04_lookup_TopHits_matchingSettings.RData",
                             "../results/05_b_Cojo_Select_Results.txt",
                             "../../../10_metaGWAS/PHENOTYPE/08_credSets/gwasresults_V5/synopsis/topliste_tabdelim/topliste_2022-11-04_credSets.txt",
                             "../results/07_",
                             "../results/09_replication_HUNT.RData",
                             "../results/11_Look_Up_GWAS_hits_eGFR_UA_only.txt",
                             "not yet finished"))

tab0

#' # Get Sup Tab 3 ####
#' ***
#' Sample Sizes & SNP Numbers, and inflation factor $\lambda$ per phenotype
#' 
#' Number of cases & controls for CKD and BUN will be added later (extract information out of autosomal publications)
#' 
#' Please note: I load the initial meta-analyses files. These files will not be shared later. We only provide statistics for SNPs valid in at least one phenotype (e.g. valid in UA but not eGFR --> listed). SNPs that are invalid in all phenotypes are not listed in the data files. However, we want to report with how many SNPs we started, and how many remained after filtering. A total of 325,770 SNPs had at least on valid association. 
#' 

statistics = list.files(path = "../data/",pattern = ".gz", recursive = TRUE)
pheno = gsub("CKDGen_ChrX_sumStat_","",statistics)
pheno = gsub(".gz","",pheno)
dummy = unlist(strsplit(pheno,"_"))

oldGWAMA = list.files(path = "../../../10_metaGWAS/",pattern = "GWASMA_", recursive = TRUE)
oldGWAMA = oldGWAMA[grepl(".RData",oldGWAMA)]
oldGWAMA = oldGWAMA[!grepl("archiv",oldGWAMA)]
oldGWAMA = oldGWAMA[!grepl("onlyEA",oldGWAMA)]
oldGWAMA

ToDoList = data.table(phenotype = pheno,
                      trait = dummy[seq(1,length(dummy),2)],
                      setting = dummy[seq(2,length(dummy),2)],
                      statistics = statistics,
                      statistic_path = paste0("../data/",statistics))
ToDoList

ToDoList2 = data.table(oldGWAMA = oldGWAMA)
ToDoList2[,setting := "ALL"]
ToDoList2[grepl("ChrX_M",oldGWAMA),setting := "MALE"]
ToDoList2[grepl("ChrX_F",oldGWAMA),setting := "FEMALE"]

ToDoList2[,pheno := gsub(".*GWASMA_","",oldGWAMA)]
ToDoList2[,pheno := gsub("_overall.*","",pheno)]
ToDoList2[,pheno := gsub("_ChrX.*","",pheno)]
ToDoList2[,pheno := gsub("_202.*","",pheno)]
ToDoList2[,pheno := gsub("uric_acid*","UA",pheno)]

ToDoList2[,phenotype := paste(pheno,setting,sep="_")]
table(is.element(ToDoList$phenotype,ToDoList2$phenotype))
matched = match(ToDoList$phenotype,ToDoList2$phenotype)
ToDoList[,oldGWAMA:=paste0("../../../10_metaGWAS/",ToDoList2[matched,oldGWAMA])]

tab3 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on phenotype ",myRow$phenotype," (",i," of ",dim(ToDoList)[1],")")
  
  GWAMA = fread(myRow$statistic_path)
  loaded = load(myRow$oldGWAMA)
  oldGWAMA = get(loaded)
  
  GWAMA = GWAMA[!is.na(beta),]
  oldGWAMA = oldGWAMA[!is.na(betaFEM),]
  n_SNPs_unfiltered = dim(oldGWAMA)[1]
  lambda_unfiltered = oldGWAMA[,median((betaFEM/seFEM)^2, na.rm = TRUE)/0.456]
  n_studies_unfiltered = max(oldGWAMA$numberStudies)
  n_samples_unfiltered = max(oldGWAMA$totalN)
  
  GWAMA = GWAMA[invalid_assoc==F,]
  if(dim(GWAMA)[1]>0){
    n_studies_filtered = max(GWAMA$numberOfStudies)
    n_samples_filtered = max(GWAMA$N)
    n_SNPs_filtered = dim(GWAMA)[1]
    lambda_filtered = GWAMA[,median((beta/SE)^2, na.rm = TRUE)/0.456]
  }else{
    n_studies_filtered = NA
    n_samples_filtered = NA
    n_SNPs_filtered = 0
    lambda_filtered = NA
  }


  res = data.table(phenotype = myRow$phenotype,
                   n_studies_unfiltered = n_studies_unfiltered,
                   n_samples_unfiltered = trunc(n_samples_unfiltered),
                   n_SNPs_unfiltered = n_SNPs_unfiltered,
                   lambda_unfiltered = round(lambda_unfiltered,3),
                   n_studies_filtered = n_studies_filtered,
                   n_samples_filtered = trunc(n_samples_filtered),
                   n_SNPs_filtered = n_SNPs_filtered,
                   lambda_filtered = round(lambda_filtered,3))
  res
}
tab3 = rbindlist(tab3)
tab3[,table(n_studies_unfiltered == n_studies_filtered)]
tab3[,table(n_samples_unfiltered == n_samples_filtered)]
tab3[,n_studies_filtered:=NULL]
tab3[,n_samples_filtered:=NULL]
setnames(tab3,"n_studies_unfiltered","n_studies")
setnames(tab3,"n_samples_unfiltered","n_samples")
tab3

#' # Get Sup Tab 4 ####
#' ***
#' Comparisons between the sexes (interaction and co-localization --> see script 02 and 03)
#' 
load("../results/04_lookup_TopHits_allOtherTraits.RData")
sexIA_eGFR = fread("../results/02_sex_ia_egfr.txt")
sexIA_UA = fread("../results/02_sex_ia_uric_acid.txt")
coloc_sexIA = fread("../results/03_coloc_sex_ia.txt")

tab4 = copy(WideTable)
tab4 = tab4[,c(1:4)]
setnames(tab4,"Cytoband","cytoband")
setnames(tab4,"rsID","indexSNP")
setnames(tab4,"topPheno","bestSetting")
tab4

filt = grepl("UA",tab4$bestSetting)
sexIA = rbind(sexIA_eGFR[!filt,],sexIA_UA[filt,])
setDT(sexIA)
sexIA[,fdr_diff := p.adjust(p_diff, method = "fdr")]

stopifnot(sexIA$SNP == tab4$indexSNP)
tab4 = cbind(tab4,sexIA[,c(3:12)])
tab4[,sexIA_FDRsig :=F]
tab4[fdr_diff<0.05,sexIA_FDRsig :=T]

coloc_sexIA = rbind(coloc_sexIA[1:7,],coloc_sexIA[7:22,])
coloc_sexIA[,region :=gsub("Region ","",locus)]
stopifnot(coloc_sexIA$region == tab4$region)
tab4 = cbind(tab4,coloc_sexIA[,c(4:9)])
coloc_sexIA = rbind(coloc_sexIA[1:7,],coloc_sexIA[7:22,])
tab4[,sexIA_coloc:="inconclusive"]
tab4[PP.H4.abf>0.75,sexIA_coloc:="both sexes - shared"]
tab4[PP.H3.abf>0.75,sexIA_coloc:="both sexes - independent"]
tab4[PP.H2.abf>0.75,sexIA_coloc:="female driven"]
tab4[PP.H1.abf>0.75,sexIA_coloc:="male driven"]
tab4

#' # Get Sup Tab 5 ####
#' ***
#' Cross-phenotype comparision (--> see script 04)
load("../results/04_lookup_TopHits_matchingSettings.RData")
tab5 = ShorterTable
tab5

#' # Get Sup Tab 6 ####
#' ***
#' Genome-wide significant & independent hit per region (--> see script 05)
#'
tab6 = fread("../results/05_b_Cojo_Select_Results.txt")
tab6[,phenotype := gsub("_all","_ALL",phenotype)]
tab6[,phenotype := gsub("_male","_MALE",phenotype)]
tab6[,phenotype := gsub("_female","_FEMALE",phenotype)]

myNames = names(tab6)[c(1:3,18,5:17)]
colsOut<-setdiff(colnames(tab6),myNames)
tab6[,get("colsOut"):=NULL]
setcolorder(tab6,myNames)
names(tab6) = c("phenotype","region","size.region","rsID","SNPID_UKBB","position","effect_allele","eaf","beta","SE","pvalue","estimated_Ne",
               "eaf_UKBB","beta_joint","SE_joint","pvalue_joint","LD_r")
filt = tab6$LD_r ==0
ld = tab6$LD_r[!filt]
x=grep("F",filt)
x2 = x+1
tab6[x,LD_r2:=ld]
tab6[x2,LD_r2:=ld]
tab6[,LD_r:=NULL]
tab6

#' # Get Sup Tabs 7 ####
#' ***
#'  Annotation of credible sets (a-e for eGFR and UA in their respective settings)
#' 
ToDoList3 = data.table(pheno = c("eGFR_ALL","eGFR_FEMALE","eGFR_MALE","UA_ALL","UA_MALE"),
                       files = c("../../../10_metaGWAS/01_eGFR_allEth_sex_combined/08_credSets/gwasresults_V5/synopsis/topliste_tabdelim/topliste_2022-11-04_credSets.txt",
                                 "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_female_V5/synopsis/topliste_tabdelim/topliste_2022-11-04_credSets_female.txt",
                                 ("../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_male_V4/synopsis/topliste_tabdelim/topliste_2022-07-22_credSets_male.txt"),
                                 ("../../../10_metaGWAS/03_uric_acid_allEth_sex_combined/08_credSets/gwasresults_V4/synopsis/topliste_tabdelim/topliste_2022-07-21_credSets.txt"),
                                 ("../../../10_metaGWAS/03_uric_acid_allEth_sex_stratified/08_credSets/gwasresults_male_V4/synopsis/topliste_tabdelim/topliste_2022-07-22_credSets_male.txt")))

tab7 = foreach(i=1:dim(ToDoList3)[1])%do%{
  #i=1
  myRow = ToDoList3[i,]
  
  tab = fread(myRow$files)
  
  stopifnot(tab$litsnp == tab$provided_name)
  
  myNames2 =  c("snp","CredSet","PostProb","SumProb","cyto","pos","tagger","r2_tagger","tagsnp","effect_allele","other_allele","eaf",
                "beta","SE","logP", "I2","nearestgenes","gene_biotype","nearestgene","Eigen","EigenPC","CADD_scaled","DANN",
                "GWAVA_Region","GWAVA_TSS","GWAVA_Unmatched","regulome_score","regulome_score_numeric","regulome_details",
                "corinfo_gwas2","cisgene","transgene","coremine_genes","hgnc4pathway","entrez4pathway","KEGG","reactome","DOSE",
                "GO","tissues")
  myNames = names(tab)[c(1:4,13,15:18,61,62,64,70:72,74,19:42)]
  stopifnot(sum(myNames == myNames2)==34)
  colsOut<-setdiff(colnames(tab),myNames)
  tab[,get("colsOut"):=NULL]
  setcolorder(tab,myNames)
  names(tab) = myNames2
  tab[,phenotype := myRow$pheno]
  tab
}
tab7 = rbindlist(tab7)
table(tab7$phenotype)
tab7

#' # Get Sup Tabs 8 ####
#' ***
#' Co-localization with eQTLs (--> see script 07)
#' 
#' Problem: no indication, which signal in the cond setting is used!!
#' 
#' --> fix this later!!
#' 
#' 
tab8_cond = data.table(read_excel("../results/07_coloc_cond_summary.xlsx"))
tab8_uncond = data.table(read_excel("../results/07_coloc_summary.xlsx"))

tab8_a = copy(tab8_uncond)
setnames(tab8_a,"locus","region")
setnames(tab8_a,"trait1","GWAMA_phenotype")
setnames(tab8_a,"trait2","tissue")
tab8_a = tab8_a[,c(1,3,2,4:10)]
tab8_a[,GWAMA_phenotype:=gsub(" [(]","_",GWAMA_phenotype)]
tab8_a[,GWAMA_phenotype:=gsub("[)]","",GWAMA_phenotype)]
tab8_a[,GWAMA_phenotype:=gsub("Uric Acid","UA",GWAMA_phenotype)]
tab8_a[,region :=gsub("Region ","",region )]

tab8_a[,dumID := paste(region,GWAMA_phenotype,sep="__")]
goodDumIDs = tab6[,paste(region,phenotype,sep="__")]
table(is.element(tab8_a$dumID,goodDumIDs))
tab8_a = tab8_a[dumID %in% goodDumIDs,]
setorder(tab8_a,region,GWAMA_phenotype,gene)
tab8_a[,dumID:=NULL]

#' # Get Sup Tabs 9 ###
#' ***
#' Replication in *HUNT* 
#' 
load("../results/09_replication_HUNT.RData")
tab9 = copy(data_Merge)
tab9

#' # Get Sup Tabs 10 ###
#' ***
#' Replication of *Graham et al*, *Kanai et al* and *Sakaue et al* results
#' 

tab10 = fread("../results/11_Look_Up_GWAS_hits_eGFR_UA_only.txt")
tab10[grepl("A cross-population atlas of genetic associations for 220",STUDY),Authors := "Sakaue et al. (2021)"]
tab10[grepl("A cross-population atlas of genetic associations for 220",STUDY),DOI := "10.1038/s41588-021-00931-x"]
tab10[grepl("Genetic analysis of quantitative traits in the Japanese",STUDY),Authors := "Kanai et al. (2018)"]
tab10[grepl("Genetic analysis of quantitative traits in the Japanese",STUDY),DOI := "10.1038/s41588-018-0047-6"]
tab10[grepl("Sex-specific and pleiotropic effects ",STUDY),Authors := "Graham et al. (2019)"]
tab10[grepl("Sex-specific and pleiotropic effects ",STUDY),DOI := "10.1038/s41467-019-09861-z"]

tab10[,`INITIAL SAMPLE SIZE` := gsub("European ancestry individuals","EUR",`INITIAL SAMPLE SIZE`)]
tab10[,`INITIAL SAMPLE SIZE` := gsub("East Asian ancestry individuals","EA",`INITIAL SAMPLE SIZE`)]
tab10[,`INITIAL SAMPLE SIZE` := gsub("European ancestry","EUR",`INITIAL SAMPLE SIZE`)]
tab10[,`INITIAL SAMPLE SIZE` := gsub("East Asian ancestry","EA",`INITIAL SAMPLE SIZE`)]
tab10[,`INITIAL SAMPLE SIZE` := gsub("Japanese ancestry individuals","JAP",`INITIAL SAMPLE SIZE`)]

names(tab10)
tab10 = tab10[,c(32,33,2,3,5:31)]
tab10[!is.na(eGFR.P),table(eGFR.P<0.05)]
tab10[!is.na(UA.P),table(UA.P<0.05)]
tab10[!is.na(UA.P),table(UA.P<0.05,eGFR.P<0.05)]
tab10


#' # Get Sup Tabs 11 ###
#' ***
#' MR-Mega
#' 
#' Not yet available
#' 

#' # Save tables ###
#' ***
tosave4 = data.table(data = c("tab0",#"tab1","tab2", 
                              "tab3", "tab4","tab5","tab6",
                              "tab7","tab8_a","tab9","tab10"), 
                     SheetNames = c("Content",#"TableS1","TableS2", 
                                    "TableS3", "TableS4","TableS5", "TableS6",
                                    "TableS7","TableS8a","TableS9","TableS10"))
excel_fn = "../tables/SupplementalTables.xlsx"
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

