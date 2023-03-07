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
#' 6) Genome-wide significant & independent hit per locus (--> see script 05)
#' 7) Annotation of credible sets (a-e for eGFR and UA in their respective settings)
#' 8) Co-localization with eQTLs (--> see script 07)
#' 9) Replication in *HUNT* (--> see script 09)
#' 10) Replication of *Graham et al*, *Kanai et al* and *Sakaue et al* results (--> see script 11)
#' 11) Summary of MR-Mega results (--> see scripts 12)**
#' 12) Annotation of additional MR-Mega results 
#' 13) Summary of Lookup of sex-biased gene expression of candidate genes (--> see script 13)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:13)),
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
                              "Summary of MR-Mega results of all index SNPs",
                              "Annotation of additional MR-Mega results",
                              "Lookup of sex-biased gene expression"),
                    Source = c("done manually",
                               "done manually",
                               "done within script ST_SupplementalTables.R",
                               "../results/02_ and ../results/03_",
                               "../results/04_lookup_TopHits_matchingSettings.RData",
                               "../results/05_b_Cojo_Select_Results.txt",
                               "../../../10_metaGWAS/PHENOTYPE/08_credSets/gwasresults_V6/synopsis/topliste_tabdelim/topliste_2022-11-24_credSets.txt",
                               "../results/07_coloc_eQTLs.RData",
                               "../results/09_replication_HUNT.RData",
                               "../results/11_Look_Up_GWAS_hits_eGFR_UA_only.txt",
                               "../results/12_SNPs_table1_in_MR_MEGA.txt",
                               "../../../12_MR-MEGA/01_annotation/gwasresults_V6/synopsis/topliste_tabdelim/topliste_2022-11-24_credSets.txt",
                               "done with script 13_Lookup-SexBiasedGeneExpression.R"))
  
  tab0
  
}

#' # Get Sup Tab 3 ####
#' ***
#' Sample Sizes & SNP Numbers, and inflation factor $\lambda$ per phenotype
#' 
#' Number of cases & controls for CKD and BUN will be added later (extract information out of autosomal publications)
#' 
#' Please note: I load the initial meta-analyses files. These files will not be shared later. We only provide statistics for SNPs valid in at least one phenotype (e.g. valid in UA but not eGFR --> listed). SNPs that are invalid in all phenotypes are not listed in the data files. However, we want to report with how many SNPs we started, and how many remained after filtering. A total of 325,770 SNPs had at least on valid association. 
#' 
{
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
                     n_SNPs_filtered = n_SNPs_filtered)
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
  
  tab3_annot = data.table(column = names(tab3),
                          description = c("Analyzed phenotype and setting",
                                          "Maximal number of participating studies",
                                          "Maximal number of individuals",
                                          "Number of SNPs before any QC was applied (number of SNPs analyzed in raw meta-analysis)",
                                          "Inflation factor lambda before any QC was applied",
                                          "Number of SNPs after QC (filtering for number of studies >=10, imputation info score >=0.8, minor allele frequency >=0.02, heterogeneity I^2 <=0.8)"))
  
  
}


#' # Get Sup Tab 4 ####
#' ***
#' Comparisons between the sexes (interaction and co-localization --> see script 02 and 03)
#' 
{
  load("../results/04_lookup_TopHits_allOtherTraits.RData")
  sexIA_eGFR = fread("../results/02_sex_ia_egfr.txt")
  sexIA_UA = fread("../results/02_sex_ia_uric_acid.txt")
  coloc_sexIA = fread("../results/03_coloc_sex_ia.txt")
  
  tab4 = copy(WideTable)
  tab4 = tab4[,c(1:4)]
  setnames(tab4,"Cytoband","cytoband")
  setnames(tab4,"rsID","indexSNP")
  setnames(tab4,"topPheno","bestSetting")
  setnames(tab4,"region","locus_NR")
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
  stopifnot(coloc_sexIA$region == tab4$locus_NR)
  tab4 = cbind(tab4,coloc_sexIA[,c(4:9)])
  coloc_sexIA = rbind(coloc_sexIA[1:7,],coloc_sexIA[7:22,])
  tab4[,sexIA_coloc:="inconclusive"]
  tab4[PP.H4.abf>0.75,sexIA_coloc:="both sexes - shared"]
  tab4[PP.H3.abf>0.75,sexIA_coloc:="both sexes - independent"]
  tab4[PP.H2.abf>0.75,sexIA_coloc:="female driven"]
  tab4[PP.H1.abf>0.75,sexIA_coloc:="male driven"]
  tab4
  
  tab4_annot = data.table(column = names(tab4),
                          description = c("Number of associated loci (1-15: eGFR, 16-22: UA)",
                                          "Genomic cytoband of index SNP",
                                          "SNP with lowest p-value in this locus",
                                          "Best phenotype and setting",
                                          "Beta estimate in males (using respective best phenotypes)",
                                          "Standard error in males (using respective best phenotypes)",
                                          "P-value in males (using respective best phenotypes)",
                                          "Beta estimate in females (using respective best phenotypes)",
                                          "Standard error in females (using respective best phenotypes)",
                                          "P-value in females (using respective best phenotypes)",
                                          "Difference of effect estimates (males - females)",
                                          "Standard error of difference (correcting for beta correlation using Sprearmans rho (eGFR: 0.1604, UA: 0.1285))",
                                          "P-value of difference",
                                          "FDR corrected p-value of difference",
                                          "TRUE/FALSE flag indicating significant sex-interaction after FDR correction",
                                          "Number of SNPs included in co-localization analysis per loci",
                                          "Posterior probability for hypothesis 0: neither trait associated",
                                          "Posterior probability for hypothesis 1: only trait 1 associated (males)",
                                          "Posterior probability for hypothesis 2: only trait 2 associated (females)",
                                          "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                          "Posterior probability for hypothesis 4: both trait associated, shared signal",
                                          "Summary of co-localization result, using 0.75 as threshold"))
  
  
}

#' # Get Sup Tab 5 ####
#' ***
#' Cross-phenotype comparision (--> see script 04)
{
  load("../results/04_lookup_TopHits_matchingSettings.RData")
  tab5_a = ShorterTable
  load("../results/04_lookup_TopHits_allSettings_logP.RData")
  tab5_b = ShorterTable_v2
  tab5_c = fread("../results/08_b_coloc_overlap.txt")
  
  names(tab5_a)[1:4] = names(tab4)[1:4]
  names(tab5_b)[1:4] = names(tab4)[1:4]
  dummy = names(tab5_a)[c(1:4,15:22)]
  tab5_a_annot = data.table(column = dummy,
                          description = c("Number of associated locus (1-15: eGFR, 16-22: UA)",
                                          "Genomic cytoband of index SNP",
                                          "SNP with lowest p-value in this locus",
                                          "Best phenotype and setting",
                                          "Sample size in phenotype indicated in column name (using respective best setting)",
                                          "Effect allele frequency in phenotype indicated in column name (using respective best setting)",
                                          "Beta estimate in phenotype indicated in column name (using respective best setting)",
                                          "Standard error in phenotype indicated in column name (using respective best setting)",
                                          "P-value in phenotype indicated in column name (using respective best setting)",
                                          "One-sided p-value in phenotype indicated in column name (using respective best setting, not calculated for eGFR and UA)",
                                          "TRUE/FALSE flag indicating significant nominal association in phenotype indicated in column name with discordant effect direction compared to eGFR (using respective best setting)",
                                          "TRUE/FALSE flag indicating significant nominal association in phenotype indicated in column name with concordant effect direction compared to UA (using respective best setting)"))
  
  tab5_b_annot = data.table(column = names(tab5_b)[1:5],
                            description = c("Number of associated locus (1-15: eGFR, 16-22: UA)",
                                            "Genomic cytoband of index SNP",
                                            "SNP with lowest p-value in this locus",
                                            "Best phenotype and setting",
                                            "-log10 transformed p-value in phenotype indicated in column name and setting ALL (one-sided p-value for BUN, CKD, UACR, MA, and gout; two-sided for eGFR and UA)"))
  
  tab5_c[,region1 := gsub("Region ","",region1)]
  tab5_c[,region2 := gsub("Region ","",region2)]
  tab5_c[,region1 := as.numeric(region1)]
  tab5_c[,region2 := as.numeric(region2)]
  setnames(tab5_c,"region1","locus1")
  setnames(tab5_c,"region2","locus2")
  
  matched = match(tab5_c$locus1,tab4$locus_NR)
  tab5_c[,cytoband := tab4[matched,cytoband]]
  tab5_c = tab5_c[,c(1,2,11,3:10)]
  
  tab5_c_annot = data.table(column = names(tab5_c),
                          description = c("Number of associated loci for eGFR",
                                          "Number of associated loci for UA",
                                          "Genomic cytoband of gene",
                                          "Analyzed phenotype and setting of eGFR",
                                          "Analyzed phenotype and setting of UA",
                                          "Number of SNPs included in co-localization analysis per locus",
                                          "Posterior probability for hypothesis 0: neither trait associated",
                                          "Posterior probability for hypothesis 1: only trait 1 associated (eGFR)",
                                          "Posterior probability for hypothesis 2: only trait 2 associated (UA)",
                                          "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                          "Posterior probability for hypothesis 4: both trait associated, shared signal"))
  
}

#' # Get Sup Tab 6 ####
#' ***
#' Genome-wide significant & independent hit per locus (--> see script 05)
#'
{
  tab6 = fread("../results/05_b_Cojo_Select_Results.txt")
  tab6_help = fread("../results/06_SNP_numbers_Credible_Sets.txt")
  
  tab6[,dumID := paste(phenotype,region,rsID,sep="::")]
  tab6_help[,dumID := paste(phenotype,region,SNP,sep="::")]
  stopifnot(sum(is.element(tab6$dumID,tab6_help$dumID))==45)
  matched = match(tab6$dumID,tab6_help$dumID)
  tab6[,CredSetSize := tab6_help[matched,CredSet.99]]
  tab6[,dumID := NULL]
  
  tab6[,phenotype := gsub("_all","_ALL",phenotype)]
  tab6[,phenotype := gsub("_male","_MALE",phenotype)]
  tab6[,phenotype := gsub("_female","_FEMALE",phenotype)]
  
  myNames = names(tab6)[c(2:3,19,18,5,1,6:17)]
  colsOut<-setdiff(colnames(tab6),myNames)
  tab6[,get("colsOut"):=NULL]
  setcolorder(tab6,myNames)
  names(tab6) = c("locus_NR","size.region","CredSetSize","rsID","SNPID_UKBB","phenotype","position","effect_allele",
                  "eaf","beta","SE","pvalue","estimated_Ne",
                  "eaf_UKBB","beta_joint","SE_joint","pvalue_joint","LD_r")
  filt = tab6$LD_r ==0
  ld = tab6$LD_r[!filt]
  x=grep("F",filt)
  x2 = x+1
  tab6[x,LD_r2:=ld]
  tab6[x2,LD_r2:=ld]
  tab6[,LD_r:=NULL]
  tab6
  
  
  tab6_annot = data.table(column = names(tab6),
                          description = c("Number of associated locus (1-15: eGFR, 16-22: UA)",
                                          "Range of tested locus (not necessary centered around the index SNP)",
                                          "Number of SNPs in 99% credible set",
                                          "Independent SNP",
                                          "SNP ID as in UKBB data (reference data in GCTA COJO analysis)",
                                          "Analyzed phenotype and setting",
                                          "Base position of independent SNP",
                                          "Effect allele in CKDGen",
                                          "Effect allele frequency in CKDGen",
                                          "Beta estimate in CKDGen",
                                          "Standard error in CKDGen",
                                          "P-value in CKDGen",
                                          "Estimated effective sample size",
                                          "Effect allele frequency in UKBB",
                                          "Beta estimate in joint analysis",
                                          "Standard error in joint analysis",
                                          "P-value in joint analysis",
                                          "LD r^2 between the two independent signals (according to UKBB)"))
  
  
  
}

#' # Get Sup Tabs 7 ####
#' ***
#'  Annotation of credible sets (a-e for eGFR and UA in their respective settings):
#' 
#' Ensembl
#' CADD score
#' Regulome score
#' GWAS catalog
#' GTEx v8 eQTLs (cis)
#' Pathways (KEGG, GO, DOSE and Reactome)
#' 
{
  ToDoList3 = data.table(pheno = c("eGFR_ALL","eGFR_FEMALE","eGFR_MALE","UA_ALL","UA_MALE"),
                         files = c(paste0(path_CS_eGFR_ALL,"topliste_2022-11-29_credSets.txt"),
                                      paste0(path_CS_eGFR_FEMALE,"topliste_2022-11-29_credSets_female.txt"),
                                      paste0(path_CS_eGFR_MALE,"topliste_2022-11-29_credSets_male.txt"),
                                      paste0(path_CS_UA_ALL,"topliste_2022-11-29_credSets.txt"),
                                      paste0(path_CS_UA_MALE,"topliste_2022-11-29_credSets_male.txt")))
  
  tab7 = foreach(i=1:dim(ToDoList3)[1])%do%{
    #i=1
    myRow = ToDoList3[i,]
    
    tab = fread(myRow$files)
    
    stopifnot(tab$litsnp == tab$provided_name)
    tab = distinct(tab)
    
    myNames2 =  c("snp","CredSet","PostProb","SumProb","cyto","pos","tagger","r2_tagger","tagsnp","effect_allele","other_allele",
                  "eaf","nSamples","info","beta","SE","logP", "I2",
                  "nearestgenes","gene_biotype","nearestgene",
                  "CADD_scaled","regulome_score","regulome_score_numeric","regulome_details",
                  "corinfo_gwas2","cisgene","coremine_genes",
                  "KEGG","reactome","DOSE","GO")
    myNames = names(tab)[c(1:4,13,15:18,61,62,64:66,70:72,74,19:21,24,29:33,35,38:41)]
    stopifnot(sum(myNames == myNames2)==24)
    colsOut<-setdiff(colnames(tab),myNames)
    tab[,get("colsOut"):=NULL]
    setcolorder(tab,myNames)
    names(tab) = myNames2
    tab[,phenotype := myRow$pheno]
    if(myRow$pheno == "eGFR_MALE"){
      tab[cyto == "Xq22.1" & CredSet == "Region1", CredSet := "Region7"]
    }
    if(myRow$pheno == "UA_ALL"){
      tab[CredSet == "Region21_SNP1", CredSet := "Region21_rs202138804"]
      tab[CredSet == "Region21_SNP2", CredSet := "Region21_rs7056552"]
      tab[CredSet == "Region22_SNP1", CredSet := "Region22_rs111884516"]
      tab[CredSet == "Region22_SNP2", CredSet := "Region22_rs4328011"]
    }
    tab
  }
  tab7 = rbindlist(tab7)
  table(tab7$phenotype)
  tab7
  tab7 = tab7[,c(33,1:32)]
  
  # update  "beta","SE","logP" with conditional statistics!
  tab7_0 = copy(tab7)
  tab7_0 = tab7_0[!grepl("Region.._",CredSet) & !grepl("Region._",CredSet)]
  
  tab7_2 = copy(tab7)
  tab7_2 = tab7_2[grepl("Region.._",CredSet) | grepl("Region._",CredSet)]
  table(tab7_2$CredSet)
  tab7_2[,table(duplicated(snp)),by=CredSet]
  
  cond = unique(tab7_2$CredSet)
  cond
  cond_stats = list.files(path = "../temp/05_c_Cojo_cond_results/",pattern = ".cma.cojo")
  cond_stats
  myCond = data.table(cond = cond[c(2,1,6,5,3,4)],
                      filename = cond_stats)
  myCond
  
  dumTab = foreach(i=1:dim(myCond)[1])%do%{
    # i=1
    myRow = myCond[i,]
    myRow
    tab = fread(paste0("../temp/05_c_Cojo_cond_results/",myRow$filename))
    table(duplicated(tab$bp))
    stopifnot(sum(duplicated(tab$bp))==0)
    
    tab7_3 = copy(tab7_2)
    tab7_3 = tab7_3[CredSet == myRow$cond,]
    stopifnot(sum(is.element(tab7_3$pos,tab$bp))==dim(tab7_3)[1])
    
    matched = match(tab7_3$pos,tab$bp)
    table(tab7_3$pos == tab$bp[matched])
    tab7_3[,beta := tab[matched,bC]]
    tab7_3[,SE := tab[matched,bC_se]]
    tab7_3[,logP := -log10(tab[matched,pC])]
    tab7_3
  }
  dumTab = rbindlist(dumTab)
  
  tab7 = rbind(tab7_0,dumTab)
  setorder(tab7,-PostProb)
  tab7[,CredSet := gsub("Region","Locus_NR_",CredSet)]
  
  tab7_annot = data.table(column = names(tab7),
                          description = c("Analyzed phenotype and setting",
                                          "SNP ID of the marker.  Mostly from dbSNP, sometimes constructed like chr1:42147691:D (chromosome, position, alleletype)",
                                          "ID of the credible set (locus number + independent signal incase of multiple signals at a locus)",
                                          "Posterior Probability of the SNP depending of the given phenotype",
                                          "Cummulative sum of Probabilities per credible set",
                                          "Cytoband with 850 resolution",
                                          "Position of SNP in basepairs on the chromosome (hg19)",
                                          "SNP that tags the marker in column snp",
                                          "LD r-square value between SNP and tagSNP",
                                          "TRUE/FALSE flag indicating marker in SNP is tagger",
                                          "Allele for which effect sizes were calculated",
                                          "Other allele",
                                          "Effect allele frquency",
                                          "Weighted Imputation info score",
                                          "Sample size per SNP",
                                          "Beta estimate (conditional in case of locus 9, 21, and 22 in setting ALL)",
                                          "Standard error (conditional in case of locus 9, 21, and 22 in setting ALL)",
                                          "-log10 transformed p-value (conditional in case of locus 9, 21, and 22 in setting ALL)",
                                          "Heterogenetity I-squared",
                                          "HGNC (Ensembl) symbols of the nearest genes as specified in the settings with information on distance to the marker in SNP including functional relevance if within or within the flanking 5 kb of a gene",
                                          "Functional relevance of the gene and its validation level  for proximate genes(Ensembl)",
                                          "HGNC (Ensembl) full name of nearest gene that actually has got a full name",
                                          "PHRED-like (-10*log10(rank/total)) scaled C-score ranking a variant relative to all possible substitutions of the human genome (8.6x10^9). For details see tab 'deleteriousness'.  Like explained above, a scaled C-score of greater or# equal 10 indicates that these are predicted to be the 10% most deleterious substitutions that you can do to the human genome, a score of greater or equal 20 indicates the 1% most deleterious and so on. If you would like to apply a cutoff on deleteriousness, eg to identify potentially pathogenic variants, we would suggest to put a cutoff somewhere between 10 and 20. Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants. However, there is not a natural choice here -- it is always arbitrary. We therefore recommend integrating C-scores with other evidence and to rank your candidates for follow up rather than hard filtering.",
                                          "Deleteriousness / functional variant relevance score",
                                          "Deleteriousness / functional variant relevance score",
                                          "Deleteriousness / functional variant relevance score",
                                          "r-square value between marker in SNP and a hit from GWAS catalog for named phenotype",
                                          "cis-eQTL genes for which reported eQTL-SNP has specified min R2 with marker in SNP",
                                          "List of all genes reported for marker in SNP in nearest genes, and cis-eQTL genes",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in KEGG pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in reactome pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in disease ontology pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in GO categories pathways"))
  
}

#' # Get Sup Tabs 8 ####
#' ***
#' Co-localization with eQTLs (--> see script 07)
#' 
#' 
{
  load("../results/07_b_coloc_eQTLs.RData")
  tab8 = copy(ColocTable)
  
  names(tab8)
  setnames(tab8,"trait1","GWAMA_phenotype")
  setnames(tab8,"trait2","tissue")
  setnames(tab8,"region","locus_NR")
  
  tab8_annot = data.table(column = names(tab8),
                          description = c("Number of associated loci (1-15: eGFR, 16-22: UA)",
                                          "Genomic cytoband of gene",
                                          "Analyzed phenotype and setting",
                                          "Analyzed gene expression of GTEx or NEPTUNE",
                                          "Analyzed tissue for gene expression",
                                          "Number of SNPs included in co-localization analysis per locus",
                                          "Posterior probability for hypothesis 0: neither trait associated",
                                          "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen)",
                                          "Posterior probability for hypothesis 2: only trait 2 associated (gene expression)",
                                          "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                          "Posterior probability for hypothesis 4: both trait associated, shared signal"))
  
  
}

#' # Get Sup Tabs 9 ###
#' ***
#' Replication in *HUNT* 
#' 
{
  load("../results/09_replication_HUNT_summary.RData")
  tab9 = copy(result)
  tab9 = tab9[!grepl("rs111410539:",SNP),]
  setnames(tab9,"region","locus_NR")
  tab9_annot = data.description
  tab9_annot[1,column := "locus_NR"]
}


#' # Get Sup Tabs 10 ###
#' ***
#' Replication of *Graham et al*, *Kanai et al* and *Sakaue et al* results
#' 
{
  tab10 = fread("../results/11_Look_Up_GWAS_hits_eGFR_BUN_UA_only.txt")
  tab10_help = fread("../results/11_Look_Up_GWAS_hits_results.txt")
  
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
  
  table(tab10$`DISEASE/TRAIT`)
  names(tab10)
  
  tab10[,phenotype := "eGFR"]
  tab10[`DISEASE/TRAIT` == "Blood urea nitrogen levels",phenotype := "BUN"]
  tab10[`DISEASE/TRAIT` == "Serum uric acid levels",phenotype := "UA"]
  
  tab10[,CKDGen_N := eGFR.N]
  tab10[,CKDGen_EA := eGFR.effect_allele]
  tab10[,CKDGen_EAF := eGFR.EAF]
  tab10[,CKDGen_infoScore := eGFR.infoScore]
  tab10[,CKDGen_I2 := eGFR.I2]
  tab10[,CKDGen_beta := eGFR.beta]
  tab10[,CKDGen_SE := eGFR.SE]
  tab10[,CKDGen_P := eGFR.P]
  tab10[,CKDGen_logP := eGFR.logP]
  tab10[,CKDGen_invalid := eGFR.invalid_assoc]
  tab10[,CKDGen_reason := eGFR.reason_for_exclusion]
  
  tab10[phenotype == "UA",CKDGen_N := UA.N]
  tab10[phenotype == "UA",CKDGen_EA := UA.effect_allele]
  tab10[phenotype == "UA",CKDGen_EAF := UA.EAF]
  tab10[phenotype == "UA",CKDGen_infoScore := UA.infoScore]
  tab10[phenotype == "UA",CKDGen_I2 := UA.I2]
  tab10[phenotype == "UA",CKDGen_beta := UA.beta]
  tab10[phenotype == "UA",CKDGen_SE := UA.SE]
  tab10[phenotype == "UA",CKDGen_P := UA.P]
  tab10[phenotype == "UA",CKDGen_logP := UA.logP]
  tab10[phenotype == "UA",CKDGen_invalid := UA.invalid_assoc]
  tab10[phenotype == "UA",CKDGen_reason := UA.reason_for_exclusion]
  
  tab10[phenotype == "BUN",CKDGen_N := BUN.N]
  tab10[phenotype == "BUN",CKDGen_EA := BUN.effect_allele]
  tab10[phenotype == "BUN",CKDGen_EAF := BUN.EAF]
  tab10[phenotype == "BUN",CKDGen_infoScore := BUN.infoScore]
  tab10[phenotype == "BUN",CKDGen_I2 := BUN.I2]
  tab10[phenotype == "BUN",CKDGen_beta := BUN.beta]
  tab10[phenotype == "BUN",CKDGen_SE := BUN.SE]
  tab10[phenotype == "BUN",CKDGen_P := BUN.P]
  tab10[phenotype == "BUN",CKDGen_logP := BUN.logP]
  tab10[phenotype == "BUN",CKDGen_invalid := BUN.invalid_assoc]
  tab10[phenotype == "BUN",CKDGen_reason := BUN.reason_for_exclusion]
  
  tab10[,cytoband:= tab10_help$REGION]
  tab10[,pos_hg19:= tab10_help$eGFR.position]
  tab10[,litCanGene := tab10_help$`REPORTED GENE(S)`]
  matched = match(tab10$cytoband,tab4$cytoband)
  tab10[,locus_NR := tab4[matched,locus_NR]]
  
  names(tab10)[c(46,47,2,3,5:10,60:63,45,48:59)]
  tab10 = tab10[,c(46,47,2,3,5:10,60:63,45,48:59)]
  tab10 = tab10[!is.na(CKDGen_invalid),]
  
  tab10_annot = data.table(column = names(tab10),
                           description = c("First Author et al. of looked-up study",
                                           "Digital Object Identifier of looked-up study",
                                           "Analyzed phenotype in looked-up study",
                                           "Sample Size of looked-up study",
                                           "Base position in looked-up study (hg38)",
                                           "SNP ID and effect allele in looked-up study",
                                           "SNP ID used for matching",
                                           "Effect allele frequency in looked-up study",
                                           "P-value in looked-up study",
                                           "Odds ratio, Z-score or beta estimate in looked-up study",
                                           "Genomic cytoband of literature SNP",
                                           "Base position in our study (hg19)",
                                           "Candidate Gene as reported in GWAS Catalog",
                                           "Corresponding locus number in our analyses",
                                           "TRUE/FALSE flag indicating successful replication in CKDGen",
                                           "CKDGen phenotype used for replication",
                                           "Sample size in CKDGen",
                                           "Effect allele in CKDGen",
                                           "Effect allele frquency in CKDGen",
                                           "Weighted imputation info score in CKDGen",
                                           "Heterogeneity I^2 in CKDGen",
                                           "Beta estimate in CKDGen",
                                           "Standard error in CKDGen",
                                           "P-value in CKDGen",
                                           "-log10 transformed p-value in CKDGen",
                                           "TRUE/FALSE flag indicating validity of SNP for CKDGen (F=valid)",
                                           "Reason to exclude SNP for CKDGen"))
  
}

#' # Get Sup Tabs 11 ###
#' ***
#' MR-Mega of our index SNPx
#' 
{
  tab11 = fread("../results/12_SNPs_table1_in_MR_MEGA.txt")
  chi2_p_anc_het = pchisq(tab11$chisq_ancestry_het,df=tab11$ndf_ancestry_het,lower.tail = F)
  chi2_p_res_het = pchisq(tab11$chisq_residual_het,df=tab11$ndf_residual_het,lower.tail = F)
  chi2_p_assoc = pchisq(tab11$chisq_association,df=tab11$ndf_association,lower.tail = F)
  
  plot(tab11$`P-value_association`,chi2_p_assoc)
  abline(0,1)
  plot(tab11$`P-value_ancestry_het`,chi2_p_anc_het)
  abline(0,1)
  plot(tab11$`P-value_residual_het`,chi2_p_res_het)
  abline(0,1)
  
  tab11[,`P-value_association` := chi2_p_assoc]
  tab11[,`P-value_ancestry_het` := chi2_p_anc_het]
  tab11[,`P-value_residual_het` := chi2_p_res_het]
  
  stopifnot(tab11$rsID == tab4$indexSNP)
  tab11 = cbind(tab4$locus_NR,tab11)
  names(tab11)[1] = "locus_NR"
  tab11_annot = data.table(column = names(tab11),
                           description = c("Number of associated locus (1-15: eGFR, 16-22: UA)",
                                           "SNP ID in CKDGen",
                                           "SNP ID used for matching",
                                           "Chromosome",
                                           "Base position (hg19)",
                                           "Effect allele",
                                           "Non-effect allele",
                                           "Effect allele frequency",
                                           "Sample size in CKDGen",
                                           "Number of studies in CKDGen",
                                           "Effect direction across cohorts (+ if the effect allele effect was positive, - if negative, 0 if the effect was zero, ? if marker was not available in cohort)",
                                           "Effect of first PC of meta-regression",
                                           "Standard error of the effect of first PC of meta-regression",
                                           "Effect of second PC of meta-regression",
                                           "Standard error of the effect of second PC of meta-regression",
                                           "Effect of third PC of meta-regression",
                                           "Standard error of the effect of third PC of meta-regression",
                                           "Effect of fourth PC of meta-regression",
                                           "Standard error of the effect of fourth PC of meta-regression",
                                           "Chisq value of the association",
                                           "Number of degrees of freedom of the association",
                                           "P-value of the association",
                                           "Chisq value of the heterogeneity due to different ancestry",
                                           "Number of degrees of freedom of the heterogeneity due to different ancestry",
                                           "P-value of the heterogeneity due to different ancestry",
                                           "Chisq value of the residual heterogeneity",
                                           "Number of degrees of freedom  of the residual heterogeneity",
                                           "Pvalue of the residual heterogeneity",
                                           "Log of Bayes factor",
                                           "Reason why marker was not analysed in MR-MEGA",
                                           "TRUE/FALSE flag indicating validity of SNP for CKDGen (F=valid)",
                                           "Phenotype used in CKDGen"))
  
}

#' # Get Sup Tabs 12 ###
#' ***
#' MR-Mega of additional results
#' 
{
  tab12 = fread(paste0(path_UA_MRMEGA,"topliste_2022-11-29_MR_Mega.txt"))
  tab12_help = fread("../results/12_MR_Mega_Hits_including_SNP_from_table1.txt")
  tab12_help = tab12_help[-1,]
  
  myNames2 =  c("snp","cyto","pos","effect_allele","other_allele",
                "eaf","info","nSamples","beta","SE","logP", "I2",
                "nearestgenes","gene_biotype","nearestgene",
                "CADD_scaled","regulome_score","regulome_score_numeric","regulome_details",
                "corinfo_gwas2","cisgene","coremine_genes",
                "KEGG","reactome","DOSE","GO")
  myNames = names(tab12)[c(1,10,12,58,59,61,63,62,67:69,71,16:18,21,26:30,32,35:38)]
  stopifnot(sum(myNames == myNames2)==18)
  colsOut<-setdiff(colnames(tab12),myNames)
  tab12[,get("colsOut"):=NULL]
  setcolorder(tab12,myNames)
  names(tab12) = myNames2
  
  stopifnot(tab12$pos == tab12_help$Position)
  stopifnot(tab12$effect_allele == tab12_help$EA)
  stopifnot(tab12$nSamples == tab12_help$Nsample)
  
  tab12 = cbind(tab12[,1:12],tab12_help[,8: 27],tab12[,13:26])
  
  tab12_annot = data.table(column = names(tab12),
                          description = c("SNP ID of the marker",
                                          "Cytoband with 850 resolution",
                                          "Position of SNP in basepairs on the chromosome (hg19)",
                                          "Allele for which effect sizes were calculated",
                                          "Other allele",
                                          "Effect allele frequency",
                                          "Weighted imputation info score",
                                          "Sample size per SNP",
                                          "Beta estimate in meta-analysis",
                                          "Standard error in meta-analysis",
                                          "-log10 transformed p-value in meta-analysis",
                                          "Heterogenetity I-squared in meta-analysis",
                                          "Number of studies in CKDGen",
                                          "Effect direction across cohorts (+ if the effect allele effect was positive, - if negative, 0 if the effect was zero, ? if marker was not available in cohort)",
                                          "Effect of first PC of meta-regression in MR-MEGA",
                                          "Standard error of the effect of first PC of meta-regression in MR-MEGA",
                                          "Effect of second PC of meta-regression in MR-MEGA",
                                          "Standard error of the effect of second PC of meta-regression in MR-MEGA",
                                          "Effect of third PC of meta-regression in MR-MEGA",
                                          "Standard error of the effect of third PC of meta-regression in MR-MEGA",
                                          "Effect of fourth PC of meta-regression in MR-MEGA",
                                          "Standard error of the effect of fourth PC of meta-regression in MR-MEGA",
                                          "Chisq value of the association in MR-MEGA",
                                          "Number of degrees of freedom of the association in MR-MEGA",
                                          "P-value of the association in MR-MEGA",
                                          "Chisq value of the heterogeneity due to different ancestry in MR-MEGA",
                                          "Number of degrees of freedom of the heterogeneity due to different ancestry in MR-MEGA",
                                          "P-value of the heterogeneity due to different ancestry in MR-MEGA",
                                          "Chisq value of the residual heterogeneity in MR-MEGA",
                                          "Number of degrees of freedom  of the residual heterogeneity in MR-MEGA",
                                          "Pvalue of the residual heterogeneity in MR-MEGA",
                                          "Log of Bayes factor in MR-MEGA",
                                          "HGNC (Ensembl) symbols of the nearest genes as specified in the settings with information on distance to the marker in SNP including functional relevance if within or within the flanking 5 kb of a gene",
                                          "Functional relevance of the gene and its validation level  for proximate genes(Ensembl)",
                                          "HGNC (Ensembl) full name of nearest gene that actually has got a full name",
                                          "PHRED-like (-10*log10(rank/total)) scaled C-score ranking a variant relative to all possible substitutions of the human genome (8.6x10^9). For details see tab 'deleteriousness'.  Like explained above, a scaled C-score of greater or# equal 10 indicates that these are predicted to be the 10% most deleterious substitutions that you can do to the human genome, a score of greater or equal 20 indicates the 1% most deleterious and so on. If you would like to apply a cutoff on deleteriousness, eg to identify potentially pathogenic variants, we would suggest to put a cutoff somewhere between 10 and 20. Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants. However, there is not a natural choice here -- it is always arbitrary. We therefore recommend integrating C-scores with other evidence and to rank your candidates for follow up rather than hard filtering.",
                                          "Deleteriousness / functional variant relevance score",
                                          "Deleteriousness / functional variant relevance score",
                                          "Deleteriousness / functional variant relevance score",
                                          "r-square value between marker in SNP and a hit from GWAS catalog for named phenotype",
                                          "cis-eQTL genes for which reported eQTL-SNP has specified min R2 with marker in SNP",
                                          "List of all genes reported for marker in SNP in nearest genes, and cis-eQTL genes",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in KEGG pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in reactome pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in disease ontology pathways",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in GO categories pathways"))
}

#' # Get Sup Tabs 13 ###
#' ***
#' Results of sex biased genes
#' 
{
  tab13 = data.table(candidate_genes = c("AR", "DCAF12L1", "DRP2", "EDA2R", "FAM9B", "HPRT1", "MST4"), 
                     number_tissues_higher_female = c(26, 1, 20, 0, 0, 11, 0), 
                     number_tissues_higher_male = c(2, 3, 0, 8, 0, 1, 0), 
                     kidney_cortex_sig = c("no","yes","yes","no","no","no","no"), 
                     gender_higher_expression_in_kidney_cortex = c(NA,"male","female",NA,NA,NA,NA))
  
  
  tab13_annot = data.table(column = names(tab13),
                           description = c("Name of candidate gene", 
                           "Number of tissues with a sex differential gene expression where expression is higher in females",
                           "Number of tissues with a sex differential gene expression where expression is higher in males",
                           "Is there a significant sex differential gene expression in the tissue kidney cortex?",
                           "Is the sex differential gene expression in kidney cortex higher in females or males?"))
}

  
#' # Save tables ###
#' ***
tosave4 = data.table(data = c("tab0",#"tab1","tab2", 
                              "tab3", "tab4","tab5_a","tab5_b","tab5_c","tab6","tab7","tab8",
                              "tab9","tab10","tab11","tab12", "tab13"), 
                     SheetNames = c("Content",#"TableS1","TableS2", 
                                    "TableS3", "TableS4","TableS5_a","TableS5_b","TableS5_c","TableS6","TableS7","TableS8",
                                    "TableS9","TableS10","TableS11","TableS12", "TableS13"))
excel_fn = "../tables/SupplementalTables.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

tosave4 = data.table(data = c("tab3_annot", "tab4_annot","tab5_a_annot","tab5_b_annot","tab5_c_annot","tab6_annot","tab7_annot","tab8_annot",
                              "tab9_annot","tab10_annot","tab11_annot","tab12_annot","tab13_annot"), 
                     SheetNames = c("TableS3_annot", "TableS4_annot","TableS5_a_annot","TableS5_b_annot","TableS5_c_annot","TableS6_annot","TableS7_annot",
                                    "TableS8_annot", "TableS9_annot","TableS10_annot","TableS11_annot","TableS12_annot", "TableS13_annot"))
excel_fn = "../tables/SupplementalTables_Annotation.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0,tab3,tab4,tab5_a,tab5_b,tab5_c,tab6,tab7,tab8,tab9,tab10,tab11,tab12,tab13,file = "../tables/SupplementalTables.RData")
save(tab3_annot,tab4_annot,tab5_a_annot,tab5_b_annot,tab5_c_annot,tab6_annot,tab7_annot,tab8_annot,tab9_annot,tab10_annot,tab11_annot,tab12_annot,tab13_annot,
     file = "../tables/SupplementalTables_annot.RData")

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

