#' ---
#' title: "Co-localization Part 3: Check output"
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
#' *Co-localization analysis of gene-expression quantitative trait loci*: 
#' 
#' We tested for overlapping causal variants between kidney trait associations and gene-expression quantitative trait loci (eQTLs).
#' 
#' Used eQTL databases:
#' 
#' * GTEx Analysis V8, 49 tissues, hg38-> liftover using the GTEx SNP Annotation file
#' * NEPTUNE, 2 tissues, hg19
#' 
#' We analyzed all tissues, but laid main focus in the following tissues: 
#' 
#' * GTEx adrenal gland
#' * GTEx kidney cortex
#' * GTEx muscle skeletal
#' * GTEx whole blood
#' * NEPTUNE glomerulus
#' * NEPTUNE tubulointerstitial  
#' 
#' Gene selection: 
#' 
#' * nearby genes (+/- 250 kb)
#' * cis eQTL genes (LD r^2>=0.3)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Load data ####
#' ***
load("../results/07_coloc_eQTLs.RData")
coloc = copy(ColocTable)

tab6 = fread("../results/05_b_Cojo_Select_Results.txt")
tab6[,phenotype := gsub("_all","_ALL",phenotype)]
tab6[,phenotype := gsub("_male","_MALE",phenotype)]
tab6[,phenotype := gsub("_female","_FEMALE",phenotype)]

tab6[,leadSNP := rsID]
tab6[grepl("rs111884516:",rsID),leadSNP := "rs4328011:152898261:G:A"]
tab6[grepl("rs7056552:",rsID),leadSNP := "rs202138804:133799101:AGT:A"]
tab6[grepl("rs111410539:",rsID),leadSNP := "rs181497961:106168067:G:A"]

ToDoList3 = data.table(pheno = c("eGFR_ALL","eGFR_FEMALE","eGFR_MALE","UA_ALL","UA_MALE"),
                       genelist = c("../../../10_metaGWAS/01_eGFR_allEth_sex_combined/08_credSets/gwasresults_V6/synopsis/topliste_tabdelim/proximate_genes_2022-11-24_credSets.txt",
                                    "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_female_V6/synopsis/topliste_tabdelim/proximate_genes_2022-11-24_credSets_female.txt",
                                    "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_male_V6/synopsis/topliste_tabdelim/proximate_genes_2022-11-24_credSets_male.txt",
                                    "../../../10_metaGWAS/03_uric_acid_allEth_sex_combined/08_credSets/gwasresults_V6/synopsis/topliste_tabdelim/proximate_genes_2022-11-24_credSets.txt",
                                    "../../../10_metaGWAS/03_uric_acid_allEth_sex_stratified/08_credSets/gwasresults_male_V6/synopsis/topliste_tabdelim/proximate_genes_2022-11-24_credSets_male.txt"),
                       eQTLlist = c("../../../10_metaGWAS/01_eGFR_allEth_sex_combined/08_credSets/gwasresults_V6/synopsis/topliste_tabdelim/eqtlinfo_2022-11-24_credSets.txt",
                                    "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_female_V6/synopsis/topliste_tabdelim/eqtlinfo_2022-11-24_credSets_female.txt",
                                    "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_male_V6/synopsis/topliste_tabdelim/eqtlinfo_2022-11-24_credSets_male.txt",
                                    "../../../10_metaGWAS/03_uric_acid_allEth_sex_combined/08_credSets/gwasresults_V6/synopsis/topliste_tabdelim/eqtlinfo_2022-11-24_credSets.txt",
                                    "../../../10_metaGWAS/03_uric_acid_allEth_sex_stratified/08_credSets/gwasresults_male_V6/synopsis/topliste_tabdelim/eqtlinfo_2022-11-24_credSets_male.txt"))


tab7 = foreach(i=1:dim(ToDoList3)[1])%do%{
  #i=1
  myRow = ToDoList3[i,]
  
  genes = fread(myRow$genelist)
  eQTLs = fread(myRow$eQTLlist)
  
  check = copy(tab6)
  check = check[phenotype == myRow$pheno,]
  
  genes = genes[markername %in% check$rsID,]
  genes[,dumID := paste(genename,markername,sep="__")]
  genes = genes[!duplicated(dumID),]
  
  eQTLs = eQTLs[snps %in% check$rsID,]
  eQTLs = eQTLs[cistrans == "cis",]
  eQTLs = eQTLs[!is.na(genesymbol) & genesymbol != "",]
  eQTLs[,dumID := paste(genesymbol,snps,sep="__")]
  eQTLs = eQTLs[!duplicated(dumID),]
  
  res = data.table(phenotype = myRow$pheno,
                   markername = c(eQTLs$snps,genes$markername),
                   genes = c(eQTLs$genesymbol,genes$genename),
                   source = c(rep("eQTL",dim(eQTLs)[1]),rep("proxGene",dim(genes)[1])))
  matched = match(res$markername,check$rsID)
  res[,region:=check[matched,region]]
  res[,leadSNP:=check[matched,leadSNP]]
  res[,dumID := paste(genes,markername, sep="__")]
  table(duplicated(res$dumID))
  dups = res[duplicated(dumID),dumID]
  res[dumID %in% dups,source := "eQTL and proxGene"]
  res = res[!duplicated(dumID),]
  res
}
tab7 = rbindlist(tab7)
table(tab7$phenotype)
tab7[markername!=leadSNP,table(region)]

tab7[genes == "NGFRAP1", genes:="Z92846.1"]
tab7[genes == "WBP5", genes:="TCEAL9"]
tab7[genes == "FAM127A", genes:="TMEM35A"]
tab7[genes == "FAM58A", genes:="CCNQ"]
tab7[genes == "KAL1", genes:="ANOS1"]
tab7[genes == "FAM46D", genes:="TENT5D"]
tab7[genes == "CXorf57", genes:="RADX"]
tab7[genes == "KCNE1L", genes:="KCNE5"]
tab7[genes == "SEPT6", genes:="SEPTIN6"]
tab7[genes == "TMEM35", genes:="TMEM35A"]
tab7[genes == "RGAG4", genes:="RTL5"]
tab7[genes == "FAM122C", genes:="PABIR3"]
tab7[genes == "BHLHB9", genes:="GPRASP3"]
tab7[genes == "FAM122B", genes:=" PABIR2"]

candidateGenes = unique(tab7$genes)
candidateGenes = candidateGenes[candidateGenes!=""]
myGenTab<-data.table(genename=candidateGenes)

genes38 = fread("../temp/07_HGNC_Download_221124.txt")
table(is.element(myGenTab$genename, genes38$symbol))
myGenTab[!is.element(genename,genes38$symbol),]
myGenTab = myGenTab[is.element(genename,genes38$symbol),]
m1 <- match(myGenTab$gene, genes38$symbol)
genes38 = genes38[m1,]

myGenTab[, `:=`(
  ensg = genes38[, ensembl_gene_id],
  entrez = genes38[, entrez_id],
  hgnc = genes38[, hgnc_id],
  description = genes38[, name],
  type = genes38[,locus_group],
  cytoband = genes38[,location_sortable ]
)]
myGenTab

table(is.na(myGenTab$ensg))
table(is.na(myGenTab$entrez))
setorder(myGenTab,cytoband)

#' # Check data ####
#' ***
#' ## Check 1: nsnps ###
#' ***
#' number of used SNPs should be >100
check1 = coloc[,min(nsnps)]
check1
coloc[nsnps == min(nsnps)]
hist(coloc$nsnps)
hist(coloc[grepl("Kidney_Cortex_Glomerular",trait2) | grepl("Kidney_Cortex_Tubulointerstitial",trait2),nsnps],
     main = "Histogram NEPTUNE",
     xlab = "NEPTUNE overlap")
hist(coloc[!grepl("Kidney_Cortex_Glomerular",trait2) & !grepl("Kidney_Cortex_Tubulointerstitial",trait2),nsnps],
     main = "Histogram GTEx",
     xlab = "GTEx overlap")

#' ## Check 2: PP.H0 ####
#' ***
#' should not be too high 
check2 = copy(coloc)
check2 = check2[PP.H0.abf>=0.75,]
check2[,dumID := paste(region, trait1,gene,sep="::")]
check2 = check2[!duplicated(dumID),]

#' What could go wrong here? 
#' 
#' * GTEx: range +/- 1 mb around gene
#' * NEPTUNE: range +/- 500 kb around gene
#' * our selection: index SNP +/- 250 kb or eQTL (max dist 1 mb)
#' 
#' * Option 1: index SNP not in eQTL data --> regions overlap, but no signal for CKD traits (might be the case for NEPTUNE data as there are less SNPs available (filtered for MAF>3% per default and lower range than GTEx and our eQTL definition))
#' * Option 2: sex-specific association (check with Main Table 1)
#' 
dummy = check2[trait1 == "eGFR_ALL"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_ALL",]

#' **eGFR_ALL**: all outlier are in NEPTUNE data (see option 1) or not for the ALL phenotype
#' 
dummy = check2[trait1 == "eGFR_MALE"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_MALE",]

#' **eGFR_MALE**: all outlier are in NEPTUNE data (see option 1) or not for the MALE phenotype
#' 
dummy = check2[trait1 == "eGFR_FEMALE",]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_FEMALE",]

#' **eGFR_FEMALE**: all outlier are in NEPTUNE data (see option 1) or not for the FEMALE phenotype
#' 
dummy = check2[trait1 == "UA_ALL"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_ALL",]

data_GWAS = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
load("../temp/07_coloc/GTEx_v8_filtered_Adipose_Subcutaneous.RData")
data1 = copy(data0)
data1 = data1[cyto == "Xq13.1",]
data_GWAS = data_GWAS[P<1e-6,]
tab6[phenotype=="UA_ALL" & region == 17,]
data_GWAS = data_GWAS[position <= tab6[phenotype=="UA_ALL" & region == 17,bp]+500000,]
data_GWAS = data_GWAS[position >= tab6[phenotype=="UA_ALL" & region == 17,bp]-500000,]
data_GWAS
table(is.element(data_GWAS$position,data1$pos_b37))

#' **UA_ALL**: 
#' 
#' * Region 22: outlier is in NEPTUNE data (see option 1)
#' * Region 17: Index SNP & the two other associated SNPs with suggestive significance are not available in GTEx. I check the GTEx Portal and searched for the rsID but found no matches. 
#' 
dummy = check2[trait1 == "UA_MALE"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_MALE",]

#' **UA_MALE**: all outlier are in NEPTUNE data (see option 1) or not for the MALE phenotype
#' 
dummy = check2[trait1 == "UA_FEMALE",]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_FEMALE",]

#' **UA_FEMALE**: all outlier are in NEPTUNE data (see option 1) or not for the FEMALE phenotype
#' 
#' ## Check 3: PP.H2 ####
#' ***
#' should not be too high 
check3 = copy(coloc)
check3 = check3[PP.H2.abf>=0.75,]
check3[,dumID := paste(region, trait1,gene,sep="::")]
check3 = check3[!duplicated(dumID),]

dummy = check3[trait1 == "eGFR_ALL"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_ALL",]

#' **eGFR_ALL**: no problems
#' 
dummy = check3[trait1 == "eGFR_MALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_MALE",]

#' **eGFR_MALE**: no problems
#' 
dummy = check3[trait1 == "eGFR_FEMALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_FEMALE",]

#' **eGFR_FEMALE**: no problems
#' 
dummy = check3[trait1 == "UA_ALL"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_ALL",]

#' **UA_ALL**: Region 17: Index SNPs not available
#' 
dummy = check3[trait1 == "UA_MALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_MALE",]

#' **UA_MALE**: no problems
#' 
dummy = check3[trait1 == "UA_FEMALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_FEMALE",]

#' **UA_FEMALE**: no problems
#' 
#' ## Check 4: Coloc ####
#' ***
#' 
check4 = copy(coloc)
check4 = check4[PP.H4.abf>=0.75,]
length(unique(check4$gene))
length(unique(check4$region))

check4[,dumID := paste(gene,trait1,sep="::")]

dumTab6<-dcast(check4,
               formula = dumID ~ trait2,
               value.var = c("PP.H4.abf"),
               sep = "_")
names(dumTab6)
dumTab6[!is.na(`GE in Kidney_Cortex_Tubulointerstitial `),c(1,28)]
dumTab6[!is.na(`GE in Muscle_Skeletal`),c(1,31)]
dumTab6[!is.na(`GE in Whole_Blood`),c(1,45)]

candidates = dumTab6[!is.na(`GE in Kidney_Cortex_Tubulointerstitial `) | !is.na(`GE in Muscle_Skeletal`) | !is.na(`GE in Whole_Blood`) , dumID]
candidates = gsub("::.*","",candidates)
candidates = unique(candidates)
res = tab7[genes %in% candidates,]
res

candidates2 = dumTab6[, dumID]
candidates2 = gsub("::.*","",candidates2)
candidates2 = unique(candidates2)
res2 = tab7[genes %in% candidates2,]
res2

#' # Summary ####
#' ***
#' 
#' **SNP numbers**: all good
#' 
#' **PP H0, PP H2**: all good or explainable
#' 
#' **Co-localization**:
#' 
#' Best candidates in relevant tissues (kidney, muscle, blood):  
#' 
#' * ARMCX2 (UA and eGFR in ALL and MALE) - region 7 and 18
#' * CDKL5 (eGFR in ALL, MALE, and FEMALE) - region 2
#' * SLC25A5 (eGFR in ALL and MALE) - region 11
#' * NDUFB11 (eGFR in ALL) - region 3
#' * TCEAL3 (eGFR and UA in MALE) - region 8 and 19
#' * ACSL4 (eGFR in ALL, MALE and FEMALE) - region 10
#' * MORF4L2 (eGFR and UA in MALE) - region 8 and 19
#' 
#' In total, there are 38 unique genes with positive coloc in at least one tissue, corresponding to 14 of the 22 regions. 
#' 
#' # Plot ####
#' *** 
#' To be discussed...
#' 
#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

 