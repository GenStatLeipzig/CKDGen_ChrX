#' ---
#' title: "Co-localization Part 1: get eQTL data"
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

source("../SourceFile_angmar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Step 1: Get genes per region ####
#' ***
#' First I load all independent signals, as those define the relevant phenotypes.
#' 
#' Second, I load the hit annotations per phenotype and extract all cis-eQTL genes (max dist 1 MB) and all proximate genes (max dist 250 kb). 
#' 
#' Then, I add the ensembl ID, which is used in the GTEx data, and the entrez ID used in NEPTUNE. For this, I use the file downloaded from https://www.genenames.org/download/statistics-and-files/ filtered for Chr X (Total Approved Symbols, N=1941, 24.11.2022)
#' 
#' ## Step 1.1: get indep signals ####
#' 
#' * rsID == indep signal
#' * lead SNP == SNP with lowest p-value in region
#' 
tab6 = fread("../results/05_b_Cojo_Select_Results.txt")
tab6[,phenotype := gsub("_all","_ALL",phenotype)]
tab6[,phenotype := gsub("_male","_MALE",phenotype)]
tab6[,phenotype := gsub("_female","_FEMALE",phenotype)]

tab6[,leadSNP := rsID]
tab6[grepl("rs111884516:",rsID),leadSNP := "rs4328011:152898261:G:A"]
tab6[grepl("rs7056552:",rsID),leadSNP := "rs202138804:133799101:AGT:A"]
tab6[grepl("rs111410539:",rsID),leadSNP := "rs181497961:106168067:G:A"]

#' ## Step 1.2: get annotation data ####
#' 
ToDoList3 = data.table(pheno = c("eGFR_ALL","eGFR_FEMALE","eGFR_MALE","UA_ALL","UA_MALE"),
                       genelist = c(paste0(path_CS_eGFR_ALL,"proximate_genes_2022-11-29_credSets.txt"),
                                    paste0(path_CS_eGFR_FEMALE,"proximate_genes_2022-11-29_credSets_female.txt"),
                                    paste0(path_CS_eGFR_MALE,"proximate_genes_2022-11-29_credSets_male.txt"),
                                    paste0(path_CS_UA_ALL,"proximate_genes_2022-11-29_credSets.txt"),
                                    paste0(path_CS_UA_MALE,"proximate_genes_2022-11-29_credSets_male.txt")),
                       eQTLlist = c(paste0(path_CS_eGFR_ALL,"eqtlinfo_2022-11-29_credSets.txt"),
                                    paste0(path_CS_eGFR_FEMALE,"eqtlinfo_2022-11-29_credSets_female.txt"),
                                    paste0(path_CS_eGFR_MALE,"eqtlinfo_2022-11-29_credSets_male.txt"),
                                    paste0(path_CS_UA_ALL,"eqtlinfo_2022-11-29_credSets.txt"),
                                    paste0(path_CS_UA_MALE,"eqtlinfo_2022-11-29_credSets_male.txt")))

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

#' ## Step 1.3: get ENSG ####
#' 
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

#' # Get eQTLs ####
#' ***
#' 
#' ## GTEx v8 ####
#' 
myeQTLs<-dir(path = path_GTExv8, pattern = "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8")
#myeQTLs = myeQTLs[c(3,30,34,49)]

gtex8annot <- fread(paste0(path_GTExv8,"GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"))
gtex8annot

dumTab1<-foreach(i=c(1:length(myeQTLs)))%do%{
  #i=1
  myTissue<-gsub(".allpairs.txt.gz","",myeQTLs[i])
  myTissue<-gsub("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_","",myTissue)
  time1 = Sys.time()

  # Step 1: Loading data
  message("Loading eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  data0<-fread(paste0(path_GTExv8, myeQTLs[i]))
  time2 = Sys.time()
  message("          Finished loading eQTL data in ",round(difftime(time2, time1, tz,units = "min"),2)," minutes")

  # Step 2: Change GeneID and filter for candidate genes
  data0 = data0[grepl("chrX_",variant_id),]
  data0[,ENSG:=gsub(gene_id, pattern = "\\..*", replacement = "")]
  # data1 = copy(data0)
  # data0 = copy(data1)
  message("          Filtering eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  filt<-is.element(data0$ENSG,myGenTab$ensg)
  data0<-data0[filt,]
  time3 = Sys.time()
  message("          Finished filtering eQTL data in ",round(difftime(time3, time2, tz,units = "min"),2)," minutes")

  # Step 3: Get chr, pos and alleles
  message("          Harmonizing column names for ",myTissue," number ",i," of ",length(myeQTLs))
  dummy<-unlist(strsplit(data0$variant_id,"_"))
  chr<-dummy[seq(1,length(dummy),by=5)]
  chr<-gsub("chr","",chr)
  pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  other_allele<-dummy[seq(3,length(dummy),by=5)]
  effect_allele<-dummy[seq(4,length(dummy),by=5)]

  data0[,chr:=chr,]
  data0[,pos_b38:=pos,]
  data0[,other_allele:=other_allele,]
  data0[,effect_allele:=effect_allele,]
  data0[,n_samples:=round(ma_count/maf/2,1)]
  setnames(data0,"slope_se","se")
  setnames(data0,"pval_nominal","pval")
  setnames(data0,"slope","beta")

  # Step 4: Get SNP ID of hg19
  data0[,variant_id_b38:= variant_id]
  data0[,variant_id:= NULL]
  matched<-match(data0$variant_id_b38,gtex8annot$variant_id)
  table(is.na(matched))
  table(data0$variant_id_b38==gtex8annot$variant_id[matched])
  data0$variant_id_b37 <- gtex8annot[matched, variant_id_b37]
  table(is.na(data0$variant_id_b37))
  data0<-data0[!is.na(variant_id_b37)]

  dummy<-unlist(strsplit(data0$variant_id_b37,"_"))
  chr_b37<-dummy[seq(1,length(dummy),by=5)]
  filt = chr_b37 == "X"
  data0<-data0[filt,]
  dummy<-unlist(strsplit(data0$variant_id_b37,"_"))
  pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  data0[,pos_b37:=pos,]
  data0[,chrPos_b37 := paste(chr,pos_b37,sep=":")]
  data0[,chrPos_b37 := paste0("chr",chrPos_b37)]
  data0[,chrPos_b38 := paste(chr,pos_b38,sep=":")]
  data0[,chrPos_b38 := paste0("chr",chrPos_b38)]

  # Step 5: Add gene name and cytoband
  matched<-match(data0$ENSG,myGenTab$ensg)
  data0[,gene:=myGenTab[matched,genename]]
  data0[,cyto:=myGenTab[matched,cytoband]]
  data0[,tissue:=myTissue]
  data0[,SNPGene:=paste(chrPos_b37,gene,sep=":")]

  # Step 6: Check for duplicates (tri allelic positions) and filter for MAF>=1%
  data0[,table(duplicated(SNPGene),duplicated(chrPos_b37))]
  dups = data0[duplicated(data0$SNPGene),]
  data0 = data0[!is.element(SNPGene,dups$SNPGene),]
  data0 = data0[maf>=0.01,]
  x1 = dim(data0)[1]
  message("          Saving ",x1," filtered SNPs")

  # Step 7: Save eqtl data
  outfn1<-paste0("../temp/07_coloc/GTEx_v8_filtered_",myTissue,".RData")
  save(data0,file=outfn1)
  # load(outfn1)
  
  # Step 8: Return something
  dummy = data0[,.N,by=gene]
  res = data.table(tissue = myTissue,
                   genes = dummy$gene,
                   n_eQTLs = dummy$N)
  res
}
dumTab1 = rbindlist(dumTab1)
dumTab1
dumTab2<-dcast(dumTab1,
                   formula = genes ~ tissue,
                   value.var = c("n_eQTLs"),
                   sep = "_")
matched = match(myGenTab$genename,dumTab2$genes)
myGenTab = cbind(myGenTab,dumTab2[matched,])
myGenTab[,genes := NULL]

save(myGenTab,file="../results/07_a_usedGenes.RData")
load("../results/07_a_usedGenes.RData")

#' ## NEPTUNE ####
#' 
myeQTLs<-dir(path = path_NephQTL, pattern = ".cis.gz")
myeQTLs = myeQTLs[!grepl(".tbi",myeQTLs)]

dumTab1<-foreach(i=c(1:length(myeQTLs)))%do%{
  #i=1
  if(grepl("peerEQTL_G",myeQTLs[i])){
    myTissue = "Kidney_Cortex_Glomerular"
    mySampleSize = 136
  }else{
    myTissue = "Kidney_Cortex_Tubulointerstitial "
    mySampleSize = 166
  }
  time1 = Sys.time()
  
  # Step 1: Loading data
  message("Loading eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  data0<-fread(paste0(path_NephQTL,myeQTLs[i]))
  time2 = Sys.time()
  message("          Finished loading eQTL data in ",round(difftime(time2, time1, tz,units = "min"),2)," minutes")

  # Step 2: Change column names
  dummy = names(data0)[1]
  dummy = unlist(strsplit(dummy," "))
  dummy[1] = "CHR"
  names(data0) = dummy

  # Step 3: Filter for chromosome X and relevant genes
  data0 = data0[CHR == "X",]
  message("          Filtering eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  filt<-is.element(data0$ENTREZ_ID,myGenTab$entrez)
  data0<-data0[filt,]
  time3 = Sys.time()
  message("          Finished filtering eQTL data in ",round(difftime(time3, time2, tz,units = "min"),2)," minutes")
  
  # Step 4: Get necessary columns
  message("          Harmonizing column names for ",myTissue," number ",i," of ",length(myeQTLs))
  matched = match(data0$ENTREZ_ID,myGenTab$entrez)
  data0[,ENSG:=myGenTab[matched,ensg]]
  data0[,chrPos_b37 := paste(CHR,POS,sep=":")]
  setnames(data0,"CHR","chr")
  setnames(data0,"POS","pos_b37")
  setnames(data0,"BETA","beta")
  setnames(data0,"SE","se")
  setnames(data0,"P","pval")
  setnames(data0,"REF","other_allele")
  setnames(data0,"ALT","effect_allele")
  data0[,chrPos_b37 := paste0("chr",chr,":",pos_b37)]
  data0[,n_samples := mySampleSize]

  # Step 5: Add gene name and cytoband
  matched<-match(data0$ENSG,myGenTab$ensg)
  data0[,gene:=myGenTab[matched,genename]]
  data0[,cyto:=myGenTab[matched,cytoband]]
  data0[,tissue:=myTissue]
  data0[,SNPGene:=paste(chrPos_b37,gene,sep=":")]

  # Step 6: Check for duplicates (tri allelic positions) and filter for MAF>=1%
  data0[,table(duplicated(SNPGene),duplicated(chrPos_b37))]
  dups = data0[duplicated(data0$SNPGene),]
  data0 = data0[!is.element(SNPGene,dups$SNPGene),]
  x1 = dim(data0)[1]
  message("          Saving ",x1," filtered SNPs")

  # Step 7: Save eqtl data
  outfn1<-paste0("../temp/07_coloc/NEPTUNE_filtered_",myTissue,".RData")
  save(data0,file=outfn1)
  load(outfn1)
  
  # Step 8: Return something
  dummy = data0[,.N,by=gene]
  res = data.table(tissue = myTissue,
                   genes = dummy$gene,
                   n_eQTLs = dummy$N)
  res
}
dumTab1 = rbindlist(dumTab1)
dumTab1
dumTab2<-dcast(dumTab1,
               formula = genes ~ tissue,
               value.var = c("n_eQTLs"),
               sep = "_")
matched = match(myGenTab$genename,dumTab2$genes)
myGenTab = cbind(myGenTab,dumTab2[matched,])
myGenTab[,genes := NULL]

x = dim(myGenTab)[2]
mySums = rowSums(myGenTab[,c(9:x),with=F], na.rm=TRUE)
table(mySums == 0)
myGenTab = myGenTab[mySums != 0,]

save(myGenTab,file="../results/07_a_usedGenes.RData")
save(tab7,file = "../temp/07_allGenes.RData")

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

