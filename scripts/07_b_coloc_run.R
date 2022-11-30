#' ---
#' title: "Co-localization Part 2: Run coloc"
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
source("../helperFunctions/colocFunction_jp.R")
source("../helperFunctions/wait4finishedOutfile_hk.R")

#' # Wait for part a ####
#' ***
wait4finishedOutfile(outfile_fn = "07_a_coloc_get_eQTL_data.R.out",
                     waitingtime_seconds = 10,
                     erkennungsmarke = "TOTAL TIME :")

#' # Prep data ####
#' ***
load("../results/07_a_usedGenes.RData")
load("../temp/07_allGenes.RData")

#' # Run Coloc ####
#' ***
#' Loop 1: CKDGen Phenotype (k=1-6 for eGFR and UA in ALL, MALE and FEMALE)
#' 
#' Loop 2: Tissue (i=1:51 for 49 GTEx tissues and 2 NEPTUNE tissues)
#' 
#' Loop 3: Genes (all available gene per tissue, independent of top phenotype?)
#' 
ToDoList = copy(myGenTab)
filt1 = is.element(myGenTab$genename,tab7[grepl("eGFR",phenotype),genes])
table(filt1)
filt2 = is.element(myGenTab$genename,tab7[grepl("UA",phenotype),genes])
table(filt2)
ToDoList[,eGFR :=F]
ToDoList[filt1,eGFR :=T]
ToDoList[,UA :=F]
ToDoList[filt2,UA :=T]
table(ToDoList$UA,ToDoList$eGFR)

myPhenos = c("eGFR_ALL","eGFR_MALE","eGFR_FEMALE","UA_ALL","UA_MALE","UA_FEMALE")

registerDoMC(cores=20)

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
  
  myeQTLs2<-dir(path = "../temp/07_coloc/",pattern = ".RData")
  myeQTLs2 = myeQTLs2[grepl("GTEx_v8_filtered_",myeQTLs2) | grepl("NEPTUNE_filtered_",myeQTLs2)]
  
  dumTab2 = foreach(i=c(1:length(myeQTLs2)))%dopar%{
    #i=5
    loaded = load(paste0("../temp/07_coloc/",myeQTLs2[i]))
    data_eQTLs = get(loaded)
    myTissue2 = gsub("GTEx_v8_filtered_","",myeQTLs2[i])
    myTissue2 = gsub("NEPTUNE_filtered_","",myTissue2)
    myTissue2 = gsub(".RData","",myTissue2)
    message("Working on tissue ",myTissue2)
    
    # get To Do list
    dummy = data_eQTLs[,.N,by=gene]
    ToDoList2 = copy(ToDoList)
    if(myPheno2 == "eGFR"){
      ToDoList2 = ToDoList2[eGFR==T,]
    }else{
      ToDoList2 = ToDoList2[UA==T,]
    }
    dummy = dummy[gene %in% ToDoList2[,genename]]
    
    dumTab3 = foreach(j=c(1:dim(dummy)[1]))%do%{
      #j=1
      myRow = dummy[j,]
      moreInfo = copy(ToDoList2)
      moreInfo = moreInfo[genename == myRow$gene,]
      
      data_eQTL1 = copy(data_eQTLs)
      data_eQTL1 = data_eQTL1[gene ==myRow$gene, ]
      dup2 = data_eQTL1[duplicated(chrPos_b37),]
      data_eQTL1 = data_eQTL1[!is.element(chrPos_b37,dup2$chrPos_b37),]
      
      data_GWAS1 = copy(data_GWAS)
      data_GWAS1 = data_GWAS1[chrPos_b37 %in% data_eQTL1$chrPos_b37,]
      dups = data_GWAS1[duplicated(chrPos_b37),]
      data_GWAS1 = data_GWAS1[!is.element(chrPos_b37,dups$chrPos_b37),]
      
      data_eQTL1 = data_eQTL1[chrPos_b37 %in% data_GWAS1$chrPos_b37,]
      
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
      x4[,cytoband:=moreInfo$cytoband]
      x4[,gene:= myRow$gene]
      x4[,trait1:=myPheno]
      x4[,trait2:=paste0("GE in ",myTissue2)]
      x4
      
    }
    dumTab3 = rbindlist(dumTab3)
    dumTab3
  }
  ColocTable = rbindlist(dumTab2)
  ColocTable
}

ColocTable = rbindlist(dumTab)

ColocTable[,table(PP.H4.abf>=0.75)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]

#' Add information about regions
#' 
tab7[,dumID := paste(region,genes,sep ="___")]
tab7 = tab7[genes != ""]
test7 = unique(tab7$dumID)
test7 = test7[order(test7)]
table(duplicated(test7))
test8 = unlist(strsplit(test7,"___"))
region = test8[seq(1,length(test8),2)]
table(region)
genes = test8[seq(2,length(test8),2)]
dumTab5 = data.table(dumID = test7, region =region, genes = genes)
dumTab5[,region:=as.numeric(region)]
setorder(dumTab5,region,genes)
dumTab6<-dcast(dumTab5,
               formula = genes ~ region,
               value.var = c("region"),
               sep = "_")
names(dumTab6)

dumTab7 = foreach(m=1:dim(dumTab6)[1])%do%{
  #m=1
  test9 = paste(dumTab6[m,2:23],collapse = ", ")
  test9 = gsub("NA, ","",test9)
  test9 = gsub(", NA","",test9)
  res = data.table(gene = dumTab6[m,1],
                   regions = test9)
  res
}  
dumTab7 = rbindlist(dumTab7)
matched = match(ColocTable$gene,dumTab7$gene)
table(is.na(matched))
ColocTable[,region := dumTab7[matched,regions]]
ColocTable = ColocTable[,c(11,7,9,8,10,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Number of associated region (1-15: eGFR, 16-22: UA); multiple regions are possible due to overlapping eQTLs (LD r2>0.3 between index SNP of region and one eQTL of gene)",
                                        "Genomic cytoband of index SNP",
                                        "Tested CKDGen phenotype",
                                        "Tested gene",
                                        "Tested tissue",
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))



save(ColocTable, description,file="../results/07_b_coloc_eQTLs.RData")
tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/07_b_coloc_eQTLs.xlsx"
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

