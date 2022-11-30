#' ---
#' title: "Heritability and Genetic Correlation"
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
#' We want to estimate the heritability for eGFR and UA on chromosome X and compare it to the autosomale heritability of previous CKDGen publications. 
#' 
#' In a second step, we want to test for genetic correlation between eGFR and UA on chromosome X. This will only be done when both traits provide significant heritability estimates. 
#' 
#' For both, we use the GCTA software and UK Biobank genetic & phenotypic data (largest participating study).  
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Step 1: create PLINK data ####
#' ***
#' We want to use genotyped SNPs only (to save computational time - can be done with complete data set later). In addition, we only keep the samples with phenotype information. 
#' 
snp.QC = fread(paste0(path_UKBB,"2010_data/ukb_mfi_chrX_v3.txt"))
snp.QC
table(snp.QC[, V6] > 0.005, snp.QC[, V8] == 1, dnn = c("MAF > 0.5%", "genotyped (or perfect information)"))
snp.QC = snp.QC[V6 > 0.005 & V8 == 1, ]
goodSNP_fn = "../temp/10_GCTA_Heritability/genotypedSNPs.txt"
write.table(snp.QC[, V2], file = goodSNP_fn, col.names = F, row.names = F, quote = F)

myPhenoTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Phenotypes.txt"))
myIDTab_Mal = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Male_exclude.txt"))
myIDTab_Fem = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Female_exclude.txt"))
myPhenoTab = myPhenoTab[!is.na(eGFR_overall) | !is.na(uric_acid_overall),]
myPhenoTab = myPhenoTab[,c(1,2)]
myPhenoTab_males = copy(myPhenoTab)
myPhenoTab_males = myPhenoTab_males[FID %in% myIDTab_Mal$V1,]
myPhenoTab_females = copy(myPhenoTab)
myPhenoTab_females = myPhenoTab_females[FID %in% myIDTab_Fem$V1,]

maleSamples_fn = "../temp/10_GCTA_Heritability/maleSamples.txt"
femaleSamples_fn = "../temp/10_GCTA_Heritability/femaleSamples.txt"
write.table(myPhenoTab_males, maleSamples_fn, col.names = T, row.names = F, quote = F, sep = " ")
write.table(myPhenoTab_females, femaleSamples_fn, col.names = T, row.names = F, quote = F, sep = " ")

bgenFile = paste0(path_UKBB,"2010_data/ukb_imp_chrX_v3.bgen")
sampleFile = paste0(path_UKBB,"2010_data/ukb20272_imp_chrX_v3_s486631.sample")
bedFile_genotyped_m = "../temp/10_GCTA_Heritability/ukb20272_imp_chrX_v3_s486631_genotyped_males"
bedFile_genotyped_f = "../temp/10_GCTA_Heritability/ukb20272_imp_chrX_v3_s486631_genotyped_females"

mycall0_m = paste0(path_plink2, " --bgen ", bgenFile, 
                   " ref-first --sample ", sampleFile, 
                   " --oxford-single-chr 23 --make-bed --hard-call-threshold 0.49",
                   " --extract ", goodSNP_fn,
                   " --keep ",maleSamples_fn,
                   " --out ", bedFile_genotyped_m)
mycall0_m
#system(mycall0_m)

mycall0_f = paste0(path_plink2, " --bgen ", bgenFile, 
                   " ref-first --sample ", sampleFile, 
                   " --oxford-single-chr 23 --make-bed --hard-call-threshold 0.49",
                   " --extract ", goodSNP_fn,
                   " --keep ",femaleSamples_fn,
                   " --out ", bedFile_genotyped_f)
mycall0_f
#system(mycall0_f)

#' # Step 1: create GRM ####
#' ***
#' Estimate the GRM from the SNPs on the X-chromosome. Given the large data file, we will split this in three steps as suggested by the GCTA webpage:
#' 
#' gcta64 --bfile test  --make-grm-xchr  --out test_xchr
#' 
#' gcta64 --bfile test --make-grm-part 3 1 --thread-num 5 --out test
#' 
#' gcta64 --bfile test --make-grm-part 3 2 --thread-num 5 --out test
#' 
#' gcta64 --bfile test --make-grm-part 3 3 --thread-num 5 --out test

path_GRM_m = "../temp/10_GCTA_Heritability/ukb20272_imp_chrX_v3_s486631_genotyped_GRM_m"
path_GRM_f = "../temp/10_GCTA_Heritability/ukb20272_imp_chrX_v3_s486631_genotyped_GRM_f"

#' ## Partition the GRM into 3 parts - males ####
mycall1<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_m,
                " --make-grm-xchr-part 3 1 --thread-num 5",
                " --out ",path_GRM_m)
mycall1
#system(mycall1)

mycall2<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_m,
                " --make-grm-xchr-part 3 2 --thread-num 5",
                " --out ",path_GRM_m)
mycall2
#system(mycall2)

mycall3<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_m,
                " --make-grm-xchr-part 3 3 --thread-num 5",
                " --out ",path_GRM_m)
mycall3
#system(mycall3)

#' ## Merge all the parts together - males #### 
mycall4 = paste0("cat ",path_GRM_m,".part_3_*.grm.id > ",path_GRM_m,".grm.id")
mycall5 = paste0("cat ",path_GRM_m,".part_3_*.grm.bin > ",path_GRM_m,".grm.bin")
mycall6 = paste0("cat ",path_GRM_m,".part_3_*.grm.N.bin > ",path_GRM_m,".grm.N.bin")

# system(mycall4)
# system(mycall5)
# system(mycall6)

#' ## Partition the GRM into 3 parts - females ####
mycall1<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_f,
                " --make-grm-xchr-part 3 1 --thread-num 5",
                " --out ",path_GRM_f)
mycall1
#system(mycall1)

mycall2<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_f,
                " --make-grm-xchr-part 3 2 --thread-num 5",
                " --out ",path_GRM_f)
mycall2
#system(mycall2)

mycall3<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped_f,
                " --make-grm-xchr-part 3 3 --thread-num 5",
                " --out ",path_GRM_f)
mycall3
#system(mycall3)

#' ## Merge all the parts together - males #### 
mycall4 = paste0("cat ",path_GRM_f,".part_3_*.grm.id > ",path_GRM_f,".grm.id")
mycall5 = paste0("cat ",path_GRM_f,".part_3_*.grm.bin > ",path_GRM_f,".grm.bin")
mycall6 = paste0("cat ",path_GRM_f,".part_3_*.grm.N.bin > ",path_GRM_f,".grm.N.bin")

# system(mycall4)
# system(mycall5)
# system(mycall6)


#' # Step 2: get necessary phenotype files ####
#' ***
myPhenoTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Phenotypes.txt"))
myCovarTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Covariates.txt"))
myIDTab_Mal = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Male_exclude.txt"))
myIDTab_Fem = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Female_exclude.txt"))

#' ## create .phen file
#' columns are family ID, individual ID and phenotypes
#' 
#' FID, IID, eGFR_overall, uric_acid_overall
names(myPhenoTab)
myPhenoTab = myPhenoTab[,c(1,2,3,12)]
write.table(myPhenoTab, "../temp/10_GCTA_Heritability/10_UKBB.phen", row.names=F, quote=F, col.names=F)

#' ## create .covar file
#' columns are family ID, individual ID and discrete covariates
#' 
#' FID, IID, sex
myCovarTab1 = copy(myCovarTab)
myCovarTab1 = myCovarTab1[,c(1,2)]
myCovarTab1[,sex := ""]
myCovarTab1[FID %in% myIDTab_Fem$V1,sex :="F"]
myCovarTab1[FID %in% myIDTab_Mal$V1,sex :="M"]
myCovarTab1[1:20000,table(sex)]
write.table(myCovarTab1, "../temp/10_GCTA_Heritability/10_UKBB.covar", row.names=F, quote=F, col.names=F)

#' ## create .qcovar file
#' columns are family ID, individual ID and quantitative covariates
#' 
#' FID, IID, age and PCs
myCovarTab2 = copy(myCovarTab)
names(myCovarTab2)
write.table(myCovarTab1, "../temp/10_GCTA_Heritability/10_UKBB.qcovar", row.names=F, quote=F, col.names=F)

#' # Step 3: estimate heritability ####
#' ***
#' gcta64  --reml  --grm test  --pheno test.phen  --grm-adj 0  --grm-cutoff 0.05  --out test
mycall7_m<-paste0(path_gcta, 
                " --reml "," --grm ",path_GRM_m,
                " --pheno ../temp/10_GCTA_Heritability/10_UKBB.phen --mpheno 1",
                " --qcovar ../temp/10_GCTA_Heritability/10_UKBB.qcovar",
                " --out ../results/10_eGFR_heritab_m")
mycall7_m
#system(mycall7_m)

mycall7_f<-paste0(path_gcta, 
                  " --reml "," --grm ",path_GRM_f,
                  " --thread-num 5",
                  " --pheno ../temp/10_GCTA_Heritability/10_UKBB.phen --mpheno 1",
                  " --qcovar ../temp/10_GCTA_Heritability/10_UKBB.qcovar",
                  " --out ../results/10_eGFR_heritab_f")
mycall7_f
#system(mycall7_f)

mycall8_m<-paste0(path_gcta, 
                " --reml "," --grm ",path_GRM_m,
                " --pheno ../temp/10_GCTA_Heritability/10_UKBB.phen --mpheno 2",
                " --qcovar ../temp/10_GCTA_Heritability/10_UKBB.qcovar",
                " --out ../results/10_UA_heritab_m")
mycall8_m
system(mycall8_m)

mycall8_f<-paste0(path_gcta, 
                  " --reml "," --grm ",path_GRM_f,
                  " --pheno ../temp/10_GCTA_Heritability/10_UKBB.phen --mpheno 2",
                  " --qcovar ../temp/10_GCTA_Heritability/10_UKBB.qcovar",
                  " --out ../results/10_UA_heritab_f")
mycall8_f
system(mycall8_f)


#' # Step 4: load results ####
#' ***
eGFR_m<-read.table("../results/10_eGFR_heritab_m.hsq",nrows=4,header=T)
eGFR_f<-read.table("../results/10_eGFR_heritab_f.hsq",nrows=4,header=T)
UA_m<-read.table("../results/10_UA_heritab_m.hsq",nrows=4,header=T)
UA_f<-read.table("../results/10_UA_heritab_f.hsq",nrows=4,header=T)

myTab = data.table(phenotype= c("eGFR_MALE","eGFR_FEMALE","UA_MALE","UA_FEMALE"),
                   h2 = c(eGFR_m$Variance[4],eGFR_f$Variance[4],UA_m$Variance[4],UA_f$Variance[4]),
                   SE = c(eGFR_m$SE[4],eGFR_f$SE[4],UA_m$SE[4],UA_f$SE[4]),
                   pval = c(eGFR_m$Variance[9],eGFR_f$Variance[9],UA_m$Variance[9],UA_f$Variance[9]),
                   N = c(eGFR_m$Variance[10],eGFR_f$Variance[10],UA_m$Variance[10],UA_f$Variance[10]))
myTab

save(myTab,file = "../results/10_Heritability.RData")
write.table(myTab, file = "../results/10_Heritability.txt",col.names =T, row.names = F, quote = F,sep="\t")


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
