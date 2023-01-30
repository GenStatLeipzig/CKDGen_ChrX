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
#' **Important note**: UKBB GRM data was to large for our IMISE intern compute cluster. Hence, we did this analyses on **sirius** (big memory cluster consisting of 2 nodes with a total of 256 cores and 12 terabytes of consumable memory). In this script, all preparations are done
#' 
#' * get PLINK data format for the genotyped SNPs in the combined setting 
#'     * SNP filter: MAF> 0.005 & genotyped or perfect information
#'     * Sample filter: no NA for eGFR or UA
#' * create GRM for the combined setting
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(paste0(projectpath,"scripts/"))

#' # Step 1: create PLINK data ####
#' ***
#' We want to use genotyped SNPs only (to save computational time - can be done with complete data set later). In addition, we only keep the samples with phenotype information. 
#' 
#' ## Get good SNPs ####
snp.QC = fread(paste0(path_UKBB,"2010_data/ukb_mfi_chrX_v3.txt"))
snp.QC
table(snp.QC[, V6] > 0.005, 
      snp.QC[, V8] == 1, 
      dnn = c("MAF > 0.5%", "genotyped (or perfect information)"))
snp.QC = snp.QC[V6 > 0.005 & V8 == 1, ]
goodSNP_fn = "../temp/10_GCTA_Heritability/Step1_PLINK/genotypedSNPs.snplist"
write.table(snp.QC[, V2], file = goodSNP_fn, col.names = F, row.names = F, quote = F)

#' ## Get good samples ####
myPhenoTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Phenotypes.txt"))
myPhenoTab
myPhenoTab[,table(!is.na(eGFR_overall), 
                  !is.na(uric_acid_overall),
                  dnn = c("eGFR","UA"))]
myPhenoTab = myPhenoTab[!is.na(eGFR_overall) | !is.na(uric_acid_overall),]
myPhenoTab = myPhenoTab[,c(1,2)]
goodSamples_fn = "../temp/10_GCTA_Heritability/Step1_PLINK/goodIndiv.list"
write.table(myPhenoTab, goodSamples_fn, col.names = T, row.names = F, quote = F, sep = " ")

#' ## PLINK Call ####
bgenFile = paste0(path_UKBB,"2010_data/ukb_imp_chrX_v3.bgen")
sampleFile = paste0(path_UKBB,"2010_data/ukb20272_imp_chrX_v3_s486631.sample")
bedFile_genotyped = "../temp/10_GCTA_Heritability/Step1_PLINK/ukb20272_imp_chrX_v3_s486631_genotyped"

mycall0 = paste0(path_plink2, " --bgen ", bgenFile, 
                   " ref-first --sample ", sampleFile, 
                   " --oxford-single-chr 23 --make-bed --hard-call-threshold 0.49",
                   " --extract ", goodSNP_fn,
                   " --keep ",goodSamples_fn,
                   " --out ", bedFile_genotyped)
mycall0
#system(mycall0)

#' # Step 2: create GRM ####
#' ***
#' Estimate the GRM from the SNPs on the X-chromosome. Given the large data file, we will split this in three steps as suggested by the GCTA webpage.
#' 
path_GRM = "../temp/10_GCTA_Heritability/Step2_GRM/ukb20272_imp_chrX_v3_s486631_genotyped_GRM"

#' ## Partition the GRM into 3 parts ####
mycall1<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped,
                " --make-grm-xchr-part 3 1 --thread-num 10",
                " --out ",path_GRM)
mycall1
#system(mycall1)

mycall2<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped,
                " --make-grm-xchr-part 3 2 --thread-num 10",
                " --out ",path_GRM)
mycall2
#system(mycall2)

mycall3<-paste0(path_gcta,
                " --bfile ", bedFile_genotyped,
                " --make-grm-xchr-part 3 3 --thread-num 10",
                " --out ",path_GRM)
mycall3
#system(mycall3)

#' ## Merge all the parts together #### 
mycall4 = paste0("cat ",path_GRM,".part_3_*.grm.id > ",path_GRM,".grm.id")
mycall5 = paste0("cat ",path_GRM,".part_3_*.grm.bin > ",path_GRM,".grm.bin")
mycall6 = paste0("cat ",path_GRM,".part_3_*.grm.N.bin > ",path_GRM,".grm.N.bin")

# system(mycall4)
# system(mycall5)
# system(mycall6)

#' # Step 3: create input for reml ####
#' ***
#' We will need three phenotype files:
#' 
#' * .phen File, containing 6 phenotypes (eGFR and UA in combined, males and females)
#' * .covar File, containing the discrete covariate sex
#' * .qcovar File, containing the quantitative covariates age and the first 10 PCs
#' 
#' In addition, I want to generate all GCTA calls here as well and save the in a .sh file to be copied to **sirius**
#' 
#' ## Create phenotype files ####
myPhenoTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Phenotypes.txt"))
myCovarTab = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_Covariates.txt"))
myIDTab_Mal = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Male_exclude.txt"))
myIDTab_Fem = fread(paste0(path_UKBB,"GWAS/UKBB_CKDGen_GWAS_noQC_and_Female_exclude.txt"))

names(myPhenoTab)
myPhenoTab = myPhenoTab[,c(1,2,3,3,3,12,12,12)]
names(myPhenoTab)[3:8] = c("eGFR_ALL","eGFR_MALE","eGFR_FEMALE",
                           "UA_ALL","UA_MALE","UA_FEMALE")

myCovarTab1 = copy(myCovarTab)
myCovarTab1 = myCovarTab1[,c(1,2)]
myCovarTab1[,sex := ""]
myCovarTab1[FID %in% myIDTab_Fem$V1,sex :="F"]
myCovarTab1[FID %in% myIDTab_Mal$V1,sex :="M"]

myPhenoTab[myCovarTab1$sex == "F",eGFR_MALE := NA]
myPhenoTab[myCovarTab1$sex == "F",UA_MALE := NA]
myPhenoTab[myCovarTab1$sex == "M",eGFR_FEMALE := NA]
myPhenoTab[myCovarTab1$sex == "M",UA_FEMALE := NA]

#' Using all samples requires too much storage for GCTA REML analyses. 
#' 
#' Hence, we only ues 200,000 random samples, (100,000 females and males each)
#' 
set.seed(2204)
samplesEGFR_M = myPhenoTab[!is.na(eGFR_MALE),FID]
samplesEGFR_F = myPhenoTab[!is.na(eGFR_FEMALE),FID]
samplesEGFR_M2 = sample(x=samplesEGFR_M,size=100000,replace = F)
samplesEGFR_F2 = sample(x=samplesEGFR_F,size=100000,replace = F)
samplesEGFR_ALL = c(samplesEGFR_M2,samplesEGFR_F2)
myPhenoTab[FID %nin% samplesEGFR_ALL,eGFR_ALL := NA]
table(is.na(myPhenoTab$eGFR_ALL),myCovarTab1$sex)

samplesUA_M = myPhenoTab[!is.na(UA_MALE),FID]
samplesUA_F = myPhenoTab[!is.na(UA_FEMALE),FID]
samplesUA_M2 = sample(x=samplesUA_M,size=100000,replace = F)
samplesUA_F2 = sample(x=samplesUA_F,size=100000,replace = F)
samplesUA_ALL = c(samplesUA_M2,samplesUA_F2)
myPhenoTab[FID %nin% samplesUA_ALL,UA_ALL := NA]
table(is.na(myPhenoTab$UA_ALL),myCovarTab1$sex)

table(is.na(myPhenoTab$UA_ALL),is.na(myPhenoTab$eGFR_ALL))

myCovarTab2 = copy(myCovarTab)
names(myCovarTab2)

write.table(myPhenoTab, "../temp/10_GCTA_Heritability/Step3_InputREML/UKBB.phen", 
            row.names=F, quote=F, col.names=F)
write.table(myCovarTab1, "../temp/10_GCTA_Heritability/Step3_InputREML/UKBB.covar", 
            row.names=F, quote=F, col.names=F)
write.table(myCovarTab2, "../temp/10_GCTA_Heritability/Step3_InputREML/UKBB.qcovar", 
            row.names=F, quote=F, col.names=F)


#' ## Create SLURM scripts ####
#' 
#' 1) eGFR - MALE - no adjustment
#' 2) eGFR - MALE - adjusted for age and PCs
#' 3) UA - MALE - no adjustment
#' 4) UA - MALE - adjusted for age and PCs
#' 5) eGFR - FEMALE - no adjustment
#' 6) eGFR - FEMALE - adjusted for age and PCs
#' 7) UA - FEMALE - no adjustment
#' 8) UA - FEMALE - adjusted for age and PCs
#' 9) eGFR - ALL - adjusted for sex, age, and PCs
#' 10) UA - ALL - adjusted for sex, age, and PCs
#' 
slurmHeader = c("#!/bin/bash",
                "",
                "#SBATCH --time=10-0",
                "#SBATCH --mem=5TB",
                "#SBATCH --ntasks=1",
                "#SBATCH --cpus-per-task=100",
                "#SBATCH --job-name=",
                "#SBATCH --partition=sirius-long",
                "#SBATCH --exclusive",
                "",
                "module load GCTA/1.94.0beta-foss-2021b")

Sirius_pwd = "/work/users/ju423tama/"
GRM_fn = gsub(".*/","",path_GRM)

call = paste0("gcta64 --reml",
               " --grm ",Sirius_pwd,GRM_fn,
               " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 2",
               " --thread-num 100",
               " --out ",Sirius_pwd,"01_eGFR_MALE_noAdj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='01_eGFR_MALE_noAdj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/01_eGFR_MALE_noAdj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
               " --grm ",Sirius_pwd,GRM_fn,
               " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 2",
               " --qcovar ",Sirius_pwd,"UKBB.qcovar",
               " --thread-num 100",
               " --out ",Sirius_pwd,"02_eGFR_MALE_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='02_eGFR_MALE_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/02_eGFR_MALE_Adj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 5",
              " --thread-num 100",
              " --out ",Sirius_pwd,"03_UA_MALE_noAdj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='03_UA_MALE_noAdj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/03_UA_MALE_noAdj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 5",
              " --qcovar ",Sirius_pwd,"UKBB.qcovar",
              " --thread-num 100",
              " --out ",Sirius_pwd,"04_UA_MALE_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='04_UA_MALE_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/04_UA_MALE_Adj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 3",
              " --thread-num 100",
              " --out ",Sirius_pwd,"05_eGFR_FEMALE_noAdj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='05_eGFR_FEMALE_noAdj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/05_eGFR_FEMALE_noAdj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 3",
              " --qcovar ",Sirius_pwd,"UKBB.qcovar",
              " --thread-num 100",
              " --out ",Sirius_pwd,"06_eGFR_FEMALE_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='06_eGFR_FEMALE_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/06_eGFR_FEMALE_Adj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 6",
              " --thread-num 100",
              " --out ",Sirius_pwd,"07_UA_FEMALE_noAdj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='07_UA_FEMALE_noAdj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/07_UA_FEMALE_noAdj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 6",
              " --qcovar ",Sirius_pwd,"UKBB.qcovar",
              " --thread-num 100",
              " --out ",Sirius_pwd,"08_UA_FEMALE_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='08_UA_FEMALE_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/08_UA_FEMALE_Adj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 1",
              " --qcovar ",Sirius_pwd,"UKBB.qcovar",
              " --covar ",Sirius_pwd,"UKBB.covar",
              " --thread-num 100",
              " --out ",Sirius_pwd,"09_eGFR_ALL_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='09_eGFR_ALL_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/09_eGFR_ALL_Adj.sh", 
            row.names=F, quote=F, col.names=F)

call = paste0("gcta64 --reml",
              " --grm ",Sirius_pwd,GRM_fn,
              " --pheno ",Sirius_pwd,"UKBB.phen --mpheno 4",
              " --qcovar ",Sirius_pwd,"UKBB.qcovar",
              " --covar ",Sirius_pwd,"UKBB.covar",
              " --thread-num 100",
              " --out ",Sirius_pwd,"10_UA_ALL_Adj")
shell = c(slurmHeader,call)
shell = gsub("job-name=","job-name='10_UA_ALL_Adj'",shell)
write.table(shell, "../temp/10_GCTA_Heritability/Step4_SlurmScripts/10_UA_ALL_Adj.sh", 
            row.names=F, quote=F, col.names=F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
