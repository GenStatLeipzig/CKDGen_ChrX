#' ---
#' title: "Colocalisation Analysis: eGFR & Uric Acid: Male vs Female"
#' subtitle: "CKDGen Chromosome X"
#' author: "Andreas Kühnapfel"
#' date: "Last compiled on `r format(Sys.time(), "%d %B %Y")`"
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
#' I want to check if these loci co-localise between the sexes (best phenotype was mainly ALL, but coloc does not require genome-wide significant results)
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list = setdiff(ls(), "allinfo"))
time0 = Sys.time()

source("SourceFile_forostar.R")

setwd(projectpath_main)





#####################################################
########## PERFORM COLOCALISATION ANALYSIS ##########
#####################################################





#rm(list = ls())
#rm(list = setdiff(ls(), c("coloc_scripts", "coloc_script")))





##### Setup
#computer  = "angmar"
#computer  = "forostar"
#max_cores = 60





##### Load Packages
#.libPaths(paste0("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/", computer))
#library("readxl")
#library("doParallel")
#library("data.table")
#?getDTthreads
#getDTthreads()
#setDTthreads(1)
#getDTthreads()
#library("coloc")
#library("reshape2")
#library("ggplot2")
#library("ggcorrplot")





##### Common Base Data
locus         = fread("../results/01_Locus_Definitions.txt")
sumStat_1     = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
sumStat_2     = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
sumStat_3     = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
sumStat_4     = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
sumStat_5     = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
sumStat_6     = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
sumStat_union = rbind(sumStat_1, sumStat_2, sumStat_3,
                      sumStat_4, sumStat_5, sumStat_6)





##### Directories
# path_load_step10_1  = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/04_Annotation/gwasresults_1.7_nstud10_V2/"
# path_load_toplist_1 = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/04_Annotation/gwasresults_1.7_nstud10_V2/synopsis/"
path_load_toplist_1 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_stratified/04_Annotation/male/gwasresults_1.7_nstud10/synopsis/"
# path_load_cojo_1    = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/07_cojo/"
# path_load_step10_2  = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_combined/04_Annotation/gwasresults_1.7_NStudies10/"
# path_load_toplist_2 = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_combined/04_Annotation/gwasresults_1.7_NStudies10/synopsis/"
path_load_toplist_2 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_stratified/04_Annotation/female/gwasresults_1.7_nstud10/synopsis/"
# path_load_cojo_2    = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_combined/05_cojo/"
# path_load_loci      = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/_results_combined/"
# path_load_ukbb      = "J:/genstat/02_projekte/1911_CKDgen_chrX/09_meta_FileQC/"
# path_save_coloc     = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_comb_strat/09_coloc/"
path_load_toplist_3 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/04_Annotation/gwasresults_1.7_nstud10_V2/synopsis/"
path_load_toplist_4 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_stratified/04_Annotation/male/gwasresults_1.7_nstud10/synopsis/"
path_load_toplist_5 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_stratified/04_Annotation/female/gwasresults_1.7_nstud10/synopsis/"
path_load_toplist_6 = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/03_uric_acid_allEth_sex_combined/04_Annotation/gwasresults_1.7_NStudies10/synopsis/"





# ##### Files
# step10_1_files  = grep(pattern = "^step10_2_primaryANDsecondary_", dir(path_load_step10_1), value = TRUE)
# step10_2_files  = grep(pattern = "^step10_2_primaryANDsecondary_", dir(path_load_step10_2), value = TRUE)





##### Individual Data
toplist_1 = read_excel(paste0(path_load_toplist_1, "step80_annotated_toplist_validated_2022-01-24_ALL_eGFR_overall_nstud10_male.xlsx"), sheet = "topliste")
toplist_1 = as.data.frame(toplist_1)
toplist_2 = read_excel(paste0(path_load_toplist_2, "step80_annotated_toplist_validated_2022-01-24_ALL_eGFR_overall_nstud10_female.xlsx"), sheet = "topliste")
toplist_2 = as.data.frame(toplist_2)
toplist_3 = read_excel(paste0(path_load_toplist_3, "step80_annotated_toplist_validated_2022-04-26_ALL_eGFR_overall_nstud10_V2.xlsx"), sheet = "topliste")
toplist_3 = as.data.frame(toplist_3)
toplist_4 = read_excel(paste0(path_load_toplist_4, "step80_annotated_toplist_validated_2022-01-24_ALL_uric_acid_overall_nstud10_male.xlsx"), sheet = "topliste")
toplist_4 = as.data.frame(toplist_4)
toplist_5 = read_excel(paste0(path_load_toplist_5, "step80_annotated_toplist_validated_2022-01-24_ALL_uric_acid_overall_nstud10_female.xlsx"), sheet = "topliste")
toplist_5 = as.data.frame(toplist_5)
toplist_6 = read_excel(paste0(path_load_toplist_6, "step80_annotated_toplist_validated_2021-06-29_ALL_uric_acid_overall.xlsx"), sheet = "topliste")
toplist_6 = as.data.frame(toplist_6)

# cojo = read_excel(paste0(path_load_loci, "Von_Wuttke_to_Lokus_220530.xlsx"), sheet = "Independent_SNPs")
# cojo = as.data.frame(cojo)

# loci = read_excel(paste0(path_load_loci, "CKDGen_Region_SNPs_in_ALL_F_M_summary_220608_interaction_20220608.xlsx"), sheet = "Overview")
# loci = as.data.frame(loci)

# ukbb_1 = fread(paste0(path_load_ukbb, "FileQC_ALL_eGFR/step25_ID_Matching_Table.txt.gz"))
# ukbb_2 = fread(paste0(path_load_ukbb, "FileQC_ALL_uric_acid/step25_ID_Matching_Table.txt.gz"))

#dum_1         = (cojo$phenotype == "eGFR") & (cojo$subanalysis == "ALL")
#dum_2         = (cojo$phenotype == "uric acid") & (cojo$subanalysis == "ALL")
#region_list_1 = sort(unique(cojo$Region[dum_1]))
#region_list_2 = sort(unique(cojo$Region[dum_2]))

data = data.frame(snp = locus$rsID,
                  chr = locus$chromosome,
                  pos = locus$position,
                  locus_ll = locus$region_start,
                  locus_ul = locus$region_end,
                  region = locus$region,
                  trait = locus$phenotype)
# for (i in 1 : nrow(data_1)) {
#   
#   dum_i = toplist_1$SNP == data_1$snp[i]
#   
#   #data_1$pos[i] = toplist_1$"position of SNP"[dum_i]
#   
#   gene_i_nearest = toplist_1$"nearest genes (functional relevance)"[dum_i]
#   gene_i_nearest = unlist(strsplit(unlist(strsplit(gene_i_nearest, split = "; ")), split = " \\("))
#   gene_i_nearest = gene_i_nearest[seq(1, length(gene_i_nearest), by = 2)]
#   
#   gene_i_cis = toplist_1$"cis-eQTL genes"[dum_i]
#   gene_i_cis = unlist(strsplit(gene_i_cis, split = " [|] "))
#   
#   gene_i_trans = toplist_1$"trans-eQTL genes"[dum_i]
#   gene_i_trans = unlist(strsplit(gene_i_trans, split = " [|] "))
#   
#   data_1$gene[i] = paste(sort(unique(c(gene_i_nearest, gene_i_cis, gene_i_trans))), collapse = "|")
# }
#data_1$pos = as.numeric(data_1$pos)

# data_2 = data.frame(snp = locus$rsID[dum_2],
#                     chr = locus$chromosome[dum_2],
#                     pos = locus$position[dum_2],
#                     locus_ll = locus$region_start[dum_2],
#                     locus_ul = locus$region_end[dum_2],
#                     region = locus$region[dum_2],
#                     trait = "eGFR_female")
# for (i in 1 : nrow(data_2)) {
#   
#   dum_i = toplist_2$SNP == data_2$snp[i]
#   
#   #data_2$pos[i] = toplist_2$"position of SNP"[dum_i]
#   
#   gene_i_nearest = toplist_2$"nearest genes (functional relevance)"[dum_i]
#   gene_i_nearest = unlist(strsplit(unlist(strsplit(gene_i_nearest, split = "; ")), split = " \\("))
#   gene_i_nearest = gene_i_nearest[seq(1, length(gene_i_nearest), by = 2)]
#   
#   gene_i_cis = toplist_2$"cis-eQTL genes"[dum_i]
#   gene_i_cis = unlist(strsplit(gene_i_cis, split = " [|] "))
#   
#   gene_i_trans = toplist_2$"trans-eQTL genes"[dum_i]
#   gene_i_trans = unlist(strsplit(gene_i_trans, split = " [|] "))
#   
#   data_2$gene[i] = paste(sort(unique(c(gene_i_nearest, gene_i_cis, gene_i_trans))), collapse = "|")
# }
#data_2$pos = as.numeric(data_2$pos)





# ##### Lift-Over: Here not necessary but present to use general code!
# 
# ### Load step 10 data, perform liftover (Hg19 -> Hg38), and filter for particular variables
# #step10_1_files = grep(pattern = paste0("^step10_2_primaryANDsecondary_", data_1$trait[1]), step10_1_files, value = TRUE)
# load(paste0(path_load_step10_1, step10_1_files))
# step10_1       = erg1
# rm(erg1)
# 
# step10_1_lifted               = step10_1
# step10_1_lifted$chr_pos_oa_ea = paste(step10_1_lifted$chr, step10_1_lifted$pos, step10_1_lifted$other_allele, step10_1_lifted$effect_allele, sep = ":")
# 
# data_1_lifted = data_1
# 
# for (i in 1 : nrow(data_1)) {
#   
#   for (j in 1 : nrow(data_1)) {
#     
#     if(j != i) {
#       
#       if (((data_1_lifted$locus_ll[j] <= data_1_lifted$locus_ll[i]) & (data_1_lifted$locus_ul[j] >= data_1_lifted$locus_ll[i])) |
#           ((data_1_lifted$locus_ll[j] <= data_1_lifted$locus_ul[i]) & (data_1_lifted$locus_ul[j] >= data_1_lifted$locus_ul[i])) |
#           ((data_1_lifted$locus_ll[j] >= data_1_lifted$locus_ll[i]) & (data_1_lifted$locus_ul[j] <= data_1_lifted$locus_ul[i])) |
#           ((data_1_lifted$locus_ll[j] <= data_1_lifted$locus_ll[i]) & (data_1_lifted$locus_ul[j] >= data_1_lifted$locus_ul[i]))) {
#         
#         message(paste0("Intersection of Region ", data_1_lifted$region[i], " and Region ", data_1_lifted$region[j], "."))
#       }
#     }
#   }
# }
# 
# 
# 
# #step10_2_files = grep(pattern = paste0("^step10_2_primaryANDsecondary_", data_2$trait[i]), step10_2_files, value = TRUE)
# load(paste0(path_load_step10_2, step10_2_files))
# step10_2       = erg1
# rm(erg1)
# 
# step10_2_lifted               = step10_2
# step10_2_lifted$chr_pos_oa_ea = paste(step10_2_lifted$chr, step10_2_lifted$pos, step10_2_lifted$other_allele, step10_2_lifted$effect_allele, sep = ":")
# 
# data_2_lifted = data_2
# 
# for (i in 1 : nrow(data_2)) {
#   
#   for (j in 1 : nrow(data_2)) {
#     
#     if(j != i) {
#       
#       if (((data_2_lifted$locus_ll[j] <= data_2_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] >= data_2_lifted$locus_ll[i])) |
#           ((data_2_lifted$locus_ll[j] <= data_2_lifted$locus_ul[i]) & (data_2_lifted$locus_ul[j] >= data_2_lifted$locus_ul[i])) |
#           ((data_2_lifted$locus_ll[j] >= data_2_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] <= data_2_lifted$locus_ul[i])) |
#           ((data_2_lifted$locus_ll[j] <= data_2_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] >= data_2_lifted$locus_ul[i]))) {
#         
#         message(paste0("Intersection of Region ", data_2_lifted$region[i], " and Region ", data_2_lifted$region[j], "."))
#       }
#     }
#   }
# }





##### Colocalisation Analysis
# registerDoParallel(cores = min(detectCores(), max_cores))
# loop = foreach (j = 1 : length(eqtl_files)) %dopar% {



# plot(0, 0,
#      xlim = c((min(loci$region_start) - 1000000), (max(loci$region_stop) + 1000000)),
#      ylim = c(0, (nrow(loci) + 1)))
# for (i in nrow(loci)) {
#   
#   lines(c(loci$region_start[i], loci$region_stop[i]), c(i, i))
# }

# index_1 = c()
# index_2 = c()
# for (i in 1 : nrow(data_1)) {
#   
#   for (j in 1 : nrow(data_2)) {
#     
#     if (((data_2_lifted$locus_ll[j] <= data_1_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] >= data_1_lifted$locus_ll[i])) |
#         ((data_2_lifted$locus_ll[j] <= data_1_lifted$locus_ul[i]) & (data_2_lifted$locus_ul[j] >= data_1_lifted$locus_ul[i])) |
#         ((data_2_lifted$locus_ll[j] >= data_1_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] <= data_1_lifted$locus_ul[i])) |
#         ((data_2_lifted$locus_ll[j] <= data_1_lifted$locus_ll[i]) & (data_2_lifted$locus_ul[j] >= data_1_lifted$locus_ul[i]))) {
#       
#       message(paste0("Intersection of Region ", data_1_lifted$region[i], " in eGFR data and Region ", data_2_lifted$region[j], " in Uric Acid data."))
#       
#       index_1 = c(index_1, i)
#       index_2 = c(index_2, j)
#     }
#   }
# }
# 
# 
# 
# dum_cojo_cond_1 = c("",
#                     "",
#                     "",
#                     "",
#                     "",
#                     "",
#                     "",
#                     "")
# 
# 
# 
# dum_cojo_cond_2 = c("rs34884874",
#                     "rs34884874",
#                     "rs72616819",
#                     "rs112708523",
#                     "rs202138804",
#                     "rs7056552",
#                     "rs782581283",
#                     "rs4328011")
# 
# 
# 
# data_1_lifted = data_1_lifted[index_1, ]
# data_2_lifted = data_2_lifted[index_2, ]
data_1_lifted = data
data_2_lifted = data
data_4_lifted = data
data_5_lifted = data
coloc         = list()
for (i in 1 : 14) { # nrow(data_1_lifted) == nrow(data_2_lifted) == nrow(data_4_lifted) == nrow(data_5_lifted)
  
  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data ->
  step10_1_lifted = sumStat_union[sumStat_union$phenotype == "eGFR_MALE"]
  
  filt_1_locus_i = (step10_1_lifted$chromosome == data_1_lifted$chr[i]) &
    (step10_1_lifted$position >= data_1_lifted$locus_ll[i]) &
    (step10_1_lifted$position <= data_1_lifted$locus_ul[i]) &
    (step10_1_lifted$invalid_assoc == FALSE) &
    (step10_1_lifted$infoScore >= 0.8) &
    (step10_1_lifted$MAF >= 0.02) &
    (step10_1_lifted$numberOfStudies >= 10) &
    (step10_1_lifted$I2 <= 0.85)

  step10_1_lifted_filt_1_locus_i = step10_1_lifted[filt_1_locus_i, ]

  gwas_1_i = list(pvalues = step10_1_lifted_filt_1_locus_i$P,
                  N = as.numeric(step10_1_lifted_filt_1_locus_i$N),
                  MAF = step10_1_lifted_filt_1_locus_i$MAF,
                  beta = step10_1_lifted_filt_1_locus_i$beta,
                  varbeta = (step10_1_lifted_filt_1_locus_i$SE)^2,
                  type = rep("quant", (sum(filt_1_locus_i, na.rm = TRUE) + sum(is.na(filt_1_locus_i)))),
                  snp = step10_1_lifted_filt_1_locus_i$rsID)
  ### <-
  
  # ### Step/Option 2: Use cojo-slct Data ->
  # dum_m = ukbb_1$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_M_20210210[match(step10_1_lifted_filt_1_locus_i$id_meta, ukbb_1$markerID)]
  # dum_f = ukbb_1$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_F_20210210[match(step10_1_lifted_filt_1_locus_i$id_meta, ukbb_1$markerID)]
  # dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  # for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
  # 
  #   if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
  # 
  #     dum[l] = dum_f[l]
  #   }
  #   else {
  # 
  #     if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
  # 
  #       dum[l] = dum_m[l]
  #     }
  #     else {
  # 
  #       if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
  # 
  #         if (dum_m[l] == dum_f[l]) { # Case 2a
  # 
  #           dum[l] = dum_m[l]
  #         }
  #         else { # Case 2b
  # 
  #           message(paste0("Discrepancy in first data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
  #         }
  #       }
  #     }
  #   }
  # }
  # step10_1_lifted_filt_1_locus_i$id_ukbb = dum
  # 
  # # grep_i = grep(pattern = paste0("Region_", data_1_lifted$region[i], ".cma.cojo"),
  # #               dir(paste0(path_load_cojo_1, "results_select_V4/")),
  # #               value = TRUE)
  # # 
  # # cojo_slct_i = read.table(paste0(path_load_cojo_1, "results_select_V4/", grep_i),
  # #                          header = TRUE,
  # #                          sep = "\t")
  # # 
  # # matched                                  = match(step10_1_lifted_filt_1_locus_i$id_ukbb, cojo_slct_i$SNP)
  # # step10_1_lifted_filt_1_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # # step10_1_lifted_filt_1_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # # step10_1_lifted_filt_1_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # # gwas_1_i[["pvalues"]]                    = step10_1_lifted_filt_1_locus_i$p_slct
  # # gwas_1_i[["beta"]]                       = step10_1_lifted_filt_1_locus_i$beta_slct
  # # gwas_1_i[["varbeta"]]                    = (step10_1_lifted_filt_1_locus_i$se_slct)^2
  # ### <-
  # 
  # ### Step/Option 3: Use cojo-cond Data ->
  # grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_1_lifted$snp[i], ":"))[1], ".cma.cojo"),
  #               dir(paste0(path_load_cojo_1, "results_cond_V4/")),
  #               value = TRUE)
  # #grep_i = grep(pattern = paste0("Cond_results_", dum_cojo_cond_1[i], ".cma.cojo"),
  # #              dir(paste0(path_load_cojo_1, "results_cond_V4/")),
  # #              value = TRUE)
  # 
  # if (length(grep_i) > 0) {
  #   
  #   cojo_cond_i = read.table(paste0(path_load_cojo_1, "results_cond_V4/", grep_i),
  #                            header = TRUE,
  #                            sep = "\t")
  #   
  #   matched                                  = match(step10_1_lifted_filt_1_locus_i$id_ukbb, cojo_cond_i$SNP)
  #   step10_1_lifted_filt_1_locus_i$beta_cond = cojo_cond_i$bC[matched]
  #   step10_1_lifted_filt_1_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
  #   step10_1_lifted_filt_1_locus_i$p_cond    = cojo_cond_i$pC[matched]
  #   gwas_1_i[["pvalues"]]                    = step10_1_lifted_filt_1_locus_i$p_cond
  #   gwas_1_i[["beta"]]                       = step10_1_lifted_filt_1_locus_i$beta_cond
  #   gwas_1_i[["varbeta"]]                    = (step10_1_lifted_filt_1_locus_i$se_cond)^2
  # }
  # ### <-
  
  
  
  ### Step/Option 1: Use step10 Data ->
  step10_2_lifted = sumStat_union[sumStat_union$phenotype == "eGFR_FEMALE"]
  
  filt_2_locus_i = (step10_2_lifted$chromosome == data_2_lifted$chr[i]) &
    (step10_2_lifted$position >= data_2_lifted$locus_ll[i]) &
    (step10_2_lifted$position <= data_2_lifted$locus_ul[i]) &
    (step10_2_lifted$invalid_assoc == FALSE) &
    (step10_2_lifted$infoScore >= 0.8) &
    (step10_2_lifted$MAF >= 0.02) &
    (step10_2_lifted$numberOfStudies >= 10) &
    (step10_2_lifted$I2 <= 0.85)
  
  step10_2_lifted_filt_2_locus_i = step10_2_lifted[filt_2_locus_i, ]
  
  gwas_2_i = list(pvalues = step10_2_lifted_filt_2_locus_i$P,
                  N = as.numeric(step10_2_lifted_filt_2_locus_i$N),
                  MAF = step10_2_lifted_filt_2_locus_i$MAF,
                  beta = step10_2_lifted_filt_2_locus_i$beta,
                  varbeta = (step10_2_lifted_filt_2_locus_i$SE)^2,
                  type = rep("quant", (sum(filt_2_locus_i, na.rm = TRUE) + sum(is.na(filt_2_locus_i)))),
                  snp = step10_2_lifted_filt_2_locus_i$rsID)
  ### <-
  
  # ### Step/Option 2: Use cojo-slct Data ->
  # dum_m = ukbb_2$UKBB_EA_uric_acid_overall_HRC_UK10K_chrX_M_20210210[match(step10_2_lifted_filt_2_locus_i$id_meta, ukbb_2$markerID)]
  # dum_f = ukbb_2$UKBB_EA_uric_acid_overall_HRC_UK10K_chrX_F_20210210[match(step10_2_lifted_filt_2_locus_i$id_meta, ukbb_2$markerID)]
  # dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  # for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
  # 
  #   if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
  # 
  #     dum[l] = dum_f[l]
  #   }
  #   else {
  # 
  #     if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
  # 
  #       dum[l] = dum_m[l]
  #     }
  #     else {
  # 
  #       if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
  # 
  #         if (dum_m[l] == dum_f[l]) { # Case 2a
  # 
  #           dum[l] = dum_m[l]
  #         }
  #         else { # Case 2b
  # 
  #           message(paste0("Discrepancy in second data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
  #         }
  #       }
  #     }
  #   }
  # }
  # step10_2_lifted_filt_2_locus_i$id_ukbb = dum
  # 
  # # grep_i = grep(pattern = paste0("Region_", data_2_lifted$region[i], ".cma.cojo"),
  # #               dir(paste0(path_load_cojo_2, "results_select_V4/")),
  # #               value = TRUE)
  # # 
  # # cojo_slct_i = read.table(paste0(path_load_cojo_2, "results_select_V4/", grep_i),
  # #                          header = TRUE,
  # #                          sep = "\t")
  # # 
  # # matched                                  = match(step10_2_lifted_filt_2_locus_i$id_ukbb, cojo_slct_i$SNP)
  # # step10_2_lifted_filt_2_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # # step10_2_lifted_filt_2_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # # step10_2_lifted_filt_2_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # # gwas_2_i[["pvalues"]]                    = step10_2_lifted_filt_2_locus_i$p_slct
  # # gwas_2_i[["beta"]]                       = step10_2_lifted_filt_2_locus_i$beta_slct
  # # gwas_2_i[["varbeta"]]                    = (step10_2_lifted_filt_2_locus_i$se_slct)^2
  # ### <-
  # 
  # ### Step/Option 3: Use cojo-cond Data ->
  # #grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_2_lifted$snp[i], ":"))[1], ".cma.cojo"),
  # #              dir(paste0(path_load_cojo_2, "results_cond_V4/")),
  # #              value = TRUE)
  # grep_i = grep(pattern = paste0("Cond_results_", dum_cojo_cond_2[i], ".cma.cojo"),
  #               dir(paste0(path_load_cojo_2, "results_cond_V4/")),
  #               value = TRUE)
  # 
  # if (length(grep_i) > 0) {
  #   
  #   cojo_cond_i = read.table(paste0(path_load_cojo_2, "results_cond_V4/", grep_i),
  #                            header = TRUE,
  #                            sep = "\t")
  #   
  #   matched                                  = match(step10_2_lifted_filt_2_locus_i$id_ukbb, cojo_cond_i$SNP)
  #   step10_2_lifted_filt_2_locus_i$beta_cond = cojo_cond_i$bC[matched]
  #   step10_2_lifted_filt_2_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
  #   step10_2_lifted_filt_2_locus_i$p_cond    = cojo_cond_i$pC[matched]
  #   gwas_2_i[["pvalues"]]                    = step10_2_lifted_filt_2_locus_i$p_cond
  #   gwas_2_i[["beta"]]                       = step10_2_lifted_filt_2_locus_i$beta_cond
  #   gwas_2_i[["varbeta"]]                    = (step10_2_lifted_filt_2_locus_i$se_cond)^2
  # }
  # ### <-
  
  
  
  ### Perform Colocalisation Analyses
  # gene_1_i = unlist(strsplit(as.character(data_1_lifted$gene[i]), split = "[|]"))
  # gene_2_i = unlist(strsplit(as.character(data_2_lifted$gene[i]), split = "[|]"))
  
  coloc_ij = list()
  #for (k1 in 1 : length(gene_1_i)) {
    
    #for (k2 in 1 : length(gene_2_i)) {
      
      #if (length(eqtl_ij[[gene_i[k]]][[1]]) > 0) {
        
        dum_k_1 = is.na(gwas_1_i$pvalues) | is.na(gwas_1_i$N) | is.na(gwas_1_i$MAF) | is.na(gwas_1_i$beta) | is.na(gwas_1_i$varbeta) | is.na(gwas_1_i$type) | is.na(gwas_1_i$snp)
        dum_k_2 = is.na(gwas_2_i$pvalues) | is.na(gwas_2_i$N) | is.na(gwas_2_i$MAF) | is.na(gwas_2_i$beta) | is.na(gwas_2_i$varbeta) | is.na(gwas_2_i$type) | is.na(gwas_2_i$snp)
        
        gwas_1_filtered = gwas_1_i
        gwas_2_filtered = gwas_2_i
        for (l in 1 : length(gwas_1_filtered)) { # length(gwas_1_filtered) == length(gwas_2_filtered)
          
          gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!dum_k_1]
          gwas_2_filtered[[l]] = gwas_2_filtered[[l]][!dum_k_2]
        }
        
        dum_k = match(gwas_1_filtered[["snp"]], gwas_2_filtered[["snp"]])
        for (l in 1 : length(gwas_1_filtered)) { # length(gwas_1_filtered) == length(gwas_2_filtered)
          
          gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!is.na(dum_k)]
          gwas_2_filtered[[l]] = gwas_2_filtered[[l]][dum_k[!is.na(dum_k)]]
        }
        
        coloc_ij = coloc.abf(gwas_1_filtered, gwas_2_filtered)
      #}
    #}
  
  
  
  ### Save Colocalisation Analysis Results
  #save(coloc_ij, file = paste0(path_save_coloc, "_coloc", "_", gsub(":", "_", data_1_lifted$snp[i]), "_vs", "_", gsub(":", "_", data_2_lifted$snp[i]), ".RData"))
  #rm(coloc_ij)
  coloc[[i]] = coloc_ij
}
for (i in 15 : 20) { # nrow(data_1_lifted) == nrow(data_2_lifted) == nrow(data_4_lifted) == nrow(data_5_lifted)
  
  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data ->
  step10_1_lifted = sumStat_union[sumStat_union$phenotype == "UA_MALE"]
  
  filt_1_locus_i = (step10_1_lifted$chromosome == data_1_lifted$chr[i]) &
    (step10_1_lifted$position >= data_1_lifted$locus_ll[i]) &
    (step10_1_lifted$position <= data_1_lifted$locus_ul[i]) &
    (step10_1_lifted$invalid_assoc == FALSE) &
    (step10_1_lifted$infoScore >= 0.8) &
    (step10_1_lifted$MAF >= 0.02) &
    (step10_1_lifted$numberOfStudies >= 10) &
    (step10_1_lifted$I2 <= 0.85)
  
  step10_1_lifted_filt_1_locus_i = step10_1_lifted[filt_1_locus_i, ]
  
  gwas_1_i = list(pvalues = step10_1_lifted_filt_1_locus_i$P,
                  N = as.numeric(step10_1_lifted_filt_1_locus_i$N),
                  MAF = step10_1_lifted_filt_1_locus_i$MAF,
                  beta = step10_1_lifted_filt_1_locus_i$beta,
                  varbeta = (step10_1_lifted_filt_1_locus_i$SE)^2,
                  type = rep("quant", (sum(filt_1_locus_i, na.rm = TRUE) + sum(is.na(filt_1_locus_i)))),
                  snp = step10_1_lifted_filt_1_locus_i$rsID)
  ### <-
  
  # ### Step/Option 2: Use cojo-slct Data ->
  # dum_m = ukbb_1$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_M_20210210[match(step10_1_lifted_filt_1_locus_i$id_meta, ukbb_1$markerID)]
  # dum_f = ukbb_1$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_F_20210210[match(step10_1_lifted_filt_1_locus_i$id_meta, ukbb_1$markerID)]
  # dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  # for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
  # 
  #   if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
  # 
  #     dum[l] = dum_f[l]
  #   }
  #   else {
  # 
  #     if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
  # 
  #       dum[l] = dum_m[l]
  #     }
  #     else {
  # 
  #       if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
  # 
  #         if (dum_m[l] == dum_f[l]) { # Case 2a
  # 
  #           dum[l] = dum_m[l]
  #         }
  #         else { # Case 2b
  # 
  #           message(paste0("Discrepancy in first data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
  #         }
  #       }
  #     }
  #   }
  # }
  # step10_1_lifted_filt_1_locus_i$id_ukbb = dum
  # 
  # # grep_i = grep(pattern = paste0("Region_", data_1_lifted$region[i], ".cma.cojo"),
  # #               dir(paste0(path_load_cojo_1, "results_select_V4/")),
  # #               value = TRUE)
  # # 
  # # cojo_slct_i = read.table(paste0(path_load_cojo_1, "results_select_V4/", grep_i),
  # #                          header = TRUE,
  # #                          sep = "\t")
  # # 
  # # matched                                  = match(step10_1_lifted_filt_1_locus_i$id_ukbb, cojo_slct_i$SNP)
  # # step10_1_lifted_filt_1_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # # step10_1_lifted_filt_1_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # # step10_1_lifted_filt_1_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # # gwas_1_i[["pvalues"]]                    = step10_1_lifted_filt_1_locus_i$p_slct
  # # gwas_1_i[["beta"]]                       = step10_1_lifted_filt_1_locus_i$beta_slct
  # # gwas_1_i[["varbeta"]]                    = (step10_1_lifted_filt_1_locus_i$se_slct)^2
  # ### <-
  # 
  # ### Step/Option 3: Use cojo-cond Data ->
  # grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_1_lifted$snp[i], ":"))[1], ".cma.cojo"),
  #               dir(paste0(path_load_cojo_1, "results_cond_V4/")),
  #               value = TRUE)
  # #grep_i = grep(pattern = paste0("Cond_results_", dum_cojo_cond_1[i], ".cma.cojo"),
  # #              dir(paste0(path_load_cojo_1, "results_cond_V4/")),
  # #              value = TRUE)
  # 
  # if (length(grep_i) > 0) {
  #   
  #   cojo_cond_i = read.table(paste0(path_load_cojo_1, "results_cond_V4/", grep_i),
  #                            header = TRUE,
  #                            sep = "\t")
  #   
  #   matched                                  = match(step10_1_lifted_filt_1_locus_i$id_ukbb, cojo_cond_i$SNP)
  #   step10_1_lifted_filt_1_locus_i$beta_cond = cojo_cond_i$bC[matched]
  #   step10_1_lifted_filt_1_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
  #   step10_1_lifted_filt_1_locus_i$p_cond    = cojo_cond_i$pC[matched]
  #   gwas_1_i[["pvalues"]]                    = step10_1_lifted_filt_1_locus_i$p_cond
  #   gwas_1_i[["beta"]]                       = step10_1_lifted_filt_1_locus_i$beta_cond
  #   gwas_1_i[["varbeta"]]                    = (step10_1_lifted_filt_1_locus_i$se_cond)^2
  # }
  # ### <-
  
  
  
  ### Step/Option 1: Use step10 Data ->
  step10_2_lifted = sumStat_union[sumStat_union$phenotype == "UA_FEMALE"]
  
  filt_2_locus_i = (step10_2_lifted$chromosome == data_2_lifted$chr[i]) &
    (step10_2_lifted$position >= data_2_lifted$locus_ll[i]) &
    (step10_2_lifted$position <= data_2_lifted$locus_ul[i]) &
    (step10_2_lifted$invalid_assoc == FALSE) &
    (step10_2_lifted$infoScore >= 0.8) &
    (step10_2_lifted$MAF >= 0.02) &
    (step10_2_lifted$numberOfStudies >= 10) &
    (step10_2_lifted$I2 <= 0.85)
  
  step10_2_lifted_filt_2_locus_i = step10_2_lifted[filt_2_locus_i, ]
  
  gwas_2_i = list(pvalues = step10_2_lifted_filt_2_locus_i$P,
                  N = as.numeric(step10_2_lifted_filt_2_locus_i$N),
                  MAF = step10_2_lifted_filt_2_locus_i$MAF,
                  beta = step10_2_lifted_filt_2_locus_i$beta,
                  varbeta = (step10_2_lifted_filt_2_locus_i$SE)^2,
                  type = rep("quant", (sum(filt_2_locus_i, na.rm = TRUE) + sum(is.na(filt_2_locus_i)))),
                  snp = step10_2_lifted_filt_2_locus_i$rsID)
  ### <-
  
  # ### Step/Option 2: Use cojo-slct Data ->
  # dum_m = ukbb_2$UKBB_EA_uric_acid_overall_HRC_UK10K_chrX_M_20210210[match(step10_2_lifted_filt_2_locus_i$id_meta, ukbb_2$markerID)]
  # dum_f = ukbb_2$UKBB_EA_uric_acid_overall_HRC_UK10K_chrX_F_20210210[match(step10_2_lifted_filt_2_locus_i$id_meta, ukbb_2$markerID)]
  # dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  # for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
  # 
  #   if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
  # 
  #     dum[l] = dum_f[l]
  #   }
  #   else {
  # 
  #     if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
  # 
  #       dum[l] = dum_m[l]
  #     }
  #     else {
  # 
  #       if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
  # 
  #         if (dum_m[l] == dum_f[l]) { # Case 2a
  # 
  #           dum[l] = dum_m[l]
  #         }
  #         else { # Case 2b
  # 
  #           message(paste0("Discrepancy in second data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
  #         }
  #       }
  #     }
  #   }
  # }
  # step10_2_lifted_filt_2_locus_i$id_ukbb = dum
  # 
  # # grep_i = grep(pattern = paste0("Region_", data_2_lifted$region[i], ".cma.cojo"),
  # #               dir(paste0(path_load_cojo_2, "results_select_V4/")),
  # #               value = TRUE)
  # # 
  # # cojo_slct_i = read.table(paste0(path_load_cojo_2, "results_select_V4/", grep_i),
  # #                          header = TRUE,
  # #                          sep = "\t")
  # # 
  # # matched                                  = match(step10_2_lifted_filt_2_locus_i$id_ukbb, cojo_slct_i$SNP)
  # # step10_2_lifted_filt_2_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # # step10_2_lifted_filt_2_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # # step10_2_lifted_filt_2_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # # gwas_2_i[["pvalues"]]                    = step10_2_lifted_filt_2_locus_i$p_slct
  # # gwas_2_i[["beta"]]                       = step10_2_lifted_filt_2_locus_i$beta_slct
  # # gwas_2_i[["varbeta"]]                    = (step10_2_lifted_filt_2_locus_i$se_slct)^2
  # ### <-
  # 
  # ### Step/Option 3: Use cojo-cond Data ->
  # #grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_2_lifted$snp[i], ":"))[1], ".cma.cojo"),
  # #              dir(paste0(path_load_cojo_2, "results_cond_V4/")),
  # #              value = TRUE)
  # grep_i = grep(pattern = paste0("Cond_results_", dum_cojo_cond_2[i], ".cma.cojo"),
  #               dir(paste0(path_load_cojo_2, "results_cond_V4/")),
  #               value = TRUE)
  # 
  # if (length(grep_i) > 0) {
  #   
  #   cojo_cond_i = read.table(paste0(path_load_cojo_2, "results_cond_V4/", grep_i),
  #                            header = TRUE,
  #                            sep = "\t")
  #   
  #   matched                                  = match(step10_2_lifted_filt_2_locus_i$id_ukbb, cojo_cond_i$SNP)
  #   step10_2_lifted_filt_2_locus_i$beta_cond = cojo_cond_i$bC[matched]
  #   step10_2_lifted_filt_2_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
  #   step10_2_lifted_filt_2_locus_i$p_cond    = cojo_cond_i$pC[matched]
  #   gwas_2_i[["pvalues"]]                    = step10_2_lifted_filt_2_locus_i$p_cond
  #   gwas_2_i[["beta"]]                       = step10_2_lifted_filt_2_locus_i$beta_cond
  #   gwas_2_i[["varbeta"]]                    = (step10_2_lifted_filt_2_locus_i$se_cond)^2
  # }
  # ### <-
  
  
  
  ### Perform Colocalisation Analyses
  # gene_1_i = unlist(strsplit(as.character(data_1_lifted$gene[i]), split = "[|]"))
  # gene_2_i = unlist(strsplit(as.character(data_2_lifted$gene[i]), split = "[|]"))
  
  coloc_ij = list()
  #for (k1 in 1 : length(gene_1_i)) {
  
  #for (k2 in 1 : length(gene_2_i)) {
  
  #if (length(eqtl_ij[[gene_i[k]]][[1]]) > 0) {
  
  dum_k_1 = is.na(gwas_1_i$pvalues) | is.na(gwas_1_i$N) | is.na(gwas_1_i$MAF) | is.na(gwas_1_i$beta) | is.na(gwas_1_i$varbeta) | is.na(gwas_1_i$type) | is.na(gwas_1_i$snp)
  dum_k_2 = is.na(gwas_2_i$pvalues) | is.na(gwas_2_i$N) | is.na(gwas_2_i$MAF) | is.na(gwas_2_i$beta) | is.na(gwas_2_i$varbeta) | is.na(gwas_2_i$type) | is.na(gwas_2_i$snp)
  
  gwas_1_filtered = gwas_1_i
  gwas_2_filtered = gwas_2_i
  for (l in 1 : length(gwas_1_filtered)) { # length(gwas_1_filtered) == length(gwas_2_filtered)
    
    gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!dum_k_1]
    gwas_2_filtered[[l]] = gwas_2_filtered[[l]][!dum_k_2]
  }
  
  dum_k = match(gwas_1_filtered[["snp"]], gwas_2_filtered[["snp"]])
  for (l in 1 : length(gwas_1_filtered)) { # length(gwas_1_filtered) == length(gwas_2_filtered)
    
    gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!is.na(dum_k)]
    gwas_2_filtered[[l]] = gwas_2_filtered[[l]][dum_k[!is.na(dum_k)]]
  }
  
  coloc_ij = coloc.abf(gwas_1_filtered, gwas_2_filtered)
  #}
  #}
  
  
  
  ### Save Colocalisation Analysis Results
  #save(coloc_ij, file = paste0(path_save_coloc, "_coloc", "_", gsub(":", "_", data_1_lifted$snp[i]), "_vs", "_", gsub(":", "_", data_2_lifted$snp[i]), ".RData"))
  #rm(coloc_ij)
  coloc[[i]] = coloc_ij
}




# ##### Plot Results (Correlation Plot)
# res_h3    = data.frame(region = paste0("Region ", data_1_lifted$region, " (", data_1_lifted$snp, ")", " vs ", "Region ", data_2_lifted$region, " (", data_2_lifted$snp, ")"),
#                        pp = NA)
# res_h4    = data.frame(region = paste0("Region ", data_1_lifted$region, " (", data_1_lifted$snp, ")", " vs ", "Region ", data_2_lifted$region, " (", data_2_lifted$snp, ")"),
#                        pp = NA)
# res_h3_h4 = data.frame(region = paste0("Region ", data_1_lifted$region, " (", data_1_lifted$snp, ")", " vs ", "Region ", data_2_lifted$region, " (", data_2_lifted$snp, ")"),
#                        pp = NA)
# 
# res_files          = dir(paste0(path_save_coloc))
res_summary        = data.frame()
# for (i in 1 : nrow(data_1_lifted)) { # nrow(data_1_lifted) == nrow(data_2_lifted)
#   
#   if (sum(grepl(pattern = paste0("_coloc", "_", gsub(":", "_", data_1_lifted$snp[i]), "_vs", "_", gsub(":", "_", data_2_lifted$snp[i]), ".RData"),
#                 res_files)) > 0) {
#     
#     load(paste0(path_save_coloc, "_coloc", "_", gsub(":", "_", data_1_lifted$snp[i]), "_vs", "_", gsub(":", "_", data_2_lifted$snp[i]), ".RData"))
#     res_ij = coloc_ij
#     rm(coloc_ij)
#     
#     if (length(res_ij) > 0) {
#         
#       res_h3[i, 2] = res_ij[["summary"]][5]
#       res_h4[i, 2] = res_ij[["summary"]][6]
#       
#       res_summary = rbind(res_summary, res_ij[["summary"]])
#     }
#   }
# }
for (i in 1 : length(coloc)) {
  
  res_summary = rbind(res_summary, coloc[[i]]$summary)
}
names(res_summary) = c("nsnps",
                       "PP.H0.abf",
                       "PP.H1.abf",
                       "PP.H2.abf",
                       "PP.H3.abf",
                       "PP.H4.abf")
res_summary        = data.frame(locus = paste0("Region ", data$region),
                                trait1 = c(rep("eGFR (MALE)", 14), rep("UA (MALE)", 6)),
                                trait2 = c(rep("eGFR (FEMALE)", 14), rep("UA (FEMALE)", 6)),
                                res_summary)

# for (i in 1 : nrow(res_h3_h4)) {
#   
#   if (!is.na(res_h3[i, 2]) & !is.na(res_h4[i, 2])) {
#     
#     res_h3_h4[i, 2] = ifelse((res_h4[i, 2] >= res_h3[i, 2]), res_h4[i, 2], (-res_h3[i, 2]))
#   }
# }

# png(paste0(path_save_coloc, "_corrplot_20220906.png"),
#     width = 3000,
#     height = 3000,
#     res = 300)
# 
# res_matrix             = as.matrix(res_h3_h4[, -1])
# res_matrix_t           = t(res_matrix)
# rownames(res_matrix_t) = "eGFR (ALL) vs Uric Acid (ALL)"
# colnames(res_matrix_t) = paste0("Region ", data_1_lifted$region, " (", data_1_lifted$snp, ")", "\nvs\n", "Region ", data_2_lifted$region, " (", data_2_lifted$snp, ")")
# ggcp = ggcorrplot(res_matrix_t, title = "", lab = TRUE) +
#   scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = "grey", limits = c(-1, 1), name = "PP", breaks = seq(-1, 1, by = 0.25)) +
#   #scale_color_discrete() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# #res_melted = melt(res_h3_h4, id.vars = names(res_h3_h4)[1], measure.vars = names(res_h3_h4)[-1])
# # ggcp = ggplot(res_melted, aes(x = variable, y = gene, fill = value)) +
# #   geom_tile() +
# #   scale_colour_discrete() +
# #   theme(axis.title.x = element_blank(),
# #         axis.title.y = element_blank()) +
# #   ggtitle("")
# 
# print(ggcp)
# 
# dev.off()

write.table(res_summary,
            file = "../results/03_coloc_sex_ia.txt",
            sep = "\t",
            row.names = FALSE)





########## END OF FILE ##########





#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")




