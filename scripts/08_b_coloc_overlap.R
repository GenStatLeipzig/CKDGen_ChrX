#' ---
#' title: "Colocalisation Analysis: eGFR versus Uric Acid"
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
#' I want to check if common loci between both phenotypes co-localise.
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list = setdiff(ls(), "allinfo"))
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(projectpath_main)





#####################################################
########## PERFORM COLOCALISATION ANALYSIS ##########
#####################################################

##### Setup
#





##### Common Base Data
locus     = fread("../results/01_Locus_Definitions.txt")
sumStat_1 = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
sumStat_2 = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
sumStat_3 = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
sumStat_4 = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
sumStat_5 = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
sumStat_6 = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
sumStat   = list(sumStat_1, sumStat_2, sumStat_3,
                 sumStat_4, sumStat_5, sumStat_6)

data = data.frame(snp = locus$rsID,
                  chr = locus$chromosome,
                  pos = locus$position,
                  locus_ll = locus$region_start,
                  locus_ul = locus$region_end,
                  region = locus$region,
                  trait = locus$phenotype,
                  trait2 = c(rep(3, 3), 1, rep(3, 11),
                             rep(6, 7)))

index_1 = c()
index_2 = c()
for (i in 1 : 15) {
  
  for (j in 16 : 22) {
    
    if (((data$locus_ll[j] <= data$locus_ll[i]) & (data$locus_ul[j] >= data$locus_ll[i])) |
        ((data$locus_ll[j] <= data$locus_ul[i]) & (data$locus_ul[j] >= data$locus_ul[i])) |
        ((data$locus_ll[j] >= data$locus_ll[i]) & (data$locus_ul[j] <= data$locus_ul[i])) |
        ((data$locus_ll[j] <= data$locus_ll[i]) & (data$locus_ul[j] >= data$locus_ul[i]))) {
      
      message(paste0("Intersection of Region ", data$region[i], " in eGFR data and Region ", data$region[j], " in Uric Acid data."))
      
      index_1 = c(index_1, i)
      index_2 = c(index_2, j)
    }
  }
}

# Manual inclusion of female eGFR locus
index_1 = c(index_1[c(1 : 3)], 7, index_1[c(4 : 7)])
index_2 = c(index_2[c(1 : 3)], 18, index_2[c(4 : 7)])
#

data_1 = data[index_1, ]
data_2 = data[index_2, ]

# Manual inclusion of female eGFR locus
sumStat_2_egfr_female = sumStat_2[(sumStat_2$rsID == "rs149995096:100479327:C:T"), ]
levels(data_1$snp)    = c(levels(data_1$snp), sumStat_2_egfr_female$rsID)
data_1$snp[4]         = sumStat_2_egfr_female$rsID
data_1$pos[4]         = sumStat_2_egfr_female$position
levels(data_1$trait)  = c(levels(data_1$trait), "eGFR_FEMALE")
data_1$trait[4]       = "eGFR_FEMALE"
data_1$trait2[4]      = 2
#

coloc = list()
for (i in 1 : nrow(data_1)) { # nrow(data_1) == nrow(data_2)
  
  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data ->
  step10_1_lifted = sumStat[[data_1$trait2[i]]]
  
  filt_1_locus_i = (step10_1_lifted$chromosome == data_1$chr[i]) &
    (step10_1_lifted$position >= min(data_1$locus_ll[i], data_2$locus_ll[i])) &
    (step10_1_lifted$position <= max(data_1$locus_ul[i], data_2$locus_ul[i])) &
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
  
  # ### Step/Option 2: Use cojo-cond Data ->
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
  step10_2_lifted = sumStat[[data_2$trait2[i]]]
  
  filt_2_locus_i = (step10_2_lifted$chromosome == data_2$chr[i]) &
    (step10_2_lifted$position >= min(data_1$locus_ll[i], data_2$locus_ll[i])) &
    (step10_2_lifted$position <= max(data_1$locus_ul[i], data_2$locus_ul[i])) &
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
  
  # ### Step/Option 2: Use cojo-cond Data ->
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
  coloc_ij = list()
  
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
  
  
  
  ### Save Colocalisation Analysis Results
  coloc[[i]] = coloc_ij
}





##### Write Results
res_summary = data.frame()
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
res_summary        = data.frame(region1 = paste0("Region ", data_1$region),
                                region2 = paste0("Region ", data_2$region),
                                trait1 = c("eGFR (MALE)", rep("eGFR (ALL)", 2), "eGFR (FEMALE)", rep("eGFR (ALL)", 4)),
                                trait2 = rep("Uric Acid (ALL)", 8),
                                res_summary)

# for (i in 1 : nrow(res_h3_h4)) {
#   
#   if (!is.na(res_h3[i, 2]) & !is.na(res_h4[i, 2])) {
#     
#     res_h3_h4[i, 2] = ifelse((res_h4[i, 2] >= res_h3[i, 2]), res_h4[i, 2], (-res_h3[i, 2]))
#   }
# }

write.table(res_summary,
            file = "../results/08_b_coloc_overlap.txt",
            sep = "\t",
            row.names = FALSE)





########## END OF FILE ##########





#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")




