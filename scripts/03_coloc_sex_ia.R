#' ---
#' title: "Colocalisation Analysis: eGFR & Uric Acid: Male versus Female"
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
#' I want to check if these loci co-localise between the sexes (best phenotype was mainly ALL, but coloc does not require genome-wide significant results).
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
locus         = fread("../results/01_Locus_Definitions.txt")
sumStat_1     = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
sumStat_2     = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
sumStat_3     = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
sumStat_4     = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
sumStat_5     = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
sumStat_6     = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
sumStat_union = rbind(sumStat_1, sumStat_2, sumStat_3,
                      sumStat_4, sumStat_5, sumStat_6)

data = data.frame(snp = locus$rsID,
                  chr = locus$chromosome,
                  pos = locus$position,
                  locus_ll = locus$region_start,
                  locus_ul = locus$region_end,
                  region = locus$region,
                  trait = locus$phenotype)

data_1_lifted = data
data_2_lifted = data
data_4_lifted = data
data_5_lifted = data
coloc         = list()
for (i in 1 : 15) { # nrow(data_1_lifted) == nrow(data_2_lifted) == nrow(data_4_lifted) == nrow(data_5_lifted)
  
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
for (i in 16 : 22) { # nrow(data_1_lifted) == nrow(data_2_lifted) == nrow(data_4_lifted) == nrow(data_5_lifted)
  
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
res_summary        = data.frame(locus = paste0("Region ", data$region),
                                trait1 = c(rep("eGFR (MALE)", 15), rep("Uric Acid (MALE)", 7)),
                                trait2 = c(rep("eGFR (FEMALE)", 15), rep("Uric Acid (FEMALE)", 7)),
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




