#####################################################
########## PERFORM COLOCALISATION ANALYSIS ##########
#####################################################





#rm(list = ls())
rm(list = setdiff(ls(), c("coloc_scripts", "coloc_script")))





##### Setup
computer  = "angmar"
#max_cores = 2





##### Load Packages
#.libPaths(paste0("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/", computer))
library("readxl")
library("biomaRt")
library("GenomicRanges")
library("rtracklayer")
#library("doParallel")
library("data.table")
#?getDTthreads
getDTthreads()
setDTthreads(1)
getDTthreads()
library("coloc")
#library("reshape2")
library("ggplot2")
library("ggcorrplot")





##### Directories
path_load_step10  = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/04_Annotation/gwasresults_1.7_nstud10_V2/"
path_load_toplist = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/04_Annotation/gwasresults_1.7_nstud10_V2/synopsis/"
path_load_cojo    = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/07_cojo/"
path_load_loci    = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/_results_combined/"
path_load_ukbb    = "J:/genstat/02_projekte/1911_CKDgen_chrX/09_meta_FileQC/"
path_load_gtex    = "J:/genstat/01_daten/2007_GTEx_v8/only_chromosome_23/"
path_load_nephqtl = "J:/genstat/01_daten/2203_nephqtl_matrixeqtl_scrnaseq/only_chromosome_23/"
path_save_coloc   = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/"





##### Files
step10_files  = grep(pattern = "^step10_2_primaryANDsecondary_", dir(path_load_step10), value = TRUE)
eqtl_files    = grep(pattern = ".gz_chrX.RData$", dir(path_load_gtex), value = TRUE)
eqtl_names    = unlist(strsplit(eqtl_files, split = "^GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_"))[seq(from = 2, to = (2 * length(eqtl_files)), by = 2)]
eqtl_names    = unlist(strsplit(eqtl_names, split = ".allpairs.txt.gz_chrX.RData"))
neptune_files = grep(pattern = ".gz_chrX.RData$", dir(path_load_nephqtl), value = TRUE)
neptune_files
neptune_names = c("Glomerulus", "Tubulointerstitial")
neptune_names





##### Individual Data
toplist = read_excel(paste0(path_load_toplist, "step80_annotated_toplist_validated_2022-04-26_ALL_eGFR_overall_nstud10_V2.xlsx"), sheet = "topliste")
toplist = as.data.frame(toplist)

cojo = read_excel(paste0(path_load_loci, "Von_Wuttke_to_Lokus_220530.xlsx"), sheet = "Independent_SNPs")
cojo = as.data.frame(cojo)

loci = read_excel(paste0(path_load_loci, "CKDGen_Region_SNPs_in_ALL_F_M_summary_220608_interaction_20220608.xlsx"), sheet = "Overview")
loci = as.data.frame(loci)

ukbb = fread(paste0(path_load_ukbb, "FileQC_ALL_eGFR/step25_ID_Matching_Table.txt.gz"))

dum_filtered = (cojo$phenotype == "eGFR") & (cojo$subanalysis == "ALL")
region_list  = sort(unique(cojo$Region[dum_filtered]))

data         = data.frame(snp = cojo$SNP[dum_filtered],
                          chr = 23,
                          pos = NA,
                          locus_ll = loci$region_start[match(cojo$Region[dum_filtered], loci$Region)],
                          locus_ul = loci$region_stop[match(cojo$Region[dum_filtered], loci$Region)],
                          region = cojo$Region[dum_filtered],
                          cyto = cojo$Cytoband[dum_filtered],
                          trait = "eGFR_overall",
                          gene = NA,
                          ensembl = NA,
                          entrez = NA)
for (i in 1 : nrow(data)) {

  dum_i = toplist$SNP == data$snp[i]

  data$pos[i] = toplist$"position of SNP"[dum_i]

  gene_i_nearest = toplist$"nearest genes (functional relevance)"[dum_i]
  gene_i_nearest = unlist(strsplit(unlist(strsplit(gene_i_nearest, split = "; ")), split = " \\("))
  gene_i_nearest = gene_i_nearest[seq(1, length(gene_i_nearest), by = 2)]

  gene_i_cis = toplist$"cis-eQTL genes"[dum_i]
  gene_i_cis = unlist(strsplit(gene_i_cis, split = " [|] "))

  gene_i_trans = toplist$"trans-eQTL genes"[dum_i]
  gene_i_trans = unlist(strsplit(gene_i_trans, split = " [|] "))

  data$gene[i] = paste(sort(unique(c(gene_i_nearest, gene_i_cis, gene_i_trans))), collapse = "|")
}
data$pos     = as.numeric(data$pos)

########## Use "R_biomaRt.R" locally!
genes = unlist(strsplit(data$gene, split = "[|]"))
write.table(genes,
            file = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_genes.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
# ...
# genes        = read.table("J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_genes.txt")
# genes        = genes[, 1]
# genes_unique = unique(genes)
# 
# 
# 
# #listMarts()
# mart     = useMart("ensembl")
# #listDatasets(mart)
# mart     = useDataset("hsapiens_gene_ensembl", mart = mart)
# #listAttributes(mart)
# #listFilters(mart)
# matching = getBM(attributes= c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "description", "gene_biotype"),
#                  filters = "external_gene_name",
#                  values = genes,
#                  mart = mart)
# 
# write.table(matching,
#             file = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_genes_ensembls.txt",
#             sep = "\t",
#             row.names = FALSE)
# ...
matching = read.table("J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_genes_ensembls.txt",
                      header = TRUE,
                      sep = "\t")
##########

for (i in 1 : nrow(data)) {
  
  dum_i = unlist(strsplit(data$gene[i], split = "[|]"))
  
  dum_ij = c()
  for (j in 1 : length(dum_i)) {
    
    dum_j = matching$ensembl_gene_id[matching$external_gene_name == dum_i[j]]
    if (length(dum_j) > 0) {
      
      dum_ij = c(dum_ij, dum_j)
    }
    else {
      
      dum_ij = c(dum_ij, NA)
    }
  }
  if (length(dum_ij) == 0) {
    
    data$ensembl[i] = NA
  }
  else {
    
    data$ensembl[i] = paste(dum_ij, collapse = "|")
  }
  
  dum_ij = c()
  for (j in 1 : length(dum_i)) {
    
    dum_j = matching$entrezgene_id[matching$external_gene_name == dum_i[j]]
    if (length(dum_j) > 0) {
      
      dum_ij = c(dum_ij, dum_j)
    }
    else {
      
      dum_ij = c(dum_ij, NA)
    }
  }
  if (length(dum_ij) == 0) {
    
    data$entrez[i] = NA
  }
  else {
    
    data$entrez[i] = paste(dum_ij, collapse = "|")
  }
}





##### Lift-Over

### Load step 10 data, perform liftover (Hg19 -> Hg38), and filter for particular variables
step10_files         = grep(pattern = paste0("^step10_2_primaryANDsecondary_", data$trait[i]), step10_files, value = TRUE)
load(paste0(path_load_step10, step10_files))
step10               = erg1
rm(erg1)
step10$chr_pos_oa_ea = paste(step10$chr, step10$pos, step10$other_allele, step10$effect_allele, sep = ":")

hg19           = step10[, c("markername", "chr", "pos")]
names(hg19)[1] = "snps"

########## Use "R_biomaRt.R" locally!
write.table(hg19,
            file = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_snps.txt",
            sep = "\t",
            row.names = FALSE)
# ...
# source("J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/00_scripts/liftOverHg19TOHg38_ak.R")
# snps          = read.table("J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_snps.txt",
#                            header = TRUE,
#                            sep = "\t")
# snps_liftover = liftOVerHg19TOHg38(snps)
# write.table(snps_liftover,
#             file = "J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_snps_liftover.txt",
#             sep = "\t",
#             row.names = FALSE)
# ...
hg38 = read.table("J:/genstat/02_projekte/1911_CKDgen_chrX/10_metaGWAS/01_eGFR_allEth_sex_combined/09_coloc/coloc_snps_liftover.txt",
                  header = TRUE,
                  sep = "\t")
##########

step10_lifted              = step10
match_hg19_hg38            = match(hg19$pos, hg38$pos)
dum_na                     = is.na(hg38$chr_hg38) | is.na(hg38$pos_hg38)
step10_lifted$chr[!dum_na] = hg38$chr_hg38[match_hg19_hg38][!dum_na]
step10_lifted$chr          = as.numeric(step10_lifted$chr)
step10_lifted$pos[!dum_na] = hg38$pos_hg38[match_hg19_hg38][!dum_na]

step10_lifted$chr_pos_oa_ea = paste(step10_lifted$chr, step10_lifted$pos, step10_lifted$other_allele, step10_lifted$effect_allele, sep = ":")

### Also for filtered toplist
data_lifted = data
for (i in 1 : length(data_lifted$pos)) {
  
  data_lifted$pos[i]      = hg38$pos_hg38[hg38$pos == data_lifted$pos[i]]
  data_lifted$locus_ll[i] = data_lifted$pos[i] - (data$pos[i] - data$locus_ll[i])
  data_lifted$locus_ul[i] = data_lifted$pos[i] + (data$locus_ul[i] - data$pos[i])
}





##### Colocalisation Analysis

### Read eQTL Data: GTEx
#registerDoParallel(cores = min(detectCores(), max_cores))
#loop_eqtl = foreach (j = 1 : length(eqtl_files)) %dopar% {
for (j in 1 : length(eqtl_files)) {

### Read eQTL Data: GTEx
load(paste0(path_load_gtex, eqtl_files[j]))
eqtl_j = eqtl_files_i_chrX
rm(eqtl_files_i_chrX)

for (i in 1 : nrow(data_lifted)) {
  
  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data ->
  filt_locus_i = (step10_lifted$chr == data_lifted$chr[i]) &
    (step10_lifted$pos >= min(data_lifted$locus_ll[i], data_lifted$locus_ll[i])) &
    (step10_lifted$pos <= max(data_lifted$locus_ul[i], data_lifted$locus_ul[i])) &
    (step10_lifted$invalid_assocs_topl == FALSE) &
    (step10_lifted$info >= 0.8) &
    (step10_lifted$maf >= 0.02)
  
  step10_lifted_filt_locus_i = step10_lifted[filt_locus_i, ]
  step10_filt_locus_i        = step10[filt_locus_i, ]
  
  gwas_i = list(pvalues = step10_lifted_filt_locus_i$p_score,
                N = as.numeric(step10_lifted_filt_locus_i$n),
                MAF = step10_lifted_filt_locus_i$maf,
                beta = step10_lifted_filt_locus_i$beta_score,
                varbeta = (step10_lifted_filt_locus_i$se_score)^2,
                type = rep("quant", nrow(step10_filt_locus_i)),
                snp = step10_lifted_filt_locus_i$chr_pos_oa_ea)
  ### <-
  
  ### Step/Option 2: Use cojo-slct Data ->
  dum_m = ukbb$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_M_20210210[match(step10_filt_locus_i$id_meta, ukbb$markerID)]
  dum_f = ukbb$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_F_20210210[match(step10_filt_locus_i$id_meta, ukbb$markerID)]
  dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
    
    if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
      
      dum[l] = dum_f[l]
    }
    else {
      
      if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
        
        dum[l] = dum_m[l]
      }
      else {
        
        if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
          
          if (dum_m[l] == dum_f[l]) { # Case 2a
            
            dum[l] = dum_m[l]
          }
          else { # Case 2b
            
            message(paste0("Discrepancy in first data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
          }
        }
      }
    }
  }
  step10_lifted_filt_locus_i$id_ukbb = dum
  
  # grep_i = grep(pattern = paste0("Region_", data_lifted$region[i], ".cma.cojo"),
  #               dir(paste0(path_load_cojo, "results_select_V4/")),
  #               value = TRUE)
  # 
  # cojo_slct_i = read.table(paste0(path_load_cojo, "results_select_V4/", grep_i),
  #                          header = TRUE,
  #                          sep = "\t")
  # 
  # matched                              = match(step10_lifted_filt_locus_i$id_ukbb, cojo_slct_i$SNP)
  # step10_lifted_filt_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # step10_lifted_filt_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # step10_lifted_filt_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # gwas_i[["pvalues"]]                  = step10_lifted_filt_locus_i$p_slct
  # gwas_i[["beta"]]                     = step10_lifted_filt_locus_i$beta_slct
  # gwas_i[["varbeta"]]                  = (step10_lifted_filt_locus_i$se_slct)^2
  ### <-
  
  ### Step/Option 3: Use cojo-cond Data ->
  grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_lifted$snp[i], ":"))[1], ".cma.cojo"),
                dir(paste0(path_load_cojo, "results_cond_V4/")),
                value = TRUE)
  
  if (length(grep_i) > 0) {
    
    cojo_cond_i = read.table(paste0(path_load_cojo, "results_cond_V4/", grep_i),
                             header = TRUE,
                             sep = "\t")
    
    matched                              = match(step10_lifted_filt_locus_i$id_ukbb, cojo_cond_i$SNP)
    step10_lifted_filt_locus_i$beta_cond = cojo_cond_i$bC[matched]
    step10_lifted_filt_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
    step10_lifted_filt_locus_i$p_cond    = cojo_cond_i$pC[matched]
    gwas_i[["pvalues"]]                  = step10_lifted_filt_locus_i$p_cond
    gwas_i[["beta"]]                     = step10_lifted_filt_locus_i$beta_cond
    gwas_i[["varbeta"]]                  = (step10_lifted_filt_locus_i$se_cond)^2
  }
  ### <-
  
  
  
  ### Filter eQTL Data for Candidate Genes
  gene_i    = unlist(strsplit(as.character(data_lifted$gene[i]), split = "[|]"))
  ensembl_i = unlist(strsplit(as.character(data_lifted$ensembl[i]), split = "[|]"))
  
  filt_ij   = rep(FALSE, nrow(eqtl_j))
  for (k in 1 : length(ensembl_i)) {
    
    filt_ijk = grepl(pattern = ensembl_i[k], eqtl_j$gene_id)
    filt_ij  = filt_ij | filt_ijk
  }
  eqtl_j_filt_ij = eqtl_j[filt_ij, ]
  
  
  
  if (nrow(eqtl_j_filt_ij) > 0) {
    
    ### Create SNP ID for eQTL Data
    dum_j = unlist(strsplit(eqtl_j_filt_ij$variant_id, split = "_"))
    
    chr = dum_j[seq(1, length(dum_j), by = 5)]
    chr = unlist(strsplit(chr, split = "chr"))[seq(from = 2, to = (2 * length(chr)), by = 2)]
    for (l in 1 : length(chr)) {
      
      if (chr[l] == "X") {
        
        chr[l] = 23
      }
    }
    chr = as.numeric(chr)
    
    pos           = as.numeric(dum_j[seq(from = 2, to = length(dum_j), by = 5)])
    other_allele  = dum_j[seq(from = 3, to = length(dum_j), by = 5)]
    effect_allele = dum_j[seq(from = 4, to = length(dum_j), by = 5)]
    
    eqtl_j_filt_ij$chr           = chr
    eqtl_j_filt_ij$pos           = pos
    eqtl_j_filt_ij$other_allele  = other_allele
    eqtl_j_filt_ij$effect_allele = effect_allele
    eqtl_j_filt_ij$chr_pos_oa_ea = paste(chr, pos, other_allele, effect_allele, sep = ":")
    
    
    
    ### Create Locus-Wise Lists as Input for Colocalisation Analyses
    eqtl_ij = list()
    for (k in 1 : length(ensembl_i)) {
      
      eqtl_j_filt_ij_k = grepl(ensembl_i[k], eqtl_j_filt_ij$gene_id)
      
      eqtl_j_k         = list(pvalues = eqtl_j_filt_ij$pval_nominal[eqtl_j_filt_ij_k], 
                              N = eqtl_j_filt_ij$ma_count[eqtl_j_filt_ij_k], 
                              MAF = eqtl_j_filt_ij$maf[eqtl_j_filt_ij_k], 
                              beta = eqtl_j_filt_ij$slope[eqtl_j_filt_ij_k], 
                              varbeta = (eqtl_j_filt_ij$slope_se[eqtl_j_filt_ij_k])^2, 
                              type = rep("quant", (sum(eqtl_j_filt_ij_k, na.rm = TRUE) + sum(is.na(eqtl_j_filt_ij_k)))),
                              snp = eqtl_j_filt_ij$chr_pos_oa_ea[eqtl_j_filt_ij_k])
      
      eqtl_ij[[gene_i[k]]] = eqtl_j_k
    }
    
    
    
    ### Perform Colocalisation Analyses
    coloc_ij = list()
    for (k in 1 : length(ensembl_i)) {
      
      if (length(eqtl_ij[[gene_i[k]]][[1]]) > 0) {
        
        dum_k_eqtl = is.na(eqtl_ij[[gene_i[k]]]$pvalues) | is.na(eqtl_ij[[gene_i[k]]]$N) | is.na(eqtl_ij[[gene_i[k]]]$MAF) | is.na(eqtl_ij[[gene_i[k]]]$beta) | is.na(eqtl_ij[[gene_i[k]]]$varbeta) | is.na(eqtl_ij[[gene_i[k]]]$type) | is.na(eqtl_ij[[gene_i[k]]]$snp)
        dum_k_gwas = is.na(gwas_i$pvalues) | is.na(gwas_i$N) | is.na(gwas_i$MAF) | is.na(gwas_i$beta) | is.na(gwas_i$varbeta) | is.na(gwas_i$type) | is.na(gwas_i$snp)
        
        eqtl_filtered = eqtl_ij[[gene_i[k]]]
        gwas_filtered = gwas_i
        for (l in 1 : length(eqtl_filtered)) { # length(eqtl_filtered) == length(gwas_filtered)
          
          eqtl_filtered[[l]] = eqtl_filtered[[l]][!dum_k_eqtl]
          gwas_filtered[[l]] = gwas_filtered[[l]][!dum_k_gwas]
        }
        
        sum_valid = length(intersect(eqtl_filtered[["snp"]], gwas_filtered[["snp"]]))
        if (sum_valid > 0) {
          
          coloc_ij_abf_k = coloc.abf(eqtl_filtered, gwas_filtered)
          coloc_ij[[gene_i[k]]] = coloc_ij_abf_k
        }
      }
    }
    
    
    
    ### Save Colocalisation Analysis Results
    save(coloc_ij, file = paste0(path_save_coloc, "coloc", "_chr_", data_lifted$chr[i], "_pos_", data_lifted$pos[i], "_trait_", data_lifted$trait[i], "_eqtl_", eqtl_names[j], ".RData"))
    rm(coloc_ij)
  }
}

rm(eqtl_j)
}



### Read eQTL Data: NephQTL
#registerDoParallel(cores = min(detectCores(), max_cores))
#loop_neptune = foreach (j = 1 : length(neptune_files)) %dopar% {
for (j in 1 : length(neptune_files)) {

load(paste0(path_load_nephqtl, neptune_files[j]))
neptune_j        = neptune_files_i_chrX
rm(neptune_files_i_chrX)
names(neptune_j) = c("chr",
                     "pos",
                     "snp",
                     "ea",
                     "oa",
                     "entrez",
                     "beta",
                     "sebeta",
                     "t",
                     "pvalues")

for (i in 1 : nrow(data)) {

  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data ->
  filt_locus_i = (step10$chr == data$chr[i]) &
    (step10$pos >= min(data$locus_ll[i], data$locus_ll[i])) &
    (step10$pos <= max(data$locus_ul[i], data$locus_ul[i])) &
    (step10$invalid_assocs_topl == FALSE) &
    (step10$info >= 0.8) &
    (step10$maf >= 0.02)
  
  step10_filt_locus_i = step10[filt_locus_i, ]
  
  gwas_i = list(pvalues = step10_filt_locus_i$p_score,
                N = as.numeric(step10_filt_locus_i$n),
                MAF = step10_filt_locus_i$maf,
                beta = step10_filt_locus_i$beta_score,
                varbeta = (step10_filt_locus_i$se_score)^2,
                type = rep("quant", nrow(step10_filt_locus_i)),
                snp = step10_filt_locus_i$chr_pos_oa_ea)
  ### <-
  
  ### Step/Option 2: Use cojo-slct Data ->
  dum_m = ukbb$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_M_20210210[match(step10_filt_locus_i$id_meta, ukbb$markerID)]
  dum_f = ukbb$UKBB_EA_eGFR_overall_HRC_UK10K_chrX_F_20210210[match(step10_filt_locus_i$id_meta, ukbb$markerID)]
  dum   = rep(NA, length(dum_m)) # length(dum_m) == length(dum_f)
  for (l in 1 : length(dum_m)) { # length(dum_m) == length(dum_f)
    
    if (is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 1a
      
      dum[l] = dum_f[l]
    }
    else {
      
      if (!is.na(dum_m[l]) & is.na(dum_f[l])) { # Case 1b
        
        dum[l] = dum_m[l]
      }
      else {
        
        if (!is.na(dum_m[l]) & !is.na(dum_f[l])) { # Case 2
          
          if (dum_m[l] == dum_f[l]) { # Case 2a
            
            dum[l] = dum_m[l]
          }
          else { # Case 2b
            
            message(paste0("Discrepancy in first data in line ", l, ": dum_m = ", dum_m[l], ", dum_f = ", dum_f[l]))
          }
        }
      }
    }
  }
  step10_filt_locus_i$id_ukbb = dum
  
  # grep_i = grep(pattern = paste0("Region_", data_lifted$region[i], ".cma.cojo"),
  #               dir(paste0(path_load_cojo, "results_select_V4/")),
  #               value = TRUE)
  # 
  # cojo_slct_i = read.table(paste0(path_load_cojo, "results_select_V4/", grep_i),
  #                          header = TRUE,
  #                          sep = "\t")
  # 
  # matched                       = match(step10_filt_locus_i$id_ukbb, cojo_slct_i$SNP)
  # step10_filt_locus_i$beta_slct = cojo_slct_i$bC[matched]
  # step10_filt_locus_i$se_slct   = cojo_slct_i$bC_se[matched]
  # step10_filt_locus_i$p_slct    = cojo_slct_i$pC[matched]
  # gwas_i[["pvalues"]]           = step10_filt_locus_i$p_slct
  # gwas_i[["beta"]]              = step10_filt_locus_i$beta_slct
  # gwas_i[["varbeta"]]           = (step10_filt_locus_i$se_slct)^2
  ### <-
  
  ### Step/Option 3: Use cojo-cond Data ->
  grep_i = grep(pattern = paste0("Cond_results_", unlist(strsplit(data_lifted$snp[i], ":"))[1], ".cma.cojo"),
                dir(paste0(path_load_cojo, "results_cond_V4/")),
                value = TRUE)
  
  if (length(grep_i) > 0) {
    
    cojo_cond_i = read.table(paste0(path_load_cojo, "results_cond_V4/", grep_i),
                             header = TRUE,
                             sep = "\t")
    
    matched                       = match(step10_filt_locus_i$id_ukbb, cojo_cond_i$SNP)
    step10_filt_locus_i$beta_cond = cojo_cond_i$bC[matched]
    step10_filt_locus_i$se_cond   = cojo_cond_i$bC_se[matched]
    step10_filt_locus_i$p_cond    = cojo_cond_i$pC[matched]
    gwas_i[["pvalues"]]           = step10_filt_locus_i$p_cond
    gwas_i[["beta"]]              = step10_filt_locus_i$beta_cond
    gwas_i[["varbeta"]]           = (step10_filt_locus_i$se_cond)^2
  }
  ### <-
  
  
  
  ### Filter eQTL Data for Candidate Genes
  gene_i   = unlist(strsplit(as.character(data$gene[i]), split = "[|]"))
  entrez_i = unlist(strsplit(as.character(data$entrez[i]), split = "[|]"))

  filt_ij   = rep(FALSE, nrow(neptune_j))
  for (k in 1 : length(entrez_i)) {

    filt_ijk = grepl(pattern = entrez_i[k], neptune_j$entrez)

    filt_ij = filt_ij | filt_ijk
  }
  neptune_j_filt_ij = neptune_j[filt_ij, ]



  if (nrow(neptune_j_filt_ij) > 0) {

    ### Create SNP ID for eQTL Data
    chr = neptune_j_filt_ij$chr
    for (l in 1 : length(chr)) {

      if (chr[l] == "X") {

        chr[l] = 23
      }
    }
    chr = as.numeric(chr)

    pos           = as.numeric(neptune_j_filt_ij$pos)
    other_allele  = neptune_j_filt_ij$oa
    effect_allele = neptune_j_filt_ij$ea

    neptune_j_filt_ij$chr           = chr
    neptune_j_filt_ij$pos           = pos
    neptune_j_filt_ij$other_allele  = other_allele
    neptune_j_filt_ij$effect_allele = effect_allele
    neptune_j_filt_ij$chr_pos_oa_ea = paste(chr, pos, other_allele, effect_allele, sep = ":")



    ### Create Locus-Wise Lists as Input for Colocalisation Analyses
    neptune_ij = list()
    for (k in 1 : length(entrez_i)) {

      neptune_j_filt_ij_k = grepl(paste0("^", entrez_i[k], "$"), neptune_j_filt_ij$entrez) # Entrez is only a number -> exact matching needed!

      neptune_j_k         = list(pvalues = neptune_j_filt_ij$pvalues[neptune_j_filt_ij_k],
                                 N = rep(ifelse((neptune_names[j] == "Glomerulus"), 136, 166), (sum(neptune_j_filt_ij_k, na.rm = TRUE) + sum(is.na(neptune_j_filt_ij_k)))),
                                 MAF = step10$maf[match(neptune_j_filt_ij$pos[neptune_j_filt_ij_k], step10$pos)],
                                 beta = neptune_j_filt_ij$beta[neptune_j_filt_ij_k],
                                 varbeta = (neptune_j_filt_ij$sebeta[neptune_j_filt_ij_k])^2,
                                 type = rep("quant", (sum(neptune_j_filt_ij_k, na.rm = TRUE) + sum(is.na(neptune_j_filt_ij_k)))),
                                 snp = neptune_j_filt_ij$chr_pos_oa_ea[neptune_j_filt_ij_k])

      neptune_ij[[gene_i[k]]] = neptune_j_k
    }



    ### Perform Colocalisation Analyses
    coloc_ij = list()
    for (k in 1 : length(entrez_i)) {

      if (length(neptune_ij[[gene_i[k]]][[1]]) > 0) {

        dum_k_neptune = is.na(neptune_ij[[gene_i[k]]]$pvalues) | is.na(neptune_ij[[gene_i[k]]]$N) | is.na(neptune_ij[[gene_i[k]]]$MAF) | is.na(neptune_ij[[gene_i[k]]]$beta) | is.na(neptune_ij[[gene_i[k]]]$varbeta) | is.na(neptune_ij[[gene_i[k]]]$type) | is.na(neptune_ij[[gene_i[k]]]$snp)
        #dum_k_neptune = is.na(neptune_ij[[gene_i[k]]]$pvalues) | is.na(neptune_ij[[gene_i[k]]]$N) | is.na(neptune_ij[[gene_i[k]]]$beta) | is.na(neptune_ij[[gene_i[k]]]$varbeta) | is.na(neptune_ij[[gene_i[k]]]$type) | is.na(neptune_ij[[gene_i[k]]]$snp)
        dum_k_gwas    = is.na(gwas_i$pvalues) | is.na(gwas_i$N) | is.na(gwas_i$MAF) | is.na(gwas_i$beta) | is.na(gwas_i$varbeta) | is.na(gwas_i$type) | is.na(gwas_i$snp)
        #dum_k_gwas    = is.na(gwas_i$pvalues) | is.na(gwas_i$N) | is.na(gwas_i$beta) | is.na(gwas_i$varbeta) | is.na(gwas_i$type) | is.na(gwas_i$snp)

        neptune_filtered = neptune_ij[[gene_i[k]]]
        gwas_filtered    = gwas_i
        for (l in 1 : length(neptune_filtered)) { # length(neptune_filtered) == length(gwas_filtered)

          neptune_filtered[[l]] = neptune_filtered[[l]][!dum_k_neptune]
          gwas_filtered[[l]]    = gwas_filtered[[l]][!dum_k_gwas]
        }
        
        sum_valid = length(intersect(neptune_filtered[["snp"]], gwas_filtered[["snp"]]))
        if (sum_valid > 0) {
          
          coloc_ij_abf_k        = coloc.abf(neptune_filtered, gwas_filtered)
          coloc_ij[[gene_i[k]]] = coloc_ij_abf_k
        }
      }
    }



    ### Save Colocalisation Analysis Results
    save(coloc_ij, file = paste0(path_save_coloc, "coloc", "_chr_", data$chr[i], "_pos_", data$pos[i], "_trait_", data$trait[i], "_eqtl_", neptune_names[j], ".RData"))
    rm(coloc_ij)
  }
}

rm(neptune_j)
}



### Plot Results (Correlation Plot)
gene_names            = unlist(strsplit(data_lifted$gene, split = "[|]"))
snp_names             = rep(data_lifted$snp, unlist(lapply(strsplit(data_lifted$gene, split = "[|]"), length)))
region_names          = rep(data_lifted$region, unlist(lapply(strsplit(data_lifted$gene, split = "[|]"), length)))
region_names_ordering = order(region_names, gene_names, decreasing = TRUE)

res_h3    = data.frame(gene = gene_names,
                       snp = snp_names)
res_h4    = data.frame(gene = gene_names,
                       snp = snp_names)
res_h3_h4 = data.frame(gene = gene_names,
                       snp = snp_names)

res_files        = dir(paste0(path_save_coloc))
res_summary_eqtl = data.frame()
locus_eqtl       = c()
gene_eqtl        = c()
trait2_eqtl      = c()
for (j in 1 : length(eqtl_files)) {
  
  res_h3                            = cbind(res_h3, NA)
  names(res_h3)[ncol(res_h3)]       = paste0("pp_", eqtl_names[j])
  res_h4                            = cbind(res_h4, NA)
  names(res_h4)[ncol(res_h4)]       = paste0("pp_", eqtl_names[j])
  res_h3_h4                         = cbind(res_h3_h4, NA)
  names(res_h3_h4)[ncol(res_h3_h4)] = paste0("pp_", eqtl_names[j])
  for (i in 1 : nrow(data_lifted)) {
    
    if (sum(grepl(pattern = paste0("coloc", "_chr_", data_lifted$chr[i], "_pos_", data_lifted$pos[i], "_trait_", data_lifted$trait[i], "_eqtl_", eqtl_names[j], ".RData"),
                  res_files)) > 0) {
      
      load(paste0(path_save_coloc, "coloc", "_chr_", data_lifted$chr[i], "_pos_", data_lifted$pos[i], "_trait_", data_lifted$trait[i], "_eqtl_", eqtl_names[j], ".RData"))
      res_ij = coloc_ij
      rm(coloc_ij)
      
      if (length(res_ij) > 0) {
        
        for (k in 1 : length(res_ij)) {
          
          dum_k = (gene_names == names(res_ij)[k]) & (snp_names == data_lifted$snp[i])
          
          res_h3[dum_k, ncol(res_h3)] = res_ij[[k]][["summary"]][5]
          res_h4[dum_k, ncol(res_h4)] = res_ij[[k]][["summary"]][6]
          
          res_summary_eqtl = rbind(res_summary_eqtl, res_ij[[k]][["summary"]])
          locus_eqtl       = c(locus_eqtl, paste0("Region ", data_lifted$region[i]))
          gene_eqtl        = c(gene_eqtl, names(res_ij)[k])
          trait2_eqtl      = c(trait2_eqtl, paste0("GTEx", " (", gsub(pattern = "_", replacement = " ", eqtl_names[j]), ")"))
        }
      }
    }
  }
}
names(res_summary_eqtl) = c("nsnps",
                            "PP.H0.abf",
                            "PP.H1.abf",
                            "PP.H2.abf",
                            "PP.H3.abf",
                            "PP.H4.abf")
res_summary_eqtl        = data.frame(locus = locus_eqtl,
                                     gene = gene_eqtl,
                                     trait1 = "eGFR (ALL)",
                                     trait2 = trait2_eqtl,
                                     res_summary_eqtl)
res_summary_neptune     = data.frame()
locus_neptune           = c()
gene_neptune            = c()
trait2_neptune          = c()
for (j in 1 : length(neptune_files)) {
  
  res_h3                            = cbind(res_h3, NA)
  names(res_h3)[ncol(res_h3)]       = paste0("pp_", neptune_names[j])
  res_h4                            = cbind(res_h4, NA)
  names(res_h4)[ncol(res_h4)]       = paste0("pp_", neptune_names[j])
  res_h3_h4                         = cbind(res_h3_h4, NA)
  names(res_h3_h4)[ncol(res_h3_h4)] = paste0("pp_", neptune_names[j])
  for (i in 1 : nrow(data)) {
    
    if (sum(grepl(pattern = paste0("coloc", "_chr_", data$chr[i], "_pos_", data$pos[i], "_trait_", data$trait[i], "_eqtl_", neptune_names[j], ".RData"),
                  res_files)) > 0) {
      
      load(paste0(path_save_coloc, "coloc", "_chr_", data$chr[i], "_pos_", data$pos[i], "_trait_", data$trait[i], "_eqtl_", neptune_names[j], ".RData"))
      res_ij = coloc_ij
      rm(coloc_ij)
      
      if (length(res_ij) > 0) {
        
        for (k in 1 : length(res_ij)) {
          
          dum_k = (gene_names == names(res_ij)[k]) & (snp_names == data$snp[i])
          
          res_h3[dum_k, ncol(res_h3)] = res_ij[[k]][["summary"]][5]
          res_h4[dum_k, ncol(res_h4)] = res_ij[[k]][["summary"]][6]
          
          res_summary_neptune = rbind(res_summary_neptune, res_ij[[k]][["summary"]])
          locus_neptune       = c(locus_neptune, paste0("Region ", data$region[i]))
          gene_neptune        = c(gene_neptune, names(res_ij)[k])
          trait2_neptune      = c(trait2_neptune, paste0("NEPTUNE", " (", gsub(pattern = "_", replacement = " ", neptune_names[j]), ")"))
        }
      }
    }
  }
}
names(res_summary_neptune) = c("nsnps",
                               "PP.H0.abf",
                               "PP.H1.abf",
                               "PP.H2.abf",
                               "PP.H3.abf",
                               "PP.H4.abf")
res_summary_neptune        = data.frame(locus = locus_neptune,
                                        gene = gene_neptune,
                                        trait1 = "eGFR (ALL)",
                                        trait2 = trait2_neptune,
                                        res_summary_neptune)

for (j in 1 : (length(eqtl_names) + length(neptune_names))) {
  
  for (i in 1 : nrow(res_h3_h4)) {
    
    if (!is.na(res_h3[i, (2 + j)]) & !is.na(res_h4[i, (2 + j)])) {
      
      res_h3_h4[i, (2 + j)] = ifelse((res_h4[i, (2 + j)] >= res_h3[i, (2 + j)]), res_h4[i, (2 + j)], (-res_h3[i, (2 + j)]))
    }
  }
}

png(paste0(path_save_coloc, "corrplot_20220906.png"),
    width = 15000,
    height = 10000,
    res = 300)

res_matrix             = as.matrix(res_h3_h4[, -c(1, 2)])
res_matrix_t           = t(res_matrix)
rownames(res_matrix_t) = c(gsub(pattern = "_", replacement = " ", eqtl_names), neptune_names)
colnames(res_matrix_t) = paste0("Region ", region_names, " (", snp_names, " | ", gene_names, ")")
ggcp = ggcorrplot(res_matrix_t[, region_names_ordering], title = "eGFR (ALL)") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = "grey", limits = c(-1, 1), name = "PP", breaks = seq(-1, 1, by = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(ggcp)

dev.off()

res_summary = rbind(res_summary_eqtl, res_summary_neptune)
write.table(res_summary,
            file = paste0(path_save_coloc, "summary_20220906.txt"),
            sep = "\t",
            row.names = FALSE)





########## END OF FILE ##########




