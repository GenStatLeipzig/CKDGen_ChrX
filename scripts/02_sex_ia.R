#' ---
#' title: "Sex Interaction"
#' subtitle: "CKDGen - Chr X"
#' author: "Andreas, Janne Pott"
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
#' I want to check for sex interaction.
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list = setdiff(ls(), "allinfo"))
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(projectpath_main)





#####################################################################################################
##### SEX-SPECIFIC ASSOCIATION ANALYSIS FOR CANDIDATE SNPS FROM CKDGEN GWAS META ANALYSIS: EGFR #####
#####################################################################################################

##### Setup
#





##### Common Base Data
locus                 = fread("../results/01_Locus_Definitions.txt")
meta_egfr_male        = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
meta_egfr_female      = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
meta_egfr_all         = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
meta_uric_acid_male   = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
meta_uric_acid_female = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
meta_uric_acid_all    = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")

locus = locus[,c(1,5,8)]
dummy = data.table(region = 7,
                   rsID = "rs149995096:100479327:C:T",
                   position = 100479327)
locus = rbind(locus,dummy)
setorder(locus,region,position)
locus


##### Analysis
for (s in 1 : 2) {

if (s == 1) {
  
  meta_male   = meta_egfr_male
  meta_female = meta_egfr_female
  meta_all    = meta_egfr_all
  
  other_meta_male   = meta_uric_acid_male
  other_meta_female = meta_uric_acid_female
  other_meta_all    = meta_uric_acid_all
}

if (s == 2) {
  
  meta_male   = meta_uric_acid_male
  meta_female = meta_uric_acid_female
  meta_all    = meta_uric_acid_all
  
  other_meta_male   = meta_egfr_male
  other_meta_female = meta_egfr_female
  other_meta_all    = meta_egfr_all
}

### Unite Results of Males and Females
meta_all_gwsig    = meta_all[(meta_all$P <= (5 * 10^(-8))), ]
meta_male_gwsig   = meta_male[(meta_male$P <= (5 * 10^(-8))), ]
meta_female_gwsig = meta_female[(meta_female$P <= (5 * 10^(-8))), ]

meta_all_gwsig_ext_other    = rbind(meta_all_gwsig, other_meta_all[(other_meta_all$P <= (5 * 10^(-8))), ])
meta_male_gwsig_ext_other   = rbind(meta_male_gwsig, other_meta_male[(other_meta_male$P <= (5 * 10^(-8))), ])
meta_female_gwsig_ext_other = rbind(meta_female_gwsig, other_meta_female[(other_meta_female$P <= (5 * 10^(-8))), ])

union        = c(meta_all_gwsig_ext_other$rsID, meta_male_gwsig_ext_other$rsID, meta_female_gwsig_ext_other$rsID)
union_unique = unique(union)

meta_male_union   = meta_male[match(union_unique, meta_male$rsID), ]
meta_female_union = meta_female[match(union_unique, meta_female$rsID), ]



### Create Table
df = data.frame(union_unique,
                ifelse((s == 1), "eGFR", "Uric Acid"),
                meta_male_union$beta,
                meta_male_union$SE,
                meta_female_union$beta,
                meta_female_union$SE,
                stringsAsFactors = FALSE)

df$p_diff = rep(NA, nrow(df))

names(df) = c("SNP",
              "Trait",
              "Male Beta",
              "Male SE",
              "Female Beta",
              "Female SE",
              "p_diff")

union_r        = c(meta_male$rsID, meta_female$rsID)
union_unique_r = unique(union_r)
r              = cor(meta_male$beta[match(union_unique_r, meta_male$rsID)],
                     meta_female$beta[match(union_unique_r, meta_female$rsID)],
                     use = "complete.obs",
                     method = "spearman")
for (i in 1 : nrow(df)) {

  z_i = (df$"Male Beta"[i] - df$"Female Beta"[i]) / sqrt(df$"Male SE"[i]^2 + df$"Female SE"[i]^2 - 2 * r * df$"Male SE"[i] * df$"Female SE"[i])
  p_i = 2 * pnorm(abs(z_i), lower.tail = FALSE)
  
  df$p_diff[i] = p_i
}

# Match Top-SNPs
df = df[match(locus$rsID, df$SNP), ]

df$fdr_diff = p.adjust(df$p_diff, method = "fdr")

col                      = rep("black", nrow(df))
col[(df$p_diff <= 0.05)] = "red"



### Plot Male versus Female
png(filename = paste0("../results/02_sex_ia_beta_plot_", ifelse((s == 1), "egfr", "uric_acid"), ".png"),
    height = 2000,
    width = 2000,
    res = 300)

bp = ggplot(data = df, aes(x = df$"Male Beta", y = df$"Female Beta")) +
       lims(x = c(ifelse((s == 1), -0.035, -0.05), ifelse((s == 1), 0.005, 0.05)), y = c(ifelse((s == 1), -0.035, -0.05), ifelse((s == 1), 0.005, 0.05))) +
       geom_hline(yintercept = 0, col = "grey") +
       geom_vline(xintercept = 0, col = "grey") +
       geom_point(size = 2, colour = col, alpha = 0.25) +
       geom_errorbar(aes(ymin = (df$"Female Beta" - 1.96 * df$"Female SE"), ymax = (df$"Female Beta" + 1.96 * df$"Female SE")), colour = col, alpha = 0.25) +
       geom_errorbarh(aes(xmin = (df$"Male Beta" - 1.96 * df$"Male SE"), xmax = (df$"Male Beta" + 1.96 * df$"Male SE")), colour = col, alpha = 0.25) +
       #geom_label_repel(aes(label = df$SNP)) +
       labs(title = ifelse((s == 1), "eGFR", "Uric Acid"), x = "Effect Size Male", y = "Effect Size Female") +
       #scale_color_manual(breaks = c("FDR <= 5%", "FDR > 5%"), values = c("red", "black")) +
       theme_bw()
print(bp)

dev.off()





### Save data
write.table(df,
            file = paste0("../results/02_sex_ia_", ifelse((s == 1), "egfr", "uric_acid"), ".txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
}





#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")




