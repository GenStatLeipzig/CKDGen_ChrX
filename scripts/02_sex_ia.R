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

source("../SourceFile_aman.R")

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

locus = locus[,c(1,5,6,8)]

dummy = data.table(region = 7,
                   rsID = "rs149995096:100479327:C:T",
                   phenotype = "eGFR_FEMALE",
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
                meta_male_union$P,
                meta_female_union$beta,
                meta_female_union$SE,
                meta_female_union$P,
                stringsAsFactors = FALSE)

df$mean_diff = rep(NA, nrow(df))
df$se_diff = rep(NA, nrow(df))
df$p_diff = rep(NA, nrow(df))

names(df) = c("SNP","Trait",
              "Male_Beta", "Male_SE", "Male_P",
              "Female_Beta", "Female_SE", "Female_P",
              "mean_diff","se_diff","p_diff")

union_r        = c(meta_male$rsID, meta_female$rsID)
union_unique_r = unique(union_r)
r              = cor(meta_male$beta[match(union_unique_r, meta_male$rsID)],
                     meta_female$beta[match(union_unique_r, meta_female$rsID)],
                     use = "complete.obs",
                     method = "spearman")
for (i in 1 : nrow(df)) {

  z_i = (df$"Male_Beta"[i] - df$"Female_Beta"[i]) / sqrt(df$"Male_SE"[i]^2 + df$"Female_SE"[i]^2 - 2 * r * df$"Male_SE"[i] * df$"Female_SE"[i])
  p_i = 2 * pnorm(abs(z_i), lower.tail = FALSE)
  diff_i = (df$"Male_Beta"[i] - df$"Female_Beta"[i])
  diff_se_i = sqrt(df$"Male_SE"[i]^2 + df$"Female_SE"[i]^2 - 2 * r * df$"Male_SE"[i] * df$"Female_SE"[i])
  
  df$mean_diff[i] = diff_i
  df$se_diff[i] = diff_se_i
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

bp = ggplot(data = df, aes(x = df$"Male_Beta", y = df$"Female_Beta")) +
       lims(x = c(ifelse((s == 1), -0.035, -0.05), ifelse((s == 1), 0.005, 0.05)), y = c(ifelse((s == 1), -0.035, -0.05), ifelse((s == 1), 0.005, 0.05))) +
       geom_hline(yintercept = 0, col = "grey") +
       geom_vline(xintercept = 0, col = "grey") +
       geom_point(size = 2, colour = col, alpha = 0.25) +
       geom_errorbar(aes(ymin = (df$"Female_Beta" - 1.96 * df$"Female_SE"), ymax = (df$"Female_Beta" + 1.96 * df$"Female_SE")), colour = col, alpha = 0.25) +
       geom_errorbarh(aes(xmin = (df$"Male_Beta" - 1.96 * df$"Male_SE"), xmax = (df$"Male_Beta" + 1.96 * df$"Male_SE")), colour = col, alpha = 0.25) +
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


#' # Supplemental Figure: Beta-Beta Plot ####
#' ***
sexIA_eGFR = fread("../results/02_sex_ia_egfr.txt")
sexIA_UA = fread("../results/02_sex_ia_uric_acid.txt")
filt = grepl("UA",locus$phenotype)
table(filt)

sexIA = rbind(sexIA_eGFR[!filt,],sexIA_UA[filt,])
setDT(sexIA)
sexIA[,fdr_diff := p.adjust(p_diff, method = "fdr")]

myPlotData<-data.table(rs_id=sexIA$SNP,
                       trait=sexIA$Trait,
                       beta_male = sexIA$Male_Beta,
                       beta_female = sexIA$Female_Beta,
                       se_male = sexIA$Male_SE,
                       se_female = sexIA$Female_SE,
                       meandiff = sexIA$mean_diff,
                       meandiff_p = sexIA$p_diff,
                       meandiff_p_FDR = sexIA$fdr_diff)
myPlotData[meandiff_p>=0.05,sig:="no"]
myPlotData[meandiff_p<0.05,sig:="yes"]
myPlotData[meandiff_p_FDR<0.05,sig:="yes (FDR 5%)"]
myPlotData[,gene2:=""]
myPlotData[meandiff_p_FDR<0.05,gene2:=gsub(":.*","",rs_id)]
myPlotData[,gene3:=""]
myPlotData[meandiff_p<0.05,gene3:=gsub(":.*","",rs_id)]
myPlotData[, sex.higherEffect := "none"]
myPlotData[(abs(beta_female) > abs(beta_male)) & sig != "no", sex.higherEffect := "female"]
myPlotData[(abs(beta_male) > abs(beta_female)) & sig != "no", sex.higherEffect := "male"]


myPlot1 = ggplot(myPlotData, aes(x=beta_male, y=beta_female, color=sex.higherEffect)) + 
  facet_wrap(~trait, scales = "free") +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = beta_female- 1.96*se_female, ymax = beta_female+ 1.96*se_female)) +
  geom_errorbarh(aes(xmin = beta_male- 1.96*se_male, xmax = beta_male + 1.96*se_male)) +
  theme_bw(base_size = 10) +
  scale_colour_manual(values=c("#B2182B","#2166AC","#000000"),
                      labels=c("female","male", "none")) +
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20)) +
  labs(x="Effect Size Male", 
       y = "Effect Size Female",
       color = "SNPs with \nsex interaction") +
  geom_label_repel(data = subset(myPlotData, sig!="no" & beta_male>-0.02),
                   aes(x=beta_male, y=beta_female, label = gene3),
                   xlim = c(-0.02,-0.005),nudge_y = 0.001)+
  geom_label_repel(data = subset(myPlotData, sig!="no" & beta_male< -0.02),
                   aes(x=beta_male, y=beta_female, label = gene3),
                   ylim = c(0,0.02),nudge_x = 0.01)+
  guides(label="none", color="none")
myPlot1

tiff(filename = "../figures/SupplementalFigure_BetaBeta_sexIA.tiff", 
     width = 3000, height = 1800, res=300, compression = 'lzw')
myPlot1
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")




