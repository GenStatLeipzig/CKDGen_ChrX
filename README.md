# CKDGen Chr X Analyses

Analyses of X-chromosomal SNPs and kidney traits

**Last Updated: 29/09/2023**

Supporting code for the following draft:

* **Working title**: X-chromosome and kidney function: Evidence from a multi-trait genetic analysis of 908,697 individuals reveals sex-specific and sex-differential findings in genes regulated by androgen-response elements
* **Short title**: CKDGen ChrX
* **Analysis Team / Repository Contributor**: Markus Scholz, Katrin Horn, Andreas Kühnapfel and Janne Pott
* **Writing Team**: Markus Scholz, Katrin Horn, Afshin Parsa, Anna Köttgen, Pascal Schlosser, Cristian Pattaro

We performed a **trans-ethnic X-chromosome-wide association study of 7 kidney traits** (estimated glomerular filtration rate (eGFR), serum uric acid (UA), urine albumin-creatinine ratio (UACR), blood urea nitrogen (BUN), chronic kidney disease (CKD), gout and microalbuminuria (MA)). **Sex-stratified and combined analyses** were performed in up to 40 studies including up to 908,697 subjects considering up to 1,032,701 SNPs. Genome-wide significant loci were tested for sex-interactions and were compared between traits. A number of secondary analyses were performed to allocate candidate genes to the discovered loci. 

For more information, please contact Markus Scholz (markus.scholz@imise.uni-leipzig.de)

# Source File

If you want to reproduce our results, you will need to customize a source file, indicating

* path to R library (please use [R Version 4.x](https://cran.r-project.org/), all necessary packages are listed in the source file)
* path to data (summary statistics will be available on zenodo)
* path to [PLINK2](https://www.cog-genomics.org/plink/2.0/)
* path to [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview)
* path to [GTEx v8 data](https://gtexportal.org/home/protectedDataAccess)
* path to [NEPTUNE](https://nephqtl.org/)
* path to [HUNT](https://www.ntnu.edu/hunt) data
* path to [UKBB](https://www.ukbiobank.ac.uk/) data

# Statistical Analyses included

1) **Locus definition**: index SNP +/- 500 kb; merging of overlapping loci (union of loci)  
2) **Interaction analysis**: difference test of the stratified associations, correcting for beta-beta correlation between males and females  
3) **Co-localization analyses of male and female signals**: test of the stratfied associations for shared causal signal 
4) **Relevance to kidney function and CKD**: Comparison of eGFR and UA index SNP effects with those in CKD, BUN, and gout (one-sided test)
5) **Identification of independent variants**: *GCTA* - conditional joint (COJO) analyses
    * Reference data set: *UKBB* (best-guessed genotypes)
    * COJO select: step-wise forward selection to indentify independent variants
    * COJO conditional: Estimation of association statistics conditional to previously selected variants (in case of multiple independent variants at a locus) 
    * LD estimation: LD between independent variants per locus using *PLINK2*
6) **Credible set analyses**: Calculation of Approximate Bayes Factors using (conditional) effect estimates and standard errors
7) **Co-localization analyses with eQTLs**: test of the associations for shared causal signal with gene-expression of nearby genes using *GTEx v8* and *NEPTUNE* data
    * Extraction of eQTL data of nearby genes in 51 tissues 
    * Colocalization analyses per phenotype, tissue and gene
    * Checks and plots
8) **Analysis of overlap of eGFR and UA signals**:
    * LD estimation: LD between eGFR and UA index SNPs using *PLINK2*
    * Co-localization analyses of eGFR and UA:  test of the associations of eGFR and UA for shared causal signal 
9) **Replication analysis**: replicate findings of eGFR in *HUNT* study (one-sided)
10) **X-chromosomal heritability**: *GCTA* - restricted maximum likelihood (REML) analyses
    * Reference data set: *UKBB* (genotyped SNPs)
    * REML: estimation of heritability for eGFR and UA
    * REML bivariate: estimation of genetic correlation between eGFR and UA 
11) **Lookup of reported variants**: replicate findings of *Graham et al*, *Kanai et al* and *Sakaue et al* (successful if nominal significance and concordant effect is observed) 
12) **Trans-ethnic meta-regression analysis**: accounting for mixed ethnicities by using MR-Mega
    * Get an overview of MR-Mega results and check results for MetaGWAS loci
    * Generate Forest plots for MR-Mega hits   

# Statistical Analyses not included 

1) **Study quality control and harmonization**: pipeline from GenStatLeipzig, not yet available on GitHub
2) **Meta-analyses**: pipeline from GenStatLeipzig, not yet available on GitHub
3) **Bioinformatic annotation**: pipeline from GenStatLeipzig, not yet available on GitHub 
4) **Lookup of reported variants**: manual evaluation of GWAS catalog data of recent publications from *Graham et al*, *Kanai et al* and *Sakaue et al*
5) **Trans-ethnic meta-regression analysis**: done by Alexander Teumers lab (Greifswald, Germany), not yet shared
6) **Assignment of candidate genes**: manual evaluation of all previous applied methods, no code available


# Figures

## Main Figures

1) **Miami-Plot** of variants associated with eGFR and UA:
    * color coding: grey, overall; blue, males; red, females; black, not genome-wide significant
    * novelty coding: candidate gene in box novel, candidate gene without box kown
    * sex interaction coding: bold italic gene names indicate loci with sex interactions
    (--> see script F1_MiamiPlot.R)
2) **Regional Association Plot of locus 7 at Xq21.1**:
    * panel with all, male, female
3) **Beta-Beta-Plot of SNP-by-Sex interaction analysis**:
    * panel with results for eGFR and UA
    (--> see script 02_sex_ia.R)
4) **Heatmap** of cross-phenotype comparison of eGFR and UA loci:
    * six phenotypes: eGFR, UA, BUN, CKD, UACR and MA
    * 23 SNPs: index SNP per locus + female-specific index SNP of region 7
    (--> see script F2_Heatmap_phenotypes.R)
5) **Results of eQTL-colocalization analysis**:
    * six phenotypes: eGFR ALL, MALE, FEMALE and UA ALL, MALE, FEMALE
    * eight genes: CDKL5, CDK16, USP11, ARMCX2, TCEAL3, MORF4L2, ACSL4, SLC25A5
    (--> see scripts starting with 07_)


## Supplemental Figures

1) **RA-Plots** of all loci for ALL, MALE, FEMALE 
2) **Bar-Plot** of X-chromosomal heritability of eGFR and UA (--> see script 10_b_Heritability.R)
3) **Beta-Beta-Plot** of validation of eGFR hits in HUNT study (--> see script 09_replication_HUNT.R)
4) **Forest Plots** for SNP rs4328011 (--> see script 12_b_MR-Mega_make_forest_plots.R)
5) **Forest Plots** of MR-MEGA findings (--> see script 12_b_MR-Mega_make_forest_plots.R)
6) **Credible sets and missense mutations** (--> see script F4_MainFigure4_PostProbByCredSetSize.R)
7) **Study design** (not done with a script)
8) **Beta-Beta-Plot** of index SNPs comparing trans-ethnic metaGWAS and Europeans only metaGWAS
9) **QQ-Plots** of phenotypes and subgroups


# Tables 

## Main Tables

1) **Results of the transethnic X chromosome-wide association analysis of eGFR and UA** (--> see script T1_MainTable1.R)  
2) **Analysis of the overlapping of eGFR and UA loci** (--> see script T2_MainTable2.R)

## Supplemental Tables

--> see script ST_SupplementTables.R

1) Description of participating studies: study design and phenotype distribution (as received from participating studies)
2) Genotyping and imputation information of participating studies (as received from participating studies)
3) **Number of data sets, samples and SNPs contributing to the different association analyses** (--> see script ST_SupplementTables.R)
4) **Comparisons between sexes**, interaction and co-localization (--> see scripts 02_sex_ia.R and 03_coloc_sex_ia.R)
5) **Cross-phenotype comparision** (--> see script 04_lookup_TopHits_otherTraits.R)
6) **Independent variants per locus and analysis group** (--> see script 05_b_Cojo_Select.R and 05_c_Cojo_Conditional.R)
7) **Annotation of 99% credible sets** (done with pipeline from GenStatLeipzig)
8) **Co-localization of genetic association signals and eQTLs** (--> see scripts starting with 07_)
9) **Validation of eGFR associations in HUNT study** (--> see script 09_replication_HUNT.R)
10) **Look-up of SNPs previously reported for UA, eGFR, creatinine and BUN** (reported variants of three GWAS of kidney traits: *Graham et al*, *Kanai et al* and *Sakaue et al*) (--> see script 11_Lookup_Candidates_SNPs.R)
11) **Results of meta-regression analysis of the 23 index SNPs** (--> see scripts 12_a_MR-Mega_results_check.R)
12) **Additional genome-wide significant associations due to meta-regression** (--> see script 12_a_MR-Meta_results_check.R)
13) **Look-up of sex-biased gene-expressions of candidate genes assigned to genetic sex-interactions** (--> see script 13_Lookup_sexbiasedGE_eQTLs_PAR_ARE.R)
14) **Co-localization of genetic association signals of CKD traits and testosterone** (--> see script 07_d_coloc_testo.R)
15) **Look-up of eQTLs in case of positive co-localization** (--> see script 13_Lookup_eQTLs_leadSNPs.R)
