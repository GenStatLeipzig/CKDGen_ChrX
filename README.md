# CKDGen Chr X Analyses

Analyses of X-chromosomal SNPs and kidney traits

**Last Updated: 07/11/2022**

Supporting code for the following draft:

* Working title: Chromosome X genetic association analysis of kidney traits in up to 908,697 subjects revealed six novel loci and sex-specific hits in genes regulated by androgen response elements
* Short title: CKDGen ChrX GWAMA
* Analysis Team / Repository Contributor: Markus Scholz, Katrin Horn, Andreas Kühnapfel and Janne Pott
* Writing Team: Markus Scholz, Afshin Parsa, Pascal Schlosser, Cristian Pattaro

We performed a **trans-ethnic X-chromosome-wide association study of 7 kidney traits** (estimated glomerular filtration rate (eGFR), serum uric acid (UA), urine albumin-creatinine ratio (UACR), blood urea nitrogen (BUN), chronic kidney disease (CKD), gout and microalbuminuria (MA)). **Sex-stratified and combined analyses** were performed in up to 46 different studies including up to 908,697 subjects considering up to 1,032,701 SNPs. Genome-wide significant loci were tested for sex-interactions and were compared between traits. A number of secondary analyses were performed to allocate candidate genes to the discovered loci. 

For more information, please contact Markus Scholz (markus.scholz@imise.uni-leipzig.de)

# Source File

If you want to reproduce our results, you will need to customize a source file, indicating

* path to R library (please use [R Version 4.x](https://cran.r-project.org/), all necessary packages are listed in the source file)
* path to data (summary statistics will be available on zenodo)
* path to [PLINK2](https://www.cog-genomics.org/plink/2.0/)
* path to [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview)
* path to [GTEx v8 data](https://gtexportal.org/home/protectedDataAccess)
* path to NEPTUNE eQTL data
* path to HUNT data
* path to UKBB data

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
7) **Co-localization analyses with eQTLs**: test of the (conditinal) associations for shared causal signal with gene-expression of nearby genes using *GTEx v8* and *NEPTUNE* data
8) **Analysis of overlap of eGFR and UA signals**:
    * LD estimation: LD between eGFR and UA index SNPs using *PLINK2*
    * Co-localization analyses of eGFR and UA:  test of the associations of eGFR and UA for shared causal signal 
9) **Replication analysis**: replicate findings of eGFR in *HUNT* study (one-sided)
10) X-chromosomal heritability: *GCTA* - restricted maximum likelihood (REML) analyses
    * Reference data set: *UKBB* (genotyped SNPs)
    * REML: estimation of heritability for eGFR and UA
    * REML bivariate: estimation of genetic correlation between eGFR and UA 
11) **Lookup of reported variants**: replicate findings of *Graham et al*, *Kanai et al* and *Sakaue et al* (successful if nominal significance and concordant effect is observed) 
12) Trans-ethnic meta-regression analysis: accounting for mixed ethnicities by using MR-Mega

# Statistical Analyses not included 

1) **Study quality control and harmonization**: pipeline from GenStatLeipzig, not yet available on github
2) **Meta-analyses**: pipeline from GenStatLeipzig, not yet available on github
3) **Bioinformatic annotation**: pipeline from GenStatLeipzig, not yet available on github 
4)  **Lookup of reported variants**: pipeline from GenStatLeipzig, not yet available on github 
5) **Trans-ethnic meta-regression analysis**: done by Alexander Teumers lab (Greifswald, Germany), not yet shared
6) **Assignment of candidate genes**: manual evaluation of all previous applied methods, no code available


# Figures (not yet completed)

## Main Figures

1) **Miami Plot** (eGFR vs UA):
    * color coding: blue, males; red, females
    * novelty coding: candidate gene in black novel, candidate gene in grey kown
2) Heatmap of p-values of index SNPs
    * 9 Phenotypes: eGFR, CKD, BUN in ALL, MALE, FEMALE
    * 16 SNPs: index SNP per locus + female-specific index SNP of region 7
3) Regional Association Plot of region 7:
    * panel with all, male, female
4) Credible Set Size vs. Posterior Probability (maybe supplement)

## Supplemental Figures

1) RA-Plots of all regions
2) **Beta-Beta-Plot** of sex interaction analysis (--> see script 02)
3) **Beta-Beta-Plot** of HUNT replication (--> see script 09)
4) Forest Plots (--> see script 12)
5) Co-localization Plots of eQTL results (--> see script 08) 

# Tables (not yet completed)

## Main Tables

1) **Genome-wide significant regions with their respective index SNPs of eGFR and UA** (--> see script T1 for scaffold of table, some information will be added later during literature review and writing process) 
2) **Overlapping regions of eGFR and UA** (--> see script T2)
3) Candidate genes and their function (manual web search for interesting genes)

## Supplemental Tables

--> see script ST

1) Description of Studies (as received from participating studies)
2) Genotyping & Imputation of Studies (as received from participating studies)
3) **Sample Sizes & SNP Numbers, and inflation factor $\lambda$ per phenotype**
4) **Comparisons between the sexes (interaction and co-localization --> see script 02 and 03)**
5) **Cross-phenotype comparision (--> see script 04)**
6) **Genome-wide significant & independent hit per region (--> see script 05)**
7) **Annotation of credible sets** (maybe seperate in a-e for eGFR and UA in their respective settings)
8) **Co-localization with eQTLs** (--> see script 07, only unconditioned so far)
9) **Replication in *HUNT* (--> see script 09)**
10) **Replication of *Graham et al*, *Kanai et al* and *Sakaue et al* results** (--> see script 11)
11) Summary of MR-Mega results (--> see scripts 12)
