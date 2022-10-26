# CKDGen Chr X Analyses

Analyses of X-chromosomal SNPs and kidney traits

**Contributor: Markus Scholz, Katrin Horn, Andreas Kühnapfel and Janne Pott**

**Last Updated: 24/10/2022**

Supporting code for the following paper:

* tba

*short description of paper idea*

*contact information*

# Source File

You will need to customize a source file, indicating

* path to R library (please use [R Version 4.x](https://cran.r-project.org/), all necessary packages are listed in the source file)
* path to data (summary statistics, eQTL data)
* path to [PLINK](https://www.cog-genomics.org/plink/2.0/)
* path to [1000 Genomes Phase 3 EUR data](https://www.internationalgenome.org/data-portal/data-collection/phase-3)
* path to [GTEx v8 Data](https://gtexportal.org/home/protectedDataAccess)
* tba

# Included (tbc)

1) Locus definition (check if programmable or just done via excel sheets - maybe shift to *not included*)
2) Interaction analysis (should be no problem, because data sets will not be filtered!) 
3) Coloc between males and females --> unconditioned 
4) Comparison of eGFR / UA hits with CKD / BUN (all phenotype in data available)
5) GCTA - cojo:
    * create input UKBB (in helperscripts)
    * select
    * cond
    * LD via PLINK2 between independent SNPs per locus
6) Credible Set Analyses
7) Coloc with eQTLs
    * all 22 unconditioned loci of script 01
    * perform  extra check for the XXX loci with independent signals using their conditioned statistics
8) Overlapping eGFR and UA loci
    * LD via PLINK2
    * Coloc --> unconditioned; check locus definition --> we only want a pairwise comparision
9) Replication in HUNT of top loci
10) GCTA - heritability (reusing input of 5a)
11) Lookup of candidate SNPs of other publications
12) MR-Mega - check filter!!

## To discuss: Main/Sub Figures:

F1) Main Figure 1: Miami Plot (UA vs eGFR, mark fem-/mal-spec) min(p_fem, p_mal, p_all), known/novel in italic
      * color 1: new gene & new locus
      * color 2: new gene & known locus
      * color 3: known gene & known locus
      * italic, bold, normal for sex-specifity
F2) Main Figure 2: RA region 7 (all, male, female as panel)
F3) Supplemental Figures:
      * RA Plots (male, female, all, all regions?)
      * beta-beta HUNT compare?
      * beta-beta male-female IA?
      * Cred Set size vs CADD variants as in publication Wuttke et al (2019) Fig 4? Cred Set size vs PP per SNP, mark the missense mutations --> maybe main figure?
      * Forest Plots? (MR-Mega color by ethnic origin)
      * Heatmap like phytopaper for eGFR - CKD/BUN comparison? all settings --> 9 traits x 15 loci Panel --> maybe main figure
      * Coloc Plot of eQTL results (UMOD plot in Wuttke?) 
      * Study design (Powerpoint --> Markus)

## To discuss: Main/Sub Tables:

T1) Main Table 1: Top-SNPs (22 loci of script 01) --> typical information of SNP annotation (not tracked in git)
T2) Main Table 2: eGFR vs UA per overlapping regions (results/combination of scripts 08) & add region start/stop
T3) Main Table 3: Candidate Genes and their function (manual search for interesting genes)
T4) Supplemental Tables:
      * Desctiption of Studies
      * Genotyping & Imputation of Studies
      * Study Sample Sizes & SNP Numbers per phenotype
      * Lambda of GWAMA
      * Cred Set Annotation for eGFR & UA in settings with genome-wide sig hits (--> 5 traits)
      * Interaction + Coloc for male-female comparison (check effective N)
      * Coloc eQTLs
      * Replication HUNT
      * CKD/BUN/Gout cross-phenotype comparison
      * Look-up of Graham, X and Y
      * MR-Mega


# Not included (reason)

1) Generate summary statistics per phenotype and setting (result of primary analyses; internal data)
2) Bioinformatic Annotation Pipeline (not yet published from GenStatLeipzig Group) 
3) Trans-ethnic meta-regression (MR-Mega - check with team in Greifswald)
4) Assignment of candidate genes (combination of before mentioned results)

# overview (to be deleted)

* scripts: alle Skripte die wir für die sekundäranalysen brauchen
* helperScripts: alle Skripte, die unsere Datenstruktur sichtbar machen und nicht getrackt werden sollen. Hier können die *not included* sachen mit abgelegt werden
* data: summary statistics, wie sie auch auf Zenodo oder LHA liegen könnten (check col names, unfiltered aber mit flag warum in unserer analyse gefiltert)
* results: alle Ergebnisse des *scripts* Ordners --> daraus kann Markus alles nehmen um das Paper zu schreiben
* temp: alle anderen evtl wichtigen datensets, die man nicht so einfach via Skript erstellen kann (Loci-Def?), aber auch Daten zu eQTL oä

**Namenskonvention**

Alle Skripte bekommen einen YAML header, wo der Autor des Skripts mit drin steht!

Alle Skripte werden nummeriert, wie hier beschrieben! (Unterpunkte wie bei GCTA mit 3_a oder 3_1 damit sortierung erkennbar!)

Plots/Tables können hier mit erzeugt werden, wenn innerhalb der Skripte 1-12 bitte ergänzen

Sourcefile basiert aktuell auf forostar, kann aber auch für jeden anderen server erstellt werden, falls nicht alle nötigen Pakete auf forostar laufen. Bitte im template alle Pakete angeben die man so braucht! (also laufend aktualisieren)
