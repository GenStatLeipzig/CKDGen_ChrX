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
7) Coloc with eQTLs --> conditioned 
8) Overlapping eGFR and UA loci
    * LD via PLINK2
    * Coloc --> unconditioned?
9) Replication in HUNT of top loci
10) GCTA - heritability (reusing input of 5 a)
11) Genetic Correlation (not yet done, compare / read publication of Yong Li)
12) Lookup of candidate genes of other publications
    
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

