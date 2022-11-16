# FullCyc2: Bacterial community dynamics explain carbon mineralization and assimilation in soils of different land-use history

This repository contains the code and some data used in the analysis of the DNA stable isotope probing study examining bacterial assimilation of carbon from 5 distinct substrates across cropland, old-field, and forest soils. 

The paper to which this analysis belongs is titled "Bacterial community dynamics explain carbon mineralization and assimilation in soils of different land-use history" (DOI:10.1111/1462-2920.16146) and can be cited as:

Barnett, S.E., Youngblut, N.D. & Buckley, D.H. (2022) Bacterial community dynamics explain carbon mineralization and assimilation in soils of different land-use history. Environmental Microbiology, 24:5230– 5247. Available from: https://doi.org/10.1111/1462-2920.16146


## Files

### R_analyses
This directory contians all R based analyses, all writen in rmarkdown (.rmd) with GitHub formatted markdown output (.md) for easy reading and figure visualization.
* **GCMS_analysis:** Analysis of CO2 mineralization data.
* **MWHRSIP_analysis:** Running multiple-window high-resolution SIP to identify 13C-labeled bacterial OTUs.
* **Incorporators:** Initial visualization and examination of the 13C-labeled bacterial OTUs (incorporators).
* **Incorporator_lifehistories:** Examination of the rRNA operon copy number of the incorporators.
* **Unfrac_analysis:** Analyses of the raw microcosm bacterial communities. DNA used here was extracted from the microcosm soils but not run through a gradient, hence "unfractionated". These microcosms had all 5 substrates added.
* **Enrichments:** Analysis of the smaller enrichment microcosm bacterial communities. These microcosms had only 1 substrate added.


### Data
This directory contains the processed data used in the analyses described above. Raw sequencing reads can be found on the NCBI SRA under BioProject PRJNA686389.
* **GCMS_calc_data.txt:** CO2 mineralization data after calculating concentrations in GCMS_analysis.
* **fullcyc2physeq.RDS:** The complete OTU dataset in phyloseq format including the metadata, taxonomy table, OTU table.
* **fullcyc2.bacteria.cogent.tree:** The complete bacterial OTU phylogenetic tree in Newick format.
* **otusn.pick.fasta:** The consens sequences of all OTUs.
* **fullcyc2_l2fc_testoutput.rds:** The Log2 fold change table output from MW-HR-SIP, indicating which OTU are 13C-labeled in each treatment (padj < 0.05).
* **final_incorp.bacteria.unique_seqs.csv:** the output of paprica indicating the predicted 16S rRNA operon copy number for each incorporator.
