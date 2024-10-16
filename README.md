# Analysis of 10X Visium, 10X Chromium and tDISCO RNA-seq data from mouse cortex post stroke (sub-acute)
These scripts detail out analysis of each dataset separately as well as the integration of 10X Visium and Chromium datasets

## Datasets
1) Visium= 4 samples (mouse brain): Day 2, Day 10 and Day 21 post-stroke samples plus a sham (day 2) 
2) 10X Chromium: scRNA-seq analysis on Glast+ astrocytes from around stroke injury area (cortex) from stroke and uninjured mouse
3) tDISCO: scRNA-seq analysis on Day 10 post-stroke mice (n=3), selected GFAP+ astrocytes at defined distances from stroke injury

### Publication:
Scott, E.Y., Safarian, N., Casasbuenas, D.L. et al. Integrating single-cell and spatially resolved transcriptomic strategies to survey the astrocyte response to stroke in male mice. Nat Commun 15, 1584 (2024). https://doi.org/10.1038/s41467-024-45821-y

MiSeq scripts for running tDISCO sequencing can be found [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/MiseqScripts/MSQv3_3FWDReads.zip)
-This Miseq script permits 3 forwards reads, where R1 = R1 (36 bp), R2 = IR1a (8 bp), R3 = IR1b (12 bp), and R4 = R2 (123 bp) (total=179 bp), with the protocol [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/MiseqScripts/Protocol_for_implementing_scripts.docx)  

10X genetables for stroke samples are [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__Stroke.tar.gz) and uninjured samples are [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__Uninjured.tar.gz)
Visium genetables for the D2 post-stroke is [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__V10A06-087-B1.tar.gz), D10 post-stroke is [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__V10A06-087-D1.tar.gz), D21 post-stroke is [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__V10A06-088-C1.tar.gz), and the D2 sham is [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/Faiz_Maryam__V10A06-088-B1.tar.gz).
The normalized tDISCO genetables are organized as the NeuN(red) versus Gfap(green) genetables [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/tDISCO/D10_norm_GvN_norm.xlsx) and the cell number and zone colour matching to the manuscript genetables [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/tDISCO/D10_Norm_Numbered.xlsx)
The proteomics, formatted similar to genetables above, can be found [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/tDISCO/Proteomics_Table_tDISCO.xlsx) and [here](https://data.cyverse.org/dav-anon/iplant/home/eyscott/tDISCO/D10_LFQ_untarg_04_09_forDEP.txt)
