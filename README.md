# Analysis of 10X Visium, 10X Chromium and tDISCO RNA-seq data from mouse cortex post stroke (sub-acute)
These scripts detail out analysis of each dataset separately as well as the integration of 10X Visium and Chromium datasets

## Datasets
1) Visium= 4 samples (mouse brain): Day 2, Day 10 and Day 21 post-stroke samples plus a sham (day 2) 
2) 10X Chromium: scRNA-seq analysis on Glast+ astrocytes from around stroke injury area (cortex) from stroke and uninjured mouse
3) tDISCO: scRNA-seq analysis on Day 10 post-stroke mice (n=3), selected GFAP+ astrocytes at defined distances from stroke injury

### Publication:
Integrating single-cell and spatially resolved transcriptomic strategies to survey astrocytes in response to stroke. Scott EY, Safarian N, Lozano Casasbuenas D, Tockovska T, Dryden M, Ali S, Peng J, Daniele E, Tripathy S, Yuzwa S, Wheeler A, Faiz M. Submitted.

MiSeq scripts for running tDISCO sequencing can be found [here](https://de.cyverse.org/data/ds/iplant/home/eyscott/MiseqScripts?type=folder&resourceId=8691cac2-7f59-11ee-a8fe-90e2ba675364)
-This Miseq script permits 3 forwards reads, where R1 = R1 (36 bp), R2 = IR1a (8 bp), R3 = IR1b (12 bp), and R4 = R2 (123 bp) (total=179 bp)  
