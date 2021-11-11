# Final Project Outline
## Title
Differential Gene Expression in TCGA within Stage A Chronic Lymphocytic Leukemia comparing Progression and Stable using DeSEQ2
## Author
Ke Wang
## Overview of project
I will identify differentially expressed genes between Chronic Lymphocytic Leukemia for progression and stable patients. 
This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. For this analysis, I'll use the TCGA cohort and have identified 205 Ht-seq counts files for tumors that fit within my cohort with 37 progression patients and 102 stable patients. Within the analysis, I will control for age and gender.
## Data
I will use the data from https://dcc.icgc.org/search. Examining clinical data, there are 139 tumor samples, and 37 are defined by me as progression patients and 102 are identified as stable patients. The classification standard is the category data in the clinical data. The specific files are available in donor.tsv.
## Milestone 1
Due Date: Thursday November 11th

I will integrate the data of eligible progression and stable patients into a tsv. Data fully loaded into vignette through HT-SEQ steps. I will complete an entire first draft of analysis analyzed through the vignette.
## Milestone 2
Due Date: Tuesday November 23rd

An initial completion of vignette. I will complete an entire first draft of analysis analyzed through the vignette, and begin to describe the analysis results.

## Deliverable
Due Date: December 3rd

A complete repository in Github with clear documentation and description of my analysis and results.
