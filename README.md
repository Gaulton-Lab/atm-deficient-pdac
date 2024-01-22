# Mouse PDAC 

Github: [https://github.com/Gaulton-Lab/mouse-pdac]

## Overview: 

This repository contains notebooks designed to run several QC metrics on the mouse PDAC data. There are 3 notebooks: merging and doublet removal, findMarkers and fGSEA, and CellChat.

This project contains two libraries consisting of two samples each. There are two genotypes: AKPC and KPC. The samples were first processed through the multiome pipeline (created by Hannah Mummey) that includes three pipelines: R Script #1, Python Script, and R Script #2. You can find the github repo for the multiome pipeline here: [https://github.com/Gaulton-Lab/multiome-pipeline]

## Required Dependencies

#### R Packages

- Seurat (v4.3.0.1)
- Signac (v1.10.0)
- EnsDb.Mmusculus.v79
- SoupX
- Harmony
- ggplot2
- ggpubr
- reticulate
- data.table
- dplyr
- GenomeInfoDb
- Cellranger
- CellChat
- Scrublet 
