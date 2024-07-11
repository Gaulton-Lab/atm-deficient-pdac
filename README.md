# Multimodal Single Cell Investigation of ATM-Deficient PDAC 

## Overview: 

This repository contains notebooks and R scripts used to process and analyze single cell sequencing data (10x Multiome) from tumor samples from KPC and AKPC (KPC + ATM knock out) mouse models of pancreatic ductal adenocarcinoma. Briefly, these notebooks were used to combine multiple single cell sequencing libraries into a single map, apply QC metrics, cluster and identify cell types, and perform downstream analyses using both the snRNA-seq and snATAC-seq. There are 5 notebooks, each containing code necessary for separate tasks: 
1. Merge libraries and perform doublet removal
2. Run differential gene expression analysis with findMarkers and fGSEA (snRNA-seq based)
3. Predict cell-cell interactions with CellChat (snRNA-seq based)
4. Identify transcription factor - target gene regulatory relationships: Infer transctiption factor-associated accessibility using ChromVAR and prepare input for GSEA of transcription factor target genes
5. Calculate and compare accessibility of TF binding motifs with ChromVAR (snATAC-seq based)
6. Construct gene regulatory networks (GRNs) with Pando (snRNA-seq and snATAC-seq based)

This project contains two libraries consisting of two pooled samples each. There are two genotypes: AKPC and KPC, each represented as one library. The samples were first processed through a 10x Multiome processing pipeline that includes three major steps: identification and removal of empty droplets, conversion of snATAC-seq data to 5kb windows, and background correction. You can find the github repo for this pipeline here: [https://github.com/Gaulton-Lab/multiome-pipeline]

## Required Dependencies

As different tools require different dependencies and sometimes even different versions of the same dependencies, we had to run different analyses in different environments. To see which dependencies were required for each notebook, see the final section with `sessionInfo()` information.

## Citation

[insert paper citation and link once it's published]
