## Custom resource files created for this project

1. Mouse PDAC celltype marker genes: best 3 marker genes per cell type (from testing larger sets from [Cross-Species Single-Cell Analysis Paper, Elyada et al. (2019)](https://aacrjournals.org/cancerdiscovery/article/9/8/1102/42174/Cross-Species-Single-Cell-Analysis-of-Pancreatic) and [Cross-Dataset Single-Cell Analysis Paper, Yang et al. (2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10137114/).
3. Genes related to TGFb pathway: manually curated from gene set resources, with some basic group annotation (target, ligand, etc).
4. Transcription factor - target gene regulatory relationships: obtained from [TRRUST 2.0](https://academic.oup.com/nar/article/46/D1/D380/4566018?login=false) and [ORegAnno 3.0](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv1203) databases and kindly provided by Medhanie A. Mulaw
5. Mouse motif to TF map: Based on chromVARmotifs `mouse_pwms_v2` PWM object, parsed TF names from the motif names.
6. Mouse motif to TF class, family, and subfamily groups map: Manually reformatted TF family information from [TFClass 2014 resource](http://www.edgar-wingender.de/muTF_classification-1.html).
7. SCREEN V3 mm10 cCREs list: GRanges object made from bed file downloaded from SCREEN website.
