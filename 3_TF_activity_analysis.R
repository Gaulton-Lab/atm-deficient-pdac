# load libraries into R
library(Seurat)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)
library(stringr)
library(presto)
library(RSQLite)
library(openxlsx)
library(dplyr)

# set up parallel computing
register(MulticoreParam(4, progressbar = TRUE))

# read seurat object
setwd("path/to/seurat/object/")
seurat <- readRDS("./AKPC_KPC_SeuratObject.rds")

# read database of TF - target gene regulatory relationships
# compiled from TRRUST 2.0 and ORegAnno 3.0 databases
# contains positively and negatively regulated target genes
load("./TF_and_targets_mm_for_gmt_all.RData")
tf.targets.up <- tfs.mm.for.gmt.all$UP
tf.targets.down <- tfs.mm.for.gmt.all$DOWN

# remove TFs that have less than 7 annotated target genes
tf.targets.up <- lapply(tf.targets.up, function(x) {
  if (length(x) > 6) { return(x) } else { return(NULL)}
})
tf.targets.up <- tf.targets.up[!sapply(tf.targets.up, is.null)]
tf.targets.down <- lapply(tf.targets.down, function(x) {
  if (length(x) > 6) { return(x) } else { return(NULL)}
})
tf.targets.down <- tf.targets.down[!sapply(tf.targets.down, is.null)]

# map name of TFs to motif IDs
tf.id.map <- read.xlsx("tf.motif.ids.xlsx")

# obtain position weight matrices of TF-binding motifs from JASPAR2020 database
matrix_set <- getMatrixByID(x = JASPAR2020, ID = tf.id.map$ID[!is.na(tf.id.map$ID)])

# generate motif object and add to Seurat object
DefaultAssay(seurat) <- "FixPeaks"
motif.matrix <- CreateMotifMatrix(features = StringToGRanges(rownames(seurat)), pwm = matrix_set, genome = 'mm10', score = TRUE, use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = matrix_set)
seurat <- SetAssayData(seurat, assay = 'FixPeaks', layer = 'motifs', new.data = motif.object)

# Calculate enrichment scores of motifs in accessible regions (chromVAR)
DefaultAssay(seurat) <- "FixPeaks"
seurat <- RunChromVAR(object = seurat, genome = BSgenome.Mmusculus.UCSC.mm10, motif.matrix = motif.matrix)

# Subset ductal cells, EMT-like cells, myCAF, iCAF and apCAF from the Seurat object
seurat.subset <- subset(seurat, idents = c("AKPC_Ductal cells", "KPC_Ductal cells", "AKPC_EMT-like cells", "KPC_EMT-like cells", "AKPC_myCAF", "KPC_myCAF", "AKPC_iCAF", "KPC_iCAF", "AKPC_apCAF", "KPC_apCAF"))

# FindAllMarkers: Identify significantly enriched motifs in each celltype (per genotype)
DefaultAssay(seurat.subset) <- "chromvar"
markers_motifs <- FindAllMarkers(seurat.subset, mean.fxn=rowMeans, only.pos = FALSE, test.use = "wilcox", fc.name = "avg_diff",logfc.threshold = -Inf, min.pct = -Inf,min.cells.feature = 0, min.cells.group = 0, return.thresh = 1)
DefaultAssay(seurat.subset) <- "FixPeaks"
# add TF name
motif.names <- markers_motifs$gene
markers_motifs$gene <- str_to_title(ConvertMotifID(seurat.subset, id = motif.names))
# write as csv
write.csv(markers_motifs, "FindAllMarkers.motifs.csv")
# filter significant motifs (p_val_adj < 0.05)
markers_motifs.sig <- markers_motifs[markers_motifs$p_val_adj < 0.05,]
table(markers_motifs.sig$cluster)

# generate a list of motifs with significantly enriched motifs per cell type (significant in AKPC, KPC or both)
motif.list.input <- list(
  Ductal = unique(markers_motifs.sig[markers_motifs.sig$cluster %in% c("AKPC_Ductal cells","KPC_Ductal cells"),]$gene),
  EMT = unique(markers_motifs.sig[markers_motifs.sig$cluster %in% c("AKPC_EMT-like cells","KPC_EMT-like cells"),]$gene),
  myCAF = unique(markers_motifs.sig[markers_motifs.sig$cluster %in% c("AKPC_myCAF","KPC_myCAF"),]$gene),
  iCAF = unique(markers_motifs.sig[markers_motifs.sig$cluster %in% c("AKPC_iCAF","KPC_iCAF"),]$gene),
  apCAF = unique(markers_motifs.sig[markers_motifs.sig$cluster %in% c("AKPC_apCAF","KPC_apCAF"),]$gene)
)

# generate a list of TF target genes per cell type
# positively regulated target genes
targets.up <- list()
for (i in 1:length(motif.list.input)){
  targets.up[[i]] <- tf.targets.up[names(tf.targets.up) %in% motif.list.input[[i]]]
}
names(targets.up) <- names(motif.list.input)
# negatively regulated target genes
targets.down <- list()
for (i in 1:length(motif.list.input)){
  targets.down[[i]] <- tf.targets.down[names(tf.targets.down) %in% motif.list.input[[i]]]
}
names(targets.down) <- names(motif.list.input)

# generate gmt files
# positively regulated target genes
gmt.list.up <- list()
for (i in 1:length(targets.up)) {
  gmt.df <- data.frame(
    GeneSet = names(targets.up[[i]]),
    Description = ".",
    Genes = sapply(targets.up[[i]], function(x) paste(x, collapse = "\t")),
    stringsAsFactors = FALSE
  )
  gmt.list.up[[i]] <- gmt.df
  write.table(gmt.df, paste0("gmt.up.",names(targets.up[i]), ".gmt.txt"),quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}
names(gmt.list.up) <- names(targets.up)
# negatively regulated target genes
gmt.list.down <- list()
for (i in 1:length(targets.down)) {
  gmt.df <- data.frame(
    GeneSet = names(targets.down[[i]]),
    Description = ".",
    Genes = sapply(targets.down[[i]], function(x) paste(x, collapse = "\t")),
    stringsAsFactors = FALSE
  )
  gmt.list.down[[i]] <- gmt.df
  write.table(gmt.df, paste0("gmt.down.",names(targets.down[i]), ".gmt.txt"),quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
names(gmt.list.down) <- names(targets.down)

# FindAllMarkers: Identify differentially expressed genes for each cell type (per genotype)
DefaultAssay(seurat.subset) <- 'SCT'
Idents(seurat.subset) <- seurat.subset$Genotype_celltype
exp.markers <- FindAllMarkers(seurat.subset, only.pos = FALSE, test.use = "wilcox", logfc.threshold = -Inf, min.pct = -Inf,min.cells.feature = 0, min.cells.group = 0, return.thresh = 1)
# write as csv
write.csv(exp.markers, "FindAllMarkers.gex.csv")

# make a list of differentially expressed genes per cell type per genotype
exp.marker.list <- list(
  Ductal.AKPC = exp.markers[exp.markers$cluster %in% c("AKPC_Ductal cells"),],
  Ductal.KPC = exp.markers[exp.markers$cluster %in% c("KPC_Ductal cells"),],
  EMT.AKPC = exp.markers[exp.markers$cluster %in% c("AKPC_EMT-like cells"),],
  EMT.KPC = exp.markers[exp.markers$cluster %in% c("KPC_EMT-like cells"),],
  myCAF.AKPC = exp.markers[exp.markers$cluster %in% c("AKPC_myCAF"),],
  myCAF.KPC = exp.markers[exp.markers$cluster %in% c("KPC_myCAF"),],
  iCAF.AKPC = exp.markers[exp.markers$cluster %in% c("AKPC_iCAF"),],
  iCAF.KPC = exp.markers[exp.markers$cluster %in% c("KPC_iCAF"),],
  apCAF.AKPC = exp.markers[exp.markers$cluster %in% c("AKPC_apCAF"),],
  apCAF.KPC = exp.markers[exp.markers$cluster %in% c("KPC_apCAF"),]
)

# rank genes according to sign(log2FC) * -log10(pvalue)
exp.marker.list <- lapply(exp.marker.list, function(x) {
  x$rank <- sign(x$avg_log2FC) * -log10(x$p_val)
  x <- x[order(x$rank, decreasing = TRUE),]
  return(x)
})

# generate rnk files for each cell type and genotype
rnk.list <- list()
for (i in seq_along(exp.marker.list)) {
  gene_symbols <- exp.marker.list[[i]]$gene
  ranks <- exp.marker.list[[i]]$rank
  max_rank <- max(ranks[is.finite(ranks)])
  min_rank <- min(ranks[is.finite(ranks)])
  ranks[is.infinite(ranks) & ranks == Inf] <- max_rank + 1
  ranks[is.infinite(ranks) & ranks == -Inf] <- min_rank -1
  df <- data.frame(Gene = gene_symbols, Rank = ranks)
  rnk.list[[i]] <- df
  write.table(df, file=paste0(names(exp.marker.list)[i], ".rnk"), sep="\t", row.names=FALSE, col.names=FALSE)
}
names(rnk.list) <- names(exp.marker.list)

##################################################################################################################
# GSEA for TF targets was performed using the Java-based GSEA application (pre-ranked mode)
##################################################################################################################


