parser = argparse::ArgumentParser(description = 'Script to integration multi rds data')
parser$add_argument('-i', dest = 'input', help = 'input directory rds files')
parser$add_argument('-o', dest = 'prefix', help = 'output prefix')
opts = parser$parse_args()

library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(future)
library(ggplot2)
library(harmony)

combined = readRDS(opts$input)
DefaultAssay(combined)<-'RNA'
combined = subset(combined,subset=nFeature_RNA>=300)
combined = subset(combined,subset=nFeature_RNA<=4000)
combined = subset(combined,subset=nCount_RNA<=10000)

combined <- NormalizeData(object = combined,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)

combined <- FindVariableFeatures(object = combined,
                                 nfeatures = 2000,
                                 selection.method = "vst")

combined <- ScaleData(object = combined,
                      verbose = FALSE)
combined <- RunPCA(object = combined,
                   npcs = 30,
                   verbose = FALSE)

combined <- RunHarmony(combined,
                       group.by.vars = "orig.ident")

combined <- combined %>% RunUMAP(reduction = "harmony", dims = 1:20) %>% FindNeighbors(reduction = "harmony", dims = 1:20)

combined <- FindClusters(combined, resolution =c(0.1,0.2,0.3,0.4,0.8,1,1.5,2,0.5))
combined <- RunTSNE(combined, reduction = "harmony", dims = 1:20)
saveRDS(combined,paste0(opts$prefix,".harmony.rds"))

pdf(paste0(opts$prefix,".umap.pdf"))
DimPlot(object = combined, reduction = "umap",label = TRUE)
dev.off()
pdf(paste0(opts$prefix,".tsne.pdf"))
DimPlot(object = combined, reduction = "tsne",label = TRUE)
dev.off()


