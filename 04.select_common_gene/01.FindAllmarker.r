library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(reticulate)
use_python("/dellfsqd2/ST_OCEAN/USER/gaoxiaomin1/software/miniforge3/envs/BT_env_R/bin/python")
anndata <- import("anndata")

args <- commandArgs(TRUE)
rds <- args[1]
prefix <- args[2]
parameter <- args[3]

myset <- readRDS(rds)
#myset <- rds

#data = anndata$read_h5ad(h5ad)
#data <- data$raw$to_adata()
#print("=== 方法1: 直接打印对象 ===")
#data
#mtx = as.matrix(data$X)
#rownames(mtx) = data$obs_names
#colnames(mtx) = data$var_names
#rownames(mtx) = as.character(py_to_r(data$obs_names$values))
#colnames(mtx) = as.character(py_to_r(data$var_names$values))
#myset <- CreateSeuratObject(counts=t(mtx),meta.data=data$obs)
#DefaultAssay(myset) <- 'RNA'

#myset <- NormalizeData(myset,normalization.method = "LogNormalize",scale.factor = 10000,assay = "RNA")
Idents(myset)<-myset@meta.data[[parameter]]
markers <- FindAllMarkers(myset, min.pct = 0.1, logfc.threshold = 0.25)
write.table(format(markers,digits=3),paste0(prefix,".AllMarkers.xls"),sep = '\t', quote = FALSE, row.names =FALSE)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing=TRUE),]
top10mk <- format(markers,digits=3) %>% group_by(cluster) %>% slice(1:30)
write.table(top10mk, paste0(prefix,'top30_Markers.xls'), sep = '\t', quote = FALSE, row.names =FALSE)
