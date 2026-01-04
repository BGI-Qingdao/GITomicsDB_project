library(Seurat)
library(harmony)
library(clustree)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggsci)
args<-commandArgs(T)

#
obj <- readRDS(args[1])
obj@meta.data <- subset(obj@meta.data, select=c("orig.ident","nCount_RNA","nFeature_RNA","species","organ","species.full","cellname.raw","x","y","z","new_x","new_y"))
#
obj <- SCTransform(obj, assay = "RNA", verbose = FALSE)
obj <- RunPCA(obj)
obj <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca")
obj <- FindNeighbors(obj, dims = 1:30,reduction = "harmony")
obj <- FindClusters(obj, verbose = FALSE, resolution = c(seq(0.1,1,0.2)))
obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony")

#
saveRDS(obj, file = paste0(args[2], ".rds"))
write.table(obj@meta.data,file=paste0(args[2], ".meta.txt"), sep ="\t", quote=FALSE, row.names =TRUE, col.names =TRUE)

#
cluster_Palette <- unique(c(pal_lancet("lanonc",alpha=1)(9), "#EFC000FF", "#8F7700FF", "#FFCD004C", "#E18727FF", "#5C88DA4C", "#EE4C97FF", "#20854EFF", pal_lancet("lanonc",alpha=0.6)(9), pal_simpsons("springfield")(16), brewer.pal(12, "Set3"), pal_igv("default")(51), pal_ucscgb("default")(26), colorRampPalette((pal_npg("nrc")(9)))(45)))

#
p1 <- clustree(obj)+coord_flip()
ggsave(paste0(args[2],".clustree.pdf"), p1, w=24, h=16)

plist <- list()
n=0
for(i in seq(0.1,1,0.2)){
        n=n+1
        res <- paste0("SCT_snn_res.",i)
        plist[[n]]<-DimPlot(obj, reduction = "umap", group.by = res, label = T, cols = cluster_Palette) #+ NoLegend()
}
ggsave(paste0(args[2],".dimplot.pdf"), patchwork::wrap_plots(plist,ncol = 5), width = 25, height = 5)
