print(paste("Start time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))

library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(DoubletFinder)
library(data.table)
library(SoupX)
library(DropletUtils)
options(future.globals.maxSize = 100000 * 1024^3)


args<-commandArgs(T)
if(length(args)!=6){
    print ("Usage:")
    print ("    Rcript soup.umap.R prefix filter_matrix raw_matrix min_nFeature_RNA  max_nFeature_RNA  doublet_rate")
    print ("Note:")
    print ("Aim create clean RDS from matrix")
    q(save = "no", status = 0, runLast = TRUE)
}

######################################
# create folder
######################################
dir.create(args[1]) #prefix
setwd(args[1])

######################################
# cluster filter
######################################
toc <- Read10X(args[2],gene.column=1) #filter_matrix
all <- toc
all <- CreateSeuratObject(all)
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)
matx <- all@meta.data

######################################
# remove background noise
######################################
tod <- Read10X(args[3],gene.column=1) #raw_matrix
tod <- tod[rownames(toc),]
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)
out = adjustCounts(sc,roundToInt=TRUE)
saveRDS(sc,"sc.rds")
DropletUtils:::write10xCounts("soupX_matrix", out,version="3")

######################################
# create clean seurat object
######################################
sample<-args[1]
object_name <- CreateSeuratObject(out,project =sample, min.cells = 3, min.features =  as.numeric(args[4]))
a=length(colnames(object_name))
print(sample)
print(a)

######################################
# basic stats & filter & visualization
######################################
object_name[["percent.mt"]] <- PercentageFeatureSet(object = object_name, pattern = "^MT-|^mt-")
if (sum(object_name@meta.data$percent.mt)== 0){
    object_name@meta.data$percent.mt[1]<-0.000001
}
pdf(paste0(sample,"raw.vln.pdf"),width=15,height=10)
VlnPlot(object_name,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
object_name <- subset(x = object_name, subset = nFeature_RNA > as.numeric(args[4]) & nFeature_RNA < as.numeric(args[5]))
a=length(colnames(object_name))
print(a)
pdf(paste0(sample,"filterGene.vln.pdf"),width=15,height=10)
VlnPlot(object_name,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
object_name <- subset(x = object_name, subset = percent.mt < 10)
pdf(paste0(sample,"filterMT.vln.pdf"),width=15,height=10)
VlnPlot(object_name,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
a=length(colnames(object_name))
print(a)

######################################
# cluster filter clean
######################################
object_name<-  SCTransform(object_name, verbose = FALSE,vars.to.regress="percent.mt")
object_name <- RunPCA(object_name, features = VariableFeatures(object = object_name))
object_name <- FindNeighbors(object_name, dims = 1:30)
object_name <- RunUMAP(object_name, dims = 1:30)
object_name <- FindClusters(object_name, resolution = 0.5)

######################################
# remove doublets
######################################
sweep.res.list_SM <- paramSweep_v3(object_name, PCs = 1:30,sct=TRUE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- object_name@meta.data$RNA_snn_res.0.5
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(as.numeric(args[6])*length(object_name@meta.data$orig.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE,sct=TRUE)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value,sct=TRUE)
object_name@meta.data[,"DF_hi.lo"] <- object_name@meta.data[,10]
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet" & object_name@meta.data[,11] == "Singlet")] <- "Doublet_lo"
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
object_name@meta.data$Doublet <-eval(parse(text = paste0("object_name@meta.data$DF.classifications_",pN_value,"_",pK_value,'_',nExp_poi)))
saveRDS(object_name,"object_name.doublet.rds")
#object_Singlet <- SubsetData(object_name,subset.name='Doublet',accept.value='Singlet')
object_Singlet <- subset(x = object_name, subset = Doublet == "Singlet")
a=length(colnames(object_name))
print(a)

######################################
# basic stats of filter clean singlet
######################################
count_expr <- as.matrix(object_Singlet@assays$RNA@counts)
n_features  <- Matrix::colSums(x = count_expr > 0)
mean(n_features)
cat(paste(sample,",Mean genes per cell after filter,", round(mean(n_features)), "\n",sep=""),file="report.csv",append=T)
cat(paste(sample,",genes after final filter,",nrow(object_Singlet), "\n",sep=""),file="report.csv", append=T)
cat(paste(sample,",cells after final filter,",ncol(object_Singlet), "\n",sep=""),file="report.csv", append=T)
pdf("Vlnplot.pdf",width=15,height=10)
VlnPlot(object_Singlet,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

######################################
# cluster filter clean singlet
######################################
object_Singlet<-SCTransform(object_Singlet, verbose = FALSE,vars.to.regress="percent.mt")
object_Singlet <- RunPCA(object = object_Singlet, npcs = , verbose = FALSE)
pdf("DimHeatmap_30pc.pdf",width=15,height=10)
DimHeatmap(object = object_Singlet, dims = 1:30, cells= 200, balanced = TRUE)
dev.off()
object_Singlet <- RunUMAP(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindNeighbors(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindClusters(object_Singlet, resolution = 0.5)
pdf("sample-umap.pdf",width=15,height=10)
DimPlot(object =object_Singlet, reduction = "umap")
dev.off()
object_Singlet=NormalizeData(object = object_Singlet, assay ="RNA",normalization.method = "LogNormalize",scale.factor = 10000)
object_Singlet=ScaleData(object_Singlet, verbose = FALSE,assay="RNA")
saveRDS(object_Singlet, file = "combined_analysis.rds")

