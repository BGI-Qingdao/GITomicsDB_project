library(Seurat)
library(monocle3)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

obj <- readRDS("../Oxy_scvi.integrated.rds")
obj <- subset(obj, subset = Oxy_subcluster %in% c("Parietal-like", "Chief-like"))

mtx <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
meta_data <- obj@meta.data
genes <- data.frame(gene_short_name = rownames(obj), row.names = rownames(obj))

# 创建 monocle3 对象
cds <- new_cell_data_set(mtx, cell_metadata = meta_data, gene_metadata = genes)

# Normalization and PCA
cds <- preprocess_cds(cds, num_dim = 50)

# 去批次（可选）
cds <- align_cds(cds, alignment_group = 'species')

# UMAP
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "Oxy_subcluster", group_label_size = 4) + ggtitle("cds.umap")
p1$layers <- p1$layers[-1]

# Seurat 的 Umap
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]

cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "Oxy_subcluster", group_label_size = 4) + ggtitle("seurat.umap")
p2$layers <- p2$layers[-1]

p <-  p1 | p2
ggsave("reduction_compare.pdf", plot = p, width = 10, height = 5)

# cluster 聚类 这一步没必要
cds <- cluster_cells(cds)
colData(cds)$assigned_cell_type <- meta_data[rownames(cds.embed), ]$Oxy_subcluster


######### 轨迹 #########
cds <- learn_graph(cds)
## 可以调整以下参数来调整分支数量
cds <- learn_graph(cds, learn_graph_control=list(minimal_branch_len=15)) # 数值越大分支越少

p <- plot_cells(cds, color_cells_by = "Oxy_subcluster", label_cell_groups = FALSE, label_branch_points = TRUE, label_groups_by_cluster = FALSE, label_leaves = TRUE,
                label_principal_points = FALSE, labels_per_group = TRUE, label_roots = TRUE, graph_label_size = 4, group_label_size = 4)
p$layers <- p$layers[-1]
ggsave("trajectory_15.pdf", plot = p, width = 10, height = 10)


######### 拟时序 #########
## 手动选择起点
cds <- order_cells(cds)
## 拟时序图
part1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE)
part1$layers <- part1$layers[-1]
ggsave("pseudotime.pdf", plot = part1, width = 6, height = 6)

saveRDS(cds, "cds.rds")

######### 下游基因分析 #########
# 筛选差异基因
cds_gene <- graph_test(cds, neighbor_graph='principal_graph', cores=100)
cds_gene <- cds_gene[order(-cds_gene$morans_I), ]
col_names <- colnames(cds_gene); col_names[5] <- "gene"
write.table(cds_gene, file = "cds_gene.txt", col.names = col_names, sep = "\t", quote = FALSE)


# gene moducle 分析
deg_ids <- row.names(subset(cds_gene, q_value < 0.05)) # 显著性基因
## find_gene_modules 函数是将这些显著基因聚类成细胞中共表达的moducle
gene_module_df <- find_gene_modules(cds[deg_ids, ], resolution = c(10^seq(-6,-1)))
write.table(gene_module_df, "gene_module_df.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## 绘制热图
cell_group_df <- tibble(cell = rownames(colData(cds)), cell_group = colData(cds)$Oxy_subcluster)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)  # 计算基因平均表达量
rownames(agg_mat) <- stringr::str_c("Module", rownames(agg_mat))

pdf("module.heatmap.pdf")
print(pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2"))
dev.off()


# gene 随拟时变化
genes_sig <- cds_gene %>% top_n(n = 100, morans_I) %>% pull(gene_short_name) %>% as.character()
## 显示 celltype
pdf("genes_in_pseudotime.pdf", width = 16, height = 4)
for(i in genes_sig){
  p1 <- plot_genes_in_pseudotime(cds[i, ], color_cells_by = 'Oxy_subcluster', min_expr = 0.2)                                      # 随细胞类型的变化
  p2 <- plot_genes_in_pseudotime(cds[i, ], color_cells_by = 'pseudotime', min_expr = 0.2)                                             # 随拟时序的变化
  p3 <- plot_cells(cds, genes = i, show_trajectory_graph=FALSE, label_cell_groups = FALSE, label_leaves = FALSE, rasterize = TRUE)    # 在UMAP图中的体现
  p4 <- plot_cells(cds, genes = i,                                                                                                    # 
    show_trajectory_graph=TRUE, trajectory_graph_color = 'black', 
    label_cell_groups = FALSE, 
    label_leaves = FALSE, label_roots = FALSE, 
    cell_size = 1, rasterize = TRUE, 
    label_branch_points = FALSE
  ) + scale_color_viridis(option = "inferno")

  p <- p1 | p2 | p3 | p4
  print(p) + plot_layout(nrow = 1, ncol = 4)
}
dev.off()
