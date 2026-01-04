library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle3)

cds <- readRDS("../cds.rds")

p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, cell_size = 0.5) + theme_void()
p$layers <- p$layers[-1]

# 取图的坐标范围
build <- ggplot_build(p)
x_range <- build$layout$panel_params[[1]]$x.range
y_range <- build$layout$panel_params[[1]]$y.range

arrow_length <- diff(x_range)/5

# 箭头的箭头样式
arrow_style <- arrow(length = unit(0.1, "inches"), type = "closed")

p2 <- p +
  # 自定义x轴箭头线（底部）
  geom_segment(aes(x = x_range[1], y = y_range[1], xend = x_range[1] + arrow_length, yend = y_range[1]),
               inherit.aes = FALSE, color = "black", size = 0.7, arrow = arrow_style) +
  # 自定义y轴箭头线（左侧）
  geom_segment(aes(x = x_range[1], y = y_range[1], xend = x_range[1], yend = y_range[1] + arrow_length),
               inherit.aes = FALSE, color = "black", size = 0.7, arrow = arrow_style) +
  
  # 添加x轴和y轴刻度线
  theme(
    title = element_text(face = "bold"),
    
    axis.ticks = element_blank(), 
    
    axis.text.x.bottom = element_blank(),
    axis.text.y.left = element_blank(),
    
    # axis.title = element_text(size = 12, face = "bold"),
    # axis.text = element_text(size = 10),
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) + 
  geom_text(aes(x = x_range[1], y = y_range[1]-0.2), label = "UMAP_1", hjust = 0, vjust = 1) +  # x轴标题靠左下方一点
  geom_text(aes(x = x_range[1] - 0.8, y = y_range[1]), label = "UMAP_2", angle = 90, hjust = 0, vjust = 1)  # y轴标题靠左下角

png("DimPlot_pseudotime.png", width = 1500, height = 1500, res = 300)
print(p2)
dev.off()
