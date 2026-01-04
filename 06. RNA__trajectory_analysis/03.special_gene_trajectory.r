library(Seurat)
library(monocle3)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

cds <- readRDS("../cds.rds")

### draw specific gene
gene_pseudotime <- function(gene){
  pdf(paste0(gene, "_pseudotime.pdf"), height = 3, width = 3)
  print(plot_genes_in_pseudotime(cds[gene, ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
          NoLegend() + 
          xlab("") + 
          ylab(""))
  dev.off()
  
  png(paste0(gene, "_pseudotime.png"), height = 2000, width = 2500, res = 1000)
  print(plot_genes_in_pseudotime(cds[gene, ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
          NoLegend() + 
          xlab("") + 
          ylab(""))
  dev.off()
}

genes <- c("JUNB", "SOCS3", "JUN", "JUND", "ATF3", "FOS", "LYZ", "PPP1R15A", "FTL", "PHLDA2", "EGR1", "CPA2", "DNAJB1", "DUSP6", "SRRM1", 
           "MDK", "EIF3F", "PLCD4", "HSP90AA1", "COX4I1", "MB", "SELENOP", "CHIA", "PGC", "CKB", "NME1", "YBX3", "FKBP11", "GAPDH", "GNMT", 
           "PGA3", "RPL22L1", "RPLP0", "CSRP1", "NDUFB3", "RPL4", "COX7B", "MPP1", "LDHB", "COX5B", "P4HB", "RPS2", "PHPT1", "CISD1", "ATP5PB", 
           "UQCR10", "MIF", "COX7A2", "ATP5F1B", "RPL6", "VDAC3", "SRSF6", "ATP4A")
for (i in genes) {
  gene_pseudotime(i)
}


### RPS3
pdf("RPS3_pseudotime.pdf", height = 6, width = 3)
print(plot_genes_in_pseudotime(cds["RPS3", ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
        NoLegend() + 
        xlab("") + 
        ylab(""))
dev.off()

png("RPS3_pseudotime.png", height = 4000, width = 2500, res = 1000)
print(plot_genes_in_pseudotime(cds["RPS3", ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
        NoLegend() + 
        xlab("") + 
        ylab(""))
dev.off()

### 
pdf("RPS3_pseudotime.pdf", height = 6, width = 3)
print(plot_genes_in_pseudotime(cds["RPS3", ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
        NoLegend() + 
        xlab("") + 
        ylab(""))
dev.off()

png("RPS3_pseudotime.png", height = 4000, width = 2500, res = 1000)
print(plot_genes_in_pseudotime(cds["RPS3", ], color_cells_by = 'pseudotime', min_expr = 0.2) + 
        NoLegend() + 
        xlab("") + 
        ylab(""))
dev.off()
