library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
args = commandArgs(T)

##
df <- read.table(args[1], head=T, row.names=1, sep="\t",check.names=F)
colnames(df)[colnames(df) == as.character(args[3]) ] <- "annotation"
prefix <- paste0(args[2], "-", args[3])
des <- c("annotation",
"ATP6+",
"B cells",
"BEST4+",
"Ciliated",
"Cycling",
"Endothelial",
"Enterocytes",
"Enteroendocrine",
"Epi-Cpl-1",
"Epi-Pan-1",
"Epi-Pse-1",
"Erythrocytes",
"Foveolar",
"Goblet",
"Goblet-like",
"ImmCycling",
"Mesenchymal",
"Mucus",
"Myeloid",
"Oxynticopeptic",
"Oxynticopeptic-like",
"T cells",
"Tuft"
)
df <- df [, names(df) %in% des]
df$id <- row.names(df)
df[is.na(df)] <- 0
write.table(df,file= paste0(args[2], ".1.weights.xls"), quote=F, sep="\t", row.names= T)
cellprops_info <- melt(df,id.vars=c("id", "annotation"))
names(cellprops_info) <- c("spot_id", "mol_niche", "name", "value")
cellprops_info$mol_niche <- paste("domain", cellprops_info$mol_niche, sep="_")
write.table(cellprops_info, file= paste0(args[2], ".2.spot_weights.xls"), quote=F, sep="\t", row.names=F)

##
ct_description <-  cellprops_info %>%
     na.omit() %>%
    select(-spot_id) %>%
    group_by(name) %>%
    nest() %>%
    mutate(wres = map(data, function(dat) {

      niches <- dat$mol_niche %>%
        unique() %>%
        set_names()

      map(niches, function(g) {

        test_data <- dat %>%
          mutate(test_group = ifelse(.data[["mol_niche"]] == g,
                                     "target", "rest")) %>%
          mutate(test_group = factor(test_group,
                                     levels = c("target", "rest")))

        wilcox.test(value ~ test_group,
                    data = test_data,
                    alternative = "greater") %>%
          broom::tidy()
      }) %>% enframe("mol_niche") %>%
        unnest()

    }))
##
wilcox_types <- ct_description %>%
    dplyr::select(wres) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::mutate(adj_pval = p.adjust(p.value,method ="fdr")) %>%
    dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
    dplyr::mutate(sign = ifelse(adj_pval < 0.001, "*", ""))

##
ct_median_desc <- ct_description %>%
    dplyr::select(data) %>%
    unnest() %>%
    group_by(name, mol_niche) %>%
    summarise(median_prop = median(value)) %>%
    mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop))

# Bind both dataframes
niche_car_df <- left_join(wilcox_types, ct_median_desc) %>% na.omit()
write_csv(niche_car_df, file = paste0(args[2], ".3.result.csv"))

# Order based on clustering of both rows and columns
ct_median_desc_mat <- niche_car_df %>% dplyr::select(name, mol_niche, scaled_median_prop) %>%
    pivot_wider(names_from = mol_niche, values_from = scaled_median_prop) %>%
    column_to_rownames("name") %>%
    as.matrix()

ct_sign_desc_mat <- niche_car_df %>% dplyr::select(name, mol_niche, sign) %>%
    pivot_wider(names_from = mol_niche, values_from = sign) %>%
    column_to_rownames("name") %>%
    as.matrix()

library(ComplexHeatmap)
color_heat <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#ff8000','#e3080b'))(100)
niche_car_plt <- Heatmap(ct_median_desc_mat,
                         column_title = args[4],
                         name = "scaled median comp.",
                         rect_gp = gpar(col = "grey", lwd = 0.1),
                         col = color_heat,
                         cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf(ct_sign_desc_mat[i, j]), x, y, gp = gpar(fontsize = 10))}
                        )
pdf(paste0(args[2], ".3.result.pdf"), height = 5, width = 7)
draw(niche_car_plt)
dev.off()
