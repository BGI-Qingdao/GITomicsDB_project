##
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(stringr)
##
meta.data <- read.table("../../../05.st/04.GenePlot/01.Spatial/01_output/04.Protopterus_annectens_intestine.markers_meta.txt", check.names=F,sep="\t")
meta.data <- subset(meta.data,ann2=="Immune")

##
color_heat <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#ff8000','#e3080b'))
fe <- "CD79B_data"
fe1 <- str_split(fe,"_",2)[[1]]
p1 <- meta.data %>% ggplot(aes_string(x = 'new_x', y = 'new_y', color = 'CD79B_data' )) +
        geom_point(shape = 16, size = 0.1) +
        theme_void() + coord_fixed(ratio=1) +
        theme(legend.position="right", strip.text = element_text(size = 18)) +
        scale_color_gradientn(colors = color_heat(100)) +
        labs(color = fe1)
ggsave("w03.04.Protopterus_annectens_intestine.CD79B_data.png", p1, width = 5, height = 5,dpi = 600, bg="white")
ggsave("w03.04.Protopterus_annectens_intestine.CD79B_data.pdf", p1, width = 5, height = 5)

#
fe <- "TCF7_data"
fe1 <- str_split(fe,"_",2)[[1]]
p1 <- meta.data %>% ggplot(aes_string(x = 'new_x', y = 'new_y', color = 'TCF7_data' )) +
        geom_point(shape = 16, size = 0.1) +
        theme_void() + coord_fixed(ratio=1) +
        theme(legend.position="right", strip.text = element_text(size = 18)) +
        scale_color_gradientn(colors = color_heat(100)) +
        labs(color = fe1)
ggsave("w03.04.Protopterus_annectens_intestine.TCF7_data.png", p1, width = 5, height = 5,dpi = 600, bg="white")
ggsave("w03.04.Protopterus_annectens_intestine.TCF7_data.pdf", p1, width = 5, height = 5)
