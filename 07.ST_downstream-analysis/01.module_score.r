library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

args <- commandArgs(T)

mymodulescore <- function(myobj, genelist){
        GenesetScore <- AddModuleScore(
                object = myobj,
                nbin = 5,
                features = list(as.vector(unlist(genelist))),
                name = "score")
        return(GenesetScore)
}

obj_input <- readRDS(args[1])

DefaultAssay(obj_input) <- "RNA"
obj_input <- NormalizeData(obj_input)

gene_input <- read.table(args[2],header=F)
gene_input <- gene_input$V1

outdir <- args[3]

score_obj <- mymodulescore(obj_input, gene_input)

mydata <- FetchData(score_obj, vars = c("new_x","new_y","score1"))

color_heat <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#ff8000','#e3080b'))

plot <- ggplot(mydata, aes(x = new_x,y = new_y, colour = score1)) +
                geom_point(shape = 16, size = 0.1) +
                theme_void() + coord_fixed(ratio=1) +
                labs(col = "Tfh Module Score") +
                scale_color_gradientn(colors = color_heat(100))
ggsave(paste0(outdir,".pdf"), plot, width= 5, height=5)
ggsave(paste0(outdir,".png"), plot, width= 5, height=5,bg="white",dpi= 600)
