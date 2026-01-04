##par
parser = argparse::ArgumentParser(description = "RCTD plot")
parser$add_argument('-i', '--input', help = 'input RDS')
parser$add_argument('-o', '--out', help = 'out directory')
parser$add_argument('-n', '--name', dest = 'name', help = 'sample name')
opts = parser$parse_args()

##lib
library(spacexr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(cowplot)

##input
myRCTD <-readRDS(opts$input)

##results
results <- myRCTD@results
spatialRNA <- myRCTD@spatialRNA

##coords
coords <- spatialRNA@coords
coords$orig.ident <- opts$name

######################################
#doublet_mode results
######################################
results_df <- results$results_df
my_table <- merge(coords,results_df,by=0)
rownames(my_table)<-my_table[,1]
my_table<-my_table[,-1]
doublet_weights <-results$weights_doublet #cell type proportions doublet_mode
colnames(doublet_weights)<-c("first_type_weight","second_type_weight")
pie_d <- merge(my_table,doublet_weights,by=0)
rownames(pie_d) <- pie_d[,1]
pie_d <- pie_d[,-1]
pie_d$first_type <- as.character(pie_d$first_type)
pie_d$second_type <- as.character(pie_d$second_type)
write.table(pie_d,file=paste0(opts$out, '/', opts$name,".doublet_result.xls"),sep="\t",quote=F)
meta.data <- pie_d

#
meta.data$id<-rownames(meta.data)
x<- meta.data[,c(5,13,15)]
colnames(x)<-c("type","weight","id")
y<-meta.data[,c(6,14,15)]
colnames(y)<-c("type","weight","id")
xy<-rbind(x,y)
m<-dcast(xy,id~type,value.var ="weight")
rownames(m)<-m$id
m <- m[,-1]
m[is.na(m)] <- 0
p <- meta.data[,c(1,2)]
#
m1 <- m
m1$x = p[match(rownames(m1),rownames(p)),]$x
m1$y = p[match(rownames(m1),rownames(p)),]$y
data_dw <- m1
write.table(data_dw,file=paste0(opts$out, '/', opts$name,".doublet_weights.xls"),sep="\t",quote=F)
#meta.data <- subset(meta.data,spot_class!="reject")

######################################
#full_mode results
######################################
norm_weights = normalize_weights(results$weights) #normalize the cell type proportions to sum to 1
norm_weights <- as.data.frame(norm_weights)
pie <- merge(coords,norm_weights,by=0)
rownames(pie)<-pie[,1]
data_fw <- pie[,-1]
len <- length(colnames(data_fw))
data_fw$max_value <- apply(data_fw[,4:len], 1, max)
data_fw$main_celltype <- apply(data_fw[,4:len], 1, function(t) colnames(data_fw[,4:len])[which.max(t)])
write.table(data_fw,file=paste0(opts$out, '/', opts$name, ".full_weights.xls"), sep="\t", quote=F)

#####################################
#color
#####################################
#cluster_Palette <- unique(c(pal_npg("nrc")(10), pal_d3("category10",alpha=0.9)(10), pal_lancet("lanonc",alpha=1)(9), pal_simpsons("spri      ngfield")(16), brewer.pal(12, "Set3"), pal_igv("default")(51), pal_ucscgb("default")(26)))
cluster_Palette <- c('#A6CEE3FF','#1F78B4FF', '#B2DF8AFF', '#33A02CFF' ,'#FB9A99FF' ,'#E31A1CFF' ,'#FDBF6FFF' ,'#FF7F00FF','#CAB2D6FF' ,      '#6A3D9AFF' ,'#FFFF99FF' ,'#B15928FF','#803620FF','#E5D2DD',  '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A'      , '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1      E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175','#FFA500FF'      , '#800080FF', '#008000FF', '#FFFF00FF', '#00FFFFFF',  '#FF00FFFF', '#8B4513FF', '#A9A9A9FF', '#DEB887FF', '#5F9EA0FF',  '#7FFF00FF', '#      D2691EFF', '#6495EDFF', '#FFF8DCFF', '#DC143CFF',  '#00FFFFFF', '#00008BFF', '#008B8BFF', '#B8860BFF', '#A0522DFF',  '#4B0082FF', '#556B      2FFF', '#FF8C00FF', '#9932CCFF', '#8B008BFF',  '#E9967AFF', '#8FBC8FFF', '#4169E1FF', '#FFA07AFF', '#DAA520FF',  '#808000FF', '#BDB76BFF      ', '#800000FF', '#483D8BFF', '#2F4F4FFF',  '#228B22FF', '#FAF0E6FF', '#FFD700FF', '#ADD8E6FF', '#F08080FF',  '#E6E6FAFF', '#708090FF', '      #FFF0F5FF', '#F0E68CFF', '#EEEEEEFF',  '#7B68EEFF', '#F5FFFAFF', '#FFE4E1FF', '#00FF00FF', '#98FB98FF',  '#FFEBCDFF', '#D3D3D3FF', '#F8F      9FAFF', '#FFB6C1FF', '#DA70D6FF')

######################################
#doublet_mode spatial plot
######################################
#first_type
plot1 <- ggplot(meta.data, aes_string(x = 'x', y = 'y', color = 'first_type')) +
                geom_point(shape = 19, size = 0.1) +
                theme_void() +
                 coord_fixed(ratio=1) +
                theme(legend.position = 'right') +
                scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1))+
                labs(title = paste0(opts$name,' doublet first')) + theme(plot.title=element_text(size=15,face ="bold"))
ggsave(paste0(opts$out, '/', opts$name, '.doublet_first.pdf'), plot1, width = 10, height = 10)

#second_type
plot2 <- ggplot(meta.data, aes_string(x = 'x', y = 'y', color = 'second_type')) +
                geom_point(shape = 19, size = 0.1) +
                theme_void() +
                coord_fixed(ratio=1) +
                theme(legend.position = 'right') +
                scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1)) +
                labs(title = paste0(opts$name,' doublet second')) + theme(plot.title=element_text(size=15,face ="bold"))
ggsave(paste0(opts$out, '/', opts$name, '.doublet_second.pdf'), plot2, width = 10, height = 10)

#scatterpies
data <- data_dw
len <- length(colnames(data))
cellname <- colnames(data)[-c(len-1,len)]
plt <- ggplot() + scatterpie::geom_scatterpie(data = data, aes(x = x, y = y), col = cellname, color = NA, pie_scale = 0.1) +
        coord_fixed(ratio = 1) +
        scale_fill_manual(values = cluster_Palette ) +
        theme_void() +
        guides(fill=guide_legend(ncol=1)) +
        labs(title = paste0(opts$name,' doublet scatterpies')) +
        theme(plot.title=element_text(size=18,face ="bold"))
#ggsave(paste0(opts$out, '/', opts$name,'.doublet_scatterpies.pdf'), plt, width = 20, height = 20)
ggsave(paste0(opts$out, '/', opts$name,'.doublet_scatterpies.png'), plt, width = 20, height = 20, dpi = 300,bg="white")

#############################################
#full plot
############################################

#weights
data <- data_fw
len <- length(colnames(data))
plots <- list()
for (i in 4:(len-2)) {
  fe <- colnames(data)[i]
  plots[[i-3]] <- data %>% ggplot(aes_string(x = 'x', y = 'y', color = data[,i])) +
        geom_point(shape = 16, size = 0.1) +
        coord_fixed(ratio=1) +
        theme_classic() +
        theme(legend.position="right", title = element_text(size=15)) +
        labs(color = "weights", title= fe) +
        scale_color_gradientn(colors = colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#ff8000','#e3080b'))(100))
}
num <- len-5
plt <- plot_grid(plotlist = plots[1:num], ncol = 6)
ggsave(paste0(opts$out, '/', opts$name, ".full_weights.pdf"), plt, width = 48, height = (num/6+1)*8)

#scatterpies
cellname <- colnames(data)[-c(1:3,len-1,len)]
plt <- ggplot() + scatterpie::geom_scatterpie(data = data, aes(x = x, y = y), col = cellname, color = NA, pie_scale = 0.1) +
        coord_fixed(ratio = 1) +
        scale_fill_manual(values = cluster_Palette ) +
        theme_void() +
        guides(fill=guide_legend(ncol=1)) +
        labs(title = paste0(opts$name,' full scatterpies')) +
        theme(plot.title=element_text(size=18,face ="bold"))
#ggsave(paste0(opts$out, '/', opts$name,'.full_scatterpies.pdf'), plt, width = 20, height = 20)
ggsave(paste0(opts$out, '/', opts$name,'.full_scatterpies.png'), plt, width = 20, height = 20, dpi = 300,bg="white")

#main_celltype
plt <- ggplot(data, aes_string(x = 'x', y = 'y', color = 'main_celltype')) +
                geom_point(shape = 19, size = 0.1) +
                theme_void() +
                coord_fixed(ratio=1) +
                guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1))+
                scale_color_manual(values = cluster_Palette) +
                labs(title = paste0(opts$name,' full main')) + theme(plot.title=element_text(size=15,face ="bold"))
ggsave(paste0(opts$out, '/', opts$name, '.full_main.pdf'), plt, width = 10, height = 10)
