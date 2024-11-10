install.packages('Seurat')
install.packages('cowplot')
install.packages("ggplot2")
library(Seurat)
library(cowplot)
library(ggplot2)



featureplot = function(data, gene){
  d1_cluster = data[,which(data@meta.data$stim == 'Homo')]
  d2_cluster = data[,which(data@meta.data$stim == 'Ctrl')]
  limit_up = max(max(Seurat::GetAssayData(d1_cluster)[gene,]),
                 max(Seurat::GetAssayData(d2_cluster)[gene,]))
  p1 = FeaturePlot(d1_cluster,features = gene,label = T)+
    scale_color_gradientn(colors = c('lightgrey', 'blue'),limits=c(0,3))+
    ggtitle(paste("Homo",gene,'_'))
  p2 = FeaturePlot(d2_cluster,features = gene,label = T)+
    scale_color_gradientn(colors = c('lightgrey', 'blue'),limits=c(0,3))+
    ggtitle(paste("Ctrl",gene,'_'))
  plot_grid(p1,p2)
}


homo_ctrl = readRDS('./Results/ctrl+homo/harmony0619/ctrl_homo0619.rds')

featureplot(homo_ctrl, 'Pax8')
VlnPlot(homo_ctrl, 'Pax8',split.by = 'stim')
DotPlot(homo_ctrl, features = c('Pax8'),split.by = "stim", cols = c("blue", "red"))
