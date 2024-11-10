library(Seurat)
packageVersion("Seurat")
update.packages("Seurat")
library(SeuratData)
library(ggplot2)
library(cowplot)
library(harmony)
library(dplyr)

#LoadData("ifnb")



hete = Read10X('./processed_data/mouse/Pax8-Hete-E13-5/outs/filtered_feature_bc_matrix/')
hete = CreateSeuratObject(hete,project = 'Hete')
hete[["percent.mt"]] = PercentageFeatureSet(hete, pattern = "^mt-") ## mouse:mt human:MT
VlnPlot(hete, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hete = subset(hete, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 5)
hete = NormalizeData(hete)
hete = FindVariableFeatures(hete)



homo = Read10X('./processed_data/mouse/Pax8-Homo-E13-5/outs/filtered_feature_bc_matrix/')
homo = CreateSeuratObject(homo,project = 'Homo',min.cells = 3)
homo[["percent.mt"]] = PercentageFeatureSet(homo, pattern = "^mt-") ## mouse:mt human:MT
VlnPlot(homo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
homo = subset(homo, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 5)
homo = NormalizeData(homo)
homo = FindVariableFeatures(homo)
mito = grep(pattern = "^mt-", x = rownames(homo), value = TRUE, ignore.case = TRUE)
homo = homo[!(rownames(homo) %in% mito),]


ctrl = Read10X('./processed_data/mouse/WT-E13-5//outs/filtered_feature_bc_matrix/')
ctrl = CreateSeuratObject(ctrl,project = 'Ctrl',min.cells = 3)
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt-") ## mouse:mt human:MT
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ctrl = subset(ctrl, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 5)
ctrl = NormalizeData(ctrl)
ctrl = FindVariableFeatures(ctrl)
mito = grep(pattern = "^mt-", x = rownames(ctrl), value = TRUE, ignore.case = TRUE)
ctrl = ctrl[!(rownames(ctrl) %in% mito),]

hete$stim = "Hete"
homo$stim = "Homo"
ctrl$stim = "Ctrl"
pax.anchors = FindIntegrationAnchors(object.list = list(ctrl, hete, homo), dims = 1:30)
pax.combined = IntegrateData(anchorset = pax.anchors, dims = 1:30)
DefaultAssay(pax.combined) = "integrated"

# Run the standard workflow for visualization and clustering
pax.combined = ScaleData(pax.combined, verbose = FALSE)
pax.combined = RunPCA(pax.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
pax.combined <- RunUMAP(pax.combined, reduction = "pca", dims = 1:20)
pax.combined <- FindNeighbors(pax.combined, reduction = "pca", dims = 1:20)
pax.combined <- FindClusters(pax.combined, resolution = 0.5)
p1 = DimPlot(pax.combined, reduction = "umap", group.by = "orig.ident")
p2 = DimPlot(pax.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
saveRDS(pax.combined, './Results/combinedData/hete+homo+ctrl_20210606.rds')


DimPlot(pax.combined, reduction = "umap", split.by = "orig.ident")
DefaultAssay(pax.combined) <- "RNA"
markers = 
for(i in 6:19){
  nk.markers = FindConservedMarkers(pax.combined, ident.1 = i,
                                    grouping.var = "orig.ident", verbose = FALSE)
  write.csv(nk.markers, paste0('./Results/combinedData/ctrl+hete+homo/conservedMarkers_',i,'.csv'))
}
head(nk.markers)

DefaultAssay(pax.combined) = "integrated"
markers = FindAllMarkers(pax.combined)
top10 = markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
write.csv(markers, paste0('./Results/combinedData/ctrl+hete+homo/allMarkers.csv'))
DoHeatmap(pax.combined, features = top10$gene) + NoLegend()


FeaturePlot(pax.combined, features = c("Wnt5a","Wnt4", "Wnt7a",
                                       "Pax2","Pax8",
                                       "Hoxa7","Hoxa10","Hoxa11"),
            label=T)



### 0607 先整合hete和homo，再整合ctrl
pax.anchors = FindIntegrationAnchors(object.list = list(homo, hete, ctrl), dims = 1:30)
pax.combined = IntegrateData(anchorset = pax.anchors, dims = 1:30)
DefaultAssay(pax.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pax.combined = ScaleData(pax.combined, verbose = FALSE)
pax.combined = RunPCA(pax.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
pax.combined <- RunUMAP(pax.combined, reduction = "pca", dims = 1:20)
pax.combined <- FindNeighbors(pax.combined, reduction = "pca", dims = 1:20)
pax.combined <- FindClusters(pax.combined, resolution = 0.5)
p1 = DimPlot(pax.combined, reduction = "umap", group.by = "orig.ident")
p2 = DimPlot(pax.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
saveRDS(pax.combined, './Results/combinedData/homo+hete+ctrl/integrated0607.rds')
DimPlot(pax.combined, reduction = "umap", split.by = "orig.ident")

for(i in 0:5){
  nk.markers = FindConservedMarkers(pax.combined, ident.1 = i,
                                    grouping.var = "orig.ident", verbose = FALSE)
  write.csv(nk.markers, paste0('./Results/combinedData/homo+hete+ctrl/conservedMarkers_',i,'.csv'))
}
#head(nk.markers)

DefaultAssay(pax.combined) = "integrated"
markers = FindAllMarkers(pax.combined)
top10 = markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
write.csv(markers, paste0('./Results/combinedData/homo+hete+ctrl/allMarkers.csv'))
DoHeatmap(pax.combined, features = top10$gene) + NoLegend()




#### harmony

harm_combined = CreateSeuratObject(counts = cbind(homo@assays$RNA@counts, 
                                                  hete@assays$RNA@counts, 
                                                  ctrl@assays$RNA@counts), 
                                   project = "harmony_integ") %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = harm_combined@var.genes, npcs = 30, verbose = FALSE)
harm_combined@meta.data$stim <- c(rep("Homo", ncol(homo@assays$RNA@counts)),
                                  rep("Hete", ncol(hete@assays$RNA@counts)), 
                                  rep("Ctrl", ncol(ctrl@assays$RNA@counts)))
p1 = DimPlot(object = harm_combined, reduction = "pca", pt.size = .1, group.by = "stim")
p1 = DimPlot(harm_combined, reduction = "umap", label = TRUE, pt.size = .1,split.by = 'stim')
harm_combined = harm_combined %>% 
  RunHarmony("stim", plot_convergence = TRUE)
p2 = DimPlot(object = harm_combined, reduction = "harmony", pt.size = .1, group.by = "stim")
plot_grid(p1,p2)
harm_combined = RunUMAP(reduction = "harmony", dims = 1:30)


harm_combined = harm_combined %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
p1 = DimPlot(harm_combined, reduction = "umap", label = F, pt.size = .1,group.by = 'stim')
p2 = DimPlot(harm_combined, reduction = "umap", label = TRUE, pt.size = .1)
plot_grid(p1,p2)


DefaultAssay(harm_combined) = "RNA"
for(i in 0:17){
  nk.markers = FindConservedMarkers(harm_combined, ident.1 = i,
                                    grouping.var = "stim", verbose = FALSE)
  write.csv(nk.markers, paste0('./Results/combinedData/harmony/conservedMarkers_',i,'.csv'))
}
#head(nk.markers)

DefaultAssay(harm_combined) = "harmony"
markers = FindAllMarkers(harm_combined)
top10 = markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
write.csv(markers, paste0('./Results/combinedData/harmony/allMarkers.csv'))
DoHeatmap(pax.combined, features = top10$gene) + NoLegend()



### 0608
# 先整合hete和homo
hete_homo = CreateSeuratObject(counts = cbind(homo@assays$RNA@counts[inter_genes,], 
                                              hete@assays$RNA@counts[inter_genes,]), 
                               project = "hete+homo") %>%
  Seurat::SCTransform(verbose = FALSE) %>%
  RunPCA(pc.genes = hete_homo@var.genes, npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30)


# 聪聪师姐说control和hete没有差异，所以做一下整合homo和ctrl的
inter_genes = intersect(rownames(homo@assays$RNA@counts),
                        rownames(ctrl@assays$RNA@counts))
homo_ctrl = CreateSeuratObject(counts = cbind(homo@assays$RNA@counts[inter_genes,], 
                                              ctrl@assays$RNA@counts[inter_genes,]), 
                               project = "homo+ctrl") %>%
  Seurat::SCTransform(verbose = FALSE) %>%
  RunPCA(pc.genes = hete_homo@var.genes, npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
homo_ctrl@meta.data$stim = c(rep("Homo", ncol(homo@assays$RNA@counts)), 
                           rep("Ctrl", ncol(ctrl@assays$RNA@counts)))
p1 = DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1,group.by = 'stim')
homo_ctrl = homo_ctrl %>% 
  RunHarmony("stim", plot_convergence = TRUE,assay.use="SCT")
homo_ctrl = RunUMAP(homo_ctrl,reduction = "harmony", dims = 1:30)
p2 = DimPlot(homo_ctrl, reduction = "umap", pt.size = .1, group.by = "stim")
plot_grid(p1,p2)


homo_ctrl = homo_ctrl %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
p1 = DimPlot(homo_ctrl, reduction = "umap", label = F, pt.size = .1,group.by = 'stim')
p2 = DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1)
plot_grid(p1,p2)

homo_ctrl = RunTSNE(homo_ctrl,reduction = "harmony", dims = 1:30)
DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1,split.by = "stim")
DimPlot(homo_ctrl, reduction = "tsne", label = TRUE, pt.size = .1,split.by = "stim")

options(future.globals.maxSize= 891289600)
DefaultAssay(homo_ctrl) = "RNA"
for(i in 0:17){
  nk.markers = FindConservedMarkers(homo_ctrl, ident.1 = i,
                                    grouping.var = "stim", verbose = FALSE)
  write.csv(nk.markers, paste0('./Results/ctrl+homo/harmony/conservedMarkers_',i,'.csv'))
}

## 0612去掉线粒体和spike-in基因
library(scater)
library(scran)
mito = grep(pattern = "^mt-", x = rownames(homo_ctrl), value = TRUE, ignore.case = TRUE)
homo_ctrl_filtered = homo_ctrl[!(rownames(homo_ctrl) %in% mito),]




## FastMNN
BiocManager::install("bachelor")
library(bachelor)


ctrl_cluster = homo_ctrl[]
FeaturePlot(homo_ctrl, features = c("Wnt7a","Pax8"),split.by = "stim",
            label=T)


tnk.markers = FindConservedMarkers(homo_ctrl, ident.1 = 8,
                                  grouping.var = "stim", verbose = FALSE)
write.csv(tnk.markers, paste0('./Results/ctrl+homo/harmony0617/conservedMarkers_',8,'.csv'))



## 20210619
homo = Read10X('./processed_data/mouse/Pax8-Homo-E13-5/outs/filtered_feature_bc_matrix/')
homo = CreateSeuratObject(homo,project = 'Homo',min.cells = 3)
homo[["percent.mt"]] = PercentageFeatureSet(homo, pattern = "^mt-") ## mouse:mt human:MT
homo = subset(homo, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 5)
mito = grep(pattern = "^mt-", x = rownames(homo), value = TRUE, ignore.case = TRUE)
homo = homo[!(rownames(homo) %in% mito),]

ctrl = Read10X('./processed_data/mouse/WT-E13-5//outs/filtered_feature_bc_matrix/')
ctrl = CreateSeuratObject(ctrl,project = 'Ctrl',min.cells = 3)
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt-") ## mouse:mt human:MT
ctrl = subset(ctrl, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 5)
mito = grep(pattern = "^mt-", x = rownames(ctrl), value = TRUE, ignore.case = TRUE)
ctrl = ctrl[!(rownames(ctrl) %in% mito),]

inter_genes = intersect(rownames(homo@assays$RNA@counts),
                        rownames(ctrl@assays$RNA@counts))
homo_ctrl = CreateSeuratObject(counts = cbind(homo@assays$RNA@counts[inter_genes,], 
                                              ctrl@assays$RNA@counts[inter_genes,]), 
                               project = "homo+ctrl") %>%
  Seurat::SCTransform(verbose = FALSE) %>%
  RunPCA(pc.genes = hete_homo@var.genes, npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
homo_ctrl@meta.data$stim = c(rep("Homo", ncol(homo@assays$RNA@counts)), 
                             rep("Ctrl", ncol(ctrl@assays$RNA@counts)))
p1 = DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1,group.by = 'stim')
homo_ctrl = homo_ctrl %>% 
  RunHarmony("stim", plot_convergence = TRUE,assay.use="SCT")
homo_ctrl = RunUMAP(homo_ctrl,reduction = "harmony", dims = 1:30)
p2 = DimPlot(homo_ctrl, reduction = "umap", pt.size = .1, group.by = "stim")
plot_grid(p1,p2)


homo_ctrl = homo_ctrl %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
p1 = DimPlot(homo_ctrl, reduction = "umap", label = F, pt.size = .1,group.by = 'stim')
p2 = DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1)
plot_grid(p1,p2)

homo_ctrl = RunTSNE(homo_ctrl,reduction = "harmony", dims = 1:30)
DimPlot(homo_ctrl, reduction = "umap", label = TRUE, pt.size = .1,split.by = "stim")
DimPlot(homo_ctrl, reduction = "tsne", label = TRUE, pt.size = .1,split.by = "stim")

options(future.globals.maxSize= 891289600)
DefaultAssay(homo_ctrl) = "RNA"
for(i in 0:16){
  nk.markers = FindConservedMarkers(homo_ctrl, ident.1 = i,
                                    grouping.var = "stim", verbose = FALSE)
  write.csv(nk.markers, paste0('./Results/ctrl+homo/harmony0619/conservedMarkers_',i,'.csv'))
}


FeaturePlot(homo_ctrl, features = c("Pax8","Wnt7a"),
            #split.by = "stim",coord.fixed=T,
            #min.cutoff = 0,max.cutoff = 2,
            coord.fixed = T,
            label=T)+scale_color_gradientn(colors = c('lightgrey', 'blue'),limits=c(0,3))

homo_cluster = homo_ctrl[,which(homo_ctrl@meta.data$stim == 'Homo')]
ctrl_cluster = homo_ctrl[,which(homo_ctrl@meta.data$stim == 'Ctrl')]
p1 = FeaturePlot(homo_cluster,features = 'Pax8',label = T)+
  scale_color_gradientn(colors = c('lightgrey', 'blue'),limits=c(0,3))
p2 = FeaturePlot(ctrl_cluster,features = 'Pax8',label = T)+
  scale_color_gradientn(colors = c('lightgrey', 'blue'),limits=c(0,3))
plot_grid(p1,p2)


saveRDS(homo_ctrl,'./Results/ctrl+homo/harmony0619/ctrl_homo0619.rds')
