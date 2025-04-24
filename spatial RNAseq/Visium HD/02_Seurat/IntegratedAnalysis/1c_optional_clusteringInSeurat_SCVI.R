
library(Seurat)
library(SeuratDisk)
library(tidyverse)
# Read in and UMAP SCVI data

options(future.globals.maxSize = 8000 * 1024^2) 

Convert("RObjects/04c_IntegratedData_scVI_default.h5ad", dest = "h5Seurat", assay = "integrated", overwrite = T)
int.obj <- LoadH5Seurat("RObjects/04c_IntegratedData_scVI_default.h5seurat", misc = F, meta.data = F)
# metadata = read.csv("RObjects/metadata_scvi_default.csv", row.names = 1)
# int.obj@meta.data = metadata


merged.obj = readRDS("RObjects/02_mergedSeurat.RDS")
meta = merged.obj@meta.data
int.obj@meta.data = meta

saveRDS(int.obj, 'RObjects/05_Integrated_scVI_metadata.rds')

int.obj = RunUMAP(int.obj, reduction = "scVI", reduction.name = "umap.scvi", dims = 1:10)

saveRDS(int.obj, 'RObjects/05b_Integrated_scVI_metadata_UMAP.rds')

p = DimPlot(int.obj, reduction = "umap.scvi", group.by = "orig.ident", label = TRUE) + 
  ggtitle("UMAP Clustering with scVI") +
  theme_minimal()
ggsave('plots/SCVI/dimPlot_origIdent.png', p, width = 8, height = 8, dpi = 'print', bg = 'white')

int.obj = FindNeighbors(int.obj, reduction = "scVI", dims = 1:10)
int.obj = FindClusters(int.obj, resolution = 1, cluster.name = "scvi_louvain_clusters_res1.0")

saveRDS(int.obj, 'RObjects/05c_Integrated_scVI_metadata_UMAP_clustered.rds')

p2 = DimPlot(int.obj, reduction = "umap.scvi", group.by = "scvi_louvain_clusters_res1.0", label = TRUE) + 
  ggtitle("UMAP Clustering with scVI") +
  theme_minimal()
ggsave('plots/SCVI/dimPlot_cluster.png', p2, width = 8, height = 8, dpi = 'print', bg = 'white')
