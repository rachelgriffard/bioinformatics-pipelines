# Seurat Visium Analysis
# Rachel Griffard-Smith
#
# Last Updated: 041525

library(rio)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)

options(future.globals.maxSize = 8000 * 1024^2) 

args = commandArgs(trailingOnly=TRUE)

sample = args[1]
print(paste0('The sample run is ', sample))

print(paste0('The class of sample is ', class(sample)))

fp = paste0('RObjects/',sample)
dir.create(fp)

fp = paste0('plots/',sample)
dir.create(fp)

fp = paste0('plots/',sample, "/Unfiltered")
dir.create(fp)

fp = paste0('plots/',sample, "/Filtered")
dir.create(fp)

hpcdir = paste0("/path/to/files/", sample, "/outs")

print(hpcdir)

obj = Load10X_Spatial(data.dir = hpcdir, bin.size = 8, assay = 'Spatial')

fp = paste0("/path/to/files/", sample, "/Removed_Extra_Tissue.csv")
barcodes = read.csv(fp)

obj = subset(obj, cells = barcodes$Barcode)

obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-") # change to MT if human
obj[["percent.ribosomal"]] = PercentageFeatureSet(obj, pattern = "^Rp") # change to RP if human

fp = paste0('RObjects/', sample,'/01_Unfiltered_', sample, '.RData')
saveRDS(obj, fp)

Assays(obj)
DefaultAssay(obj) = "Spatial.008um"

vln.plot = VlnPlot(obj, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot = SpatialFeaturePlot(obj, features = "nCount_Spatial.008um") + theme(legend.position = "right")
ggsave(vln.plot,filename = paste0('plots/', sample, '/vln_', sample,'.png'), width = 8, height = 8, dpi = 'print')
ggsave(count.plot,filename = paste0('plots/', sample, '/counts_', sample,'.png'), width = 8, height = 8, dpi = 'print')

metadata = obj@meta.data

# QC 
metadata$log10GenesPerUMI = log10(metadata$nFeature_Spatial.008um) / log10(metadata$nCount_Spatial.008um)

# log scaled
metadata %>% 
  ggplot(aes(x=nCount_Spatial.008um)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  ylab("Cell density") +
  xlab("nCount_Spatial.008um (nUMI)") +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/Density_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=nFeature_Spatial.008um)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(100,250,9000)) +
  ggtitle("nGenes per Cell") +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/nGenesPerCell_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=nCount_Spatial.008um, y=nFeature_Spatial.008um, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  # geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/nCountbynFeature', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=percent.mt)) + 
  geom_density(alpha = 0.2) + 
  # scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 10) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/pt_mito_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/GenesPerUMI', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, y=nFeature_Spatial.008um)) +
  geom_point() + 
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Unfiltered/GenesPerUMIvFeature', sample,'.png'), height = 6, width = 6, dpi = 'print')

obj$log10GenesPerUMI = log10(obj@meta.data$nFeature_Spatial.008um) / log10(obj@meta.data$nCount_Spatial.008um)

obj.filt = subset(x = obj,
                         subset = 
                           (nCount_Spatial.008um >= 20) & 
                           (nFeature_Spatial.008um >= 20) & 
                           (nFeature_Spatial.008um <= 9000) &
                           (log10GenesPerUMI >= 0.8) & 
                           (percent.mt < 10))

# QC 
metadata$log10GenesPerUMI = log10(metadata$nFeature_Spatial.008um) / log10(metadata$nCount_Spatial.008um)

# log scaled
metadata %>% 
  ggplot(aes(x=nCount_Spatial.008um)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  ylab("Cell density") +
  xlab("nCount_Spatial.008um (nUMI)") +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/Density_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=nFeature_Spatial.008um)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(100,250,9000)) +
  ggtitle("nGenes per Cell") +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/nGenesPerCell_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=nCount_Spatial.008um, y=nFeature_Spatial.008um, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  # geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/nCountbynFeature', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(x=percent.mt)) + 
  geom_density(alpha = 0.2) + 
  # scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 10) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/pt_mito_', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/GenesPerUMI', sample,'.png'), height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, y=nFeature_Spatial.008um)) +
  geom_point() + 
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave(paste0('plots/', sample, '/Filtered/GenesPerUMIvFeature', sample,'.png'), height = 6, width = 6, dpi = 'print')

DefaultAssay(obj.filt) = "Spatial.008um"
obj.filt = NormalizeData(obj.filt)

fp = paste0('RObjects/', sample,'/02_Normalized_', sample, '.RData')
saveRDS(obj.filt, fp)

DefaultAssay(obj.filt) = "Spatial.008um"
p1 = SpatialFeaturePlot(obj.filt, features = "Rorb") + ggtitle("Rorb expression (8um)")

ggsave(p1,filename = paste0('plots/', sample, '/spatial_', sample,'.png'), width = 8, height = 8, dpi = 'print')

DefaultAssay(obj.filt) = "Spatial.008um"
obj.filt = FindVariableFeatures(obj.filt)
obj.filt = ScaleData(obj.filt)

obj.filt = NormalizeData(obj.filt)
obj.filt = FindVariableFeatures(obj.filt)
obj.filt = ScaleData(obj.filt)
obj.filt = RunPCA(obj.filt)
obj.filt = FindNeighbors(obj.filt)
obj.filt = FindClusters(obj.filt, resolution = 1.0)
obj.filt = RunUMAP(obj.filt, dims = 1:30)

p2 = DimPlot(obj.filt, reduction = "umap", group.by = "seurat_clusters", label = T, pt.size = 1)
p3 = SpatialDimPlot(obj.filt, image.alpha = 0.5, pt.size.factor = 8, shape = 22)

ggsave(p2,filename = paste0('plots/', sample, '/dimPlot_', sample,'.png'), width = 8, height = 8, dpi = 'print')
ggsave(p3,filename = paste0('plots/', sample, '/spatialdimPlot_', sample,'.png'), width = 8, height = 8, dpi = 'print')

fp = paste0('RObjects/', sample,'/03_NormalizedClustered', sample, '.RData')
saveRDS(obj.filt, fp)

markers = FindAllMarkers(obj.filt, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

obj.filt_subset = ScaleData(obj.filt, assay = "Spatial.008um", features = top5$gene)
p4 = DoHeatmap(obj.filt, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
ggsave(p4,filename = paste0('plots/', sample, '/heatmapMarkers_', sample,'.png'), width = 8, height = 8, dpi = 'print')

fp = paste0('RObjects/', sample, '/fullRDS_', sample, '.RData')
save.image(fp)