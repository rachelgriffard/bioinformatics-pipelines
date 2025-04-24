# Visium HD - integration with sc atlas
# Atlas: https://pubmed.ncbi.nlm.nih.gov/39405347/
# Rachel Griffard-Smith
# Last updated: 042325

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(spacexr)

options(future.globals.maxSize = 8000 * 1024^2)

vis = readRDS('../integration/RObjects/01a_listSeurat_filtered.RDS')

ref = readRDS('reference/muto-reference-sc.rds') # path to reference sc data in seurat format

# set up reference
Idents(ref) = "celltype"
counts = ref[["RNA"]]$counts
cluster = as.factor(ref$celltype)
nUMI = ref$nCount_RNA
levels(cluster) = gsub("/", "-", levels(cluster))
cluster = droplevels(cluster)
reference = Reference(counts, cluster, nUMI)

saveRDS(ref, 'RObjects/reference_spacexr.rds')

# set up query
for (i in 1:length(vis))  { 
  sample = unique(vis[[i]]$orig.ident)
  
  coords = GetTissueCoordinates(vis[[i]][["slice1.008um"]], which = "centroids")
  query.counts = GetAssayData(vis[[i]], assay = "Spatial.008um", layer = "counts")
  rownames(coords) = coords$cell
  coords$cell = NULL
  query = SpatialRNA(coords, query.counts, colSums(query.counts))
  
  fp = paste0('RObjects/query_', sample, '_spacexr.rds')
  saveRDS(query, fp)
  
  # create the RCTD object and run with doublet detection
  RCTD = create.RCTD(query, reference, max_cores = 16)
  RCTD = run.RCTD(RCTD, doublet_mode = "doublet")
  
  fp = paste0('RObjects/predictions_', sample, '_spacexr.rds')
  saveRDS(RCTD, fp)
}