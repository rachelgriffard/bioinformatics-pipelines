# Seurat Visium Analysis - Integrate samples & QC
# Rachel Griffard-Smith
#
# Last Updated: 041525

library(rio)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(tidyverse)

options(future.globals.maxSize = 8000 * 1024^2) # to increase access to mem for Seurat in HPC

args = readLines('samplenames.txt') # one sample name per line

dir.create('integration')
dir.create('integration/plots')
dir.create('integration/RObjects')

# read in separate samples through R data

obj.list = list()

for (i in 1:length(args)) { # load samples
  fp = paste0("path/to/sample/files/", args[i], "/outs") # hpc
  #fp = paste0("/hawk/work/biostat/r816g589/collaborations/tran/032525/rawdata",  args[i], "/outs") #hawk
  obj = Load10X_Spatial(data.dir = fp, bin.size = 8, assay = 'Spatial')
   
   
  fp = paste0("path/to/sample/files/", args[i], "/Removed_Extra_Tissue.csv")
  
  # Check if the file exists
  if (file.exists(fp)) { # if created manual fidicual removal in Loupe save as Removed_Extra_Tissue.csv
    # If the file exists, read it and subset the object
    barcodes = read.csv(fp)
    obj = subset(obj, cells = barcodes$Barcode)
  } else {
    # If the file does not exist, print a message and continue
    message("File does not exist: ", fp, ". Continuing without subsetting.")
  }
  
  obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-") # change to MT if human
  obj[["percent.ribosomal"]] = PercentageFeatureSet(obj, pattern = "^Rp") # change to RP if human
  
  obj$orig.ident = rep(args[i])
  obj$group = rep(toupper(stringr::str_extract(args[i], "[^\\-]+$")))
  
  obj.list[[i]] = obj
}

saveRDS(obj.list, 'integration/RObjects/01_listSeurat.RDS')

# merge samples
merged.obj = merge(obj.list[[1]], y = c(obj.list[[2]],obj.list[[3]],obj.list[[4]],obj.list[[5]],obj.list[[6]]),
                   add.cell.ids = names(obj.list))

saveRDS(merged.obj, 'integration/RObjects/02_mergedSeurat.RDS')


metadata = merged.obj@meta.data

# QC 
metadata$log10GenesPerUMI = log10(metadata$nFeature_Spatial.008um) / log10(metadata$nCount_Spatial.008um)

metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/NCells.png', height = 6, width = 6, dpi = 'print')

metadata.long = gather(metadata, type, value, 2:5)

metadata.long %>%
  ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) + 
  geom_violin() + 
  # geom_point() +
  facet_wrap(vars(type), ncol = 2, scales = "free")  +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/Violin.png', height = 6, width = 6, dpi = 'print')

# log scaled
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_Spatial.008um, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  ylab("Cell density") +
  xlab("nCount_Spatial.008um (nUMI)") +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/CellDensity_nUMI.png', height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_Spatial.008um, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(100,250,9000)) +
  ggtitle("nGenes per Cell") +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/nGenesperCell.png', height = 6, width = 6, dpi = 'print')

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
ggsave('integration/plots/QC/Unfiltered/nCount_Spatial.008umvnFeature_Spatial.008um.png', height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  # scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 10) +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/percent_mt.png', height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/log10GenesPerUMIvDensity.png', height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, y=nFeature_Spatial.008um, color = log10(nCount_Spatial.008um))) +
  geom_point() + 
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave('integration/plots/QC/Unfiltered/log10GenesPerUMIvnFeature_Spatial.008um.png', height = 6, width = 6, dpi = 'print')

