# Seurat Visium Analysis - Filter and replot
# Rachel Griffard-Smith
#
# Last Updated: 041525

# filter based on prior script

library(rio)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(tidyverse)

options(future.globals.maxSize = 8000 * 1024^2) # to increase access to mem for Seurat in HPC

args = commandArgs(trailingOnly=TRUE)

dir.create('integration')
dir.create('integration/plots')
dir.create('integration/RObjects')

merged.obj = readRDS('integration/RObjects/02_mergedSeurat.RDS') # load merged samples

merged.obj$log10GenesPerUMI = log10(merged.obj@meta.data$nFeature_Spatial.008um) / log10(merged.obj@meta.data$nCount_Spatial.008um)

# filter based on QC plots
filtered_seurat = subset(x = merged.obj,
                         subset = 
                           (nCount_Spatial.008um >= 20) & 
                           (nFeature_Spatial.008um >= 20) & 
                           (nFeature_Spatial.008um <= 9000) &
                           (log10GenesPerUMI >= 0.8) & 
                           (percent.mt < 10))
filtered_seurat

saveRDS(filtered_seurat, 'integration/RObjects/03_filteredSeurat.RDS')

# redo plots from before
# QC 
metadata = filtered_seurat@meta.data
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
ggsave('integration/plots/QC/Filtered/NCells.png', height = 6, width = 6, dpi = 'print')

metadata.long = gather(metadata, type, value, 2:5)

metadata.long %>%
  ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) + 
  geom_violin() + 
  # geom_point() +
  facet_wrap(vars(type), ncol = 2, scales = "free")  +
  theme_bw()
ggsave('integration/plots/QC/Filtered/Violin.png', height = 6, width = 6, dpi = 'print')

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
ggsave('integration/plots/QC/Filtered/CellDensity_nUMI.png', height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_Spatial.008um, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(100,250,9000)) +
  ggtitle("nGenes per Cell") +
  theme_bw()
ggsave('integration/plots/QC/Filtered/nGenesperCell.png', height = 6, width = 6, dpi = 'print')

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
ggsave('integration/plots/QC/Filtered/nCount_Spatial.008umvnFeature_Spatial.008um.png', height = 6, width = 6, dpi = 'print')

metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  # scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 10) +
  theme_bw()
ggsave('integration/plots/QC/Filtered/percent_mt.png', height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme_bw()
ggsave('integration/plots/QC/Filtered/log10GenesPerUMIvDensity.png', height = 6, width = 6, dpi = 'print')

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, y=nFeature_Spatial.008um, color = log10(nCount_Spatial.008um))) +
  geom_point() + 
  facet_wrap(~orig.ident) +
  theme_bw()
ggsave('integration/plots/QC/Filtered/log10GenesPerUMIvnFeature_Spatial.008um.png', height = 6, width = 6, dpi = 'print')