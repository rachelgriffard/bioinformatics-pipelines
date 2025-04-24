library(optima)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(stringr)
library(patchwork)
library(rio)


all_categories = c('WT', 'HET', 'HOM', 'Missing')
colors = structure(rev(c('black', 'red', 'blue', 'gray')), names = c('WT', 'HET', 'HOM', 'Missing'))

# load('RObjects/data-loaded.RData')

sample_num = as.integer(args[1])
sample_num = sample_num + 1
samples = read.csv('samples.csv')
sample = samples[sample_num, 1]
file.name = samples[sample_num, 2]

sample = get(file.name)

paste0('Loading sample ', sample, '...')

fp = paste0('collab_Tapestri/collab_', sample, '.dna.h5')
sample = readHdf5(directory = fp,
              sample.name = sample,
              omic.type = "DNA")




fp = paste0('RObjects/collab_', sample, '.rds')
saveRDS(sample, fp)

paste0('Sample ', sample, ' loaded and saved!')

paste0('Filtering sample ', sample, '...')

sample.filtered = filterVariant(sample,
                            vaf.ref = 5,
                            vaf.hom = 95,
                            vaf.het = 20)

fp = paste0('RObjects/collab_', sample, 'filtered.rds')
saveRDS(sample.filtered, fp)

paste0('Filtered sample ', sample, 'saved!')

paste0('Creating heatmap for sample ', sample, '...')

sample.gt = data.frame(sample.filtered@gt.mtx)
rownames(sample.gt) = rownames(getDNAmtx(sample.filtered))
colnames(sample.gt) = colnames(getDNAmtx(sample.filtered))

sample.gt[sample.gt == 0] = 'WT'
sample.gt[sample.gt == 1] = 'HET'
sample.gt[sample.gt == 2] = 'HOM'
sample.gt[sample.gt == 3] = 'Missing'

fp = paste0('output/Heatmap/', sample, '_heatmap.png')
png(file = fp, width = 750, height = 750)
ComplexHeatmap::Heatmap(as.matrix(sample.gt), col = colors,
                        show_row_names = FALSE, column_names_rot = 90, column_names_gp = gpar(fontsize = 10),
                        column_title = sample, column_title_gp = gpar(fontsize = 15, fontface = 'bold'))
dev.off()