# Expression Heatmap
# Rachel Griffard
# 062824

################ Heatmap split by two conditions
library(ComplexHeatmap)

library(circlize)
# col_fun = colorRamp2(c(-2, 0, 2), c("black", "magenta", "white"))
# col_fun = colorRamp2(c(-2, 0, 2), c("#440061", "#28a99e", "yellow"))
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


ncount = read.csv(# normalized count matrix)
gene = c(# genes of interest)
  

gene_int = ncounts %>%
  subset(rownames(ncounts) %in% gene)
  
# scale for heatmap
gene_int_scaled = t(apply(gene_int, 1, scale))
  

# create annotations object
ha = HeatmapAnnotation(tmt = c('t','u','t','u','t','u'), # order this based on samples
                       annotation_legend_param = list(
                         tmt = list(
                           title = "Treatment",
                           at = c("t", "u"), # match above
                           labels = c("100pM RTA-408", "Untreated") # names included in legend
                         )),
                       col = list(tmt = c('t' = 'blue', 'u' = 'darkgreen'))) # colors ordered as appear

Heatmap(as.matrix(gene_int_scaled),
        top_annotation = ha,
        heatmap_legend_param = list(title = 'Scaled expression'),
        column_labels = rep('', 6),
        column_title = 'Metabolic gene expression',
        col = col_fun,
        show_row_dend = FALSE,
        border_gp = gpar(col='black')
)

################ Multiple heatmaps combined
ha = rowAnnotation(tmt = c('t','u','t','u','t','u'),
                   show_legend = FALSE,
                   annotation_legend_param = list(
                     tmt = list(
                       title = "Treatment",
                       at = c("t", "u")
                     )),
                   col = list(tmt = c('t' = 'blue', 'u' = 'darkgreen')))

h1_mat = t(as.matrix(gene_int_scaled)) # genes of interest as cols

h1 = ha + Heatmap(h1_mat,
                  heatmap_legend_param = list(title = 'Scaled expression'),
                  # column_labels = rep('', 19),
                  column_title = 'Metabolic gene expression',
                  col = col_fun,
                  show_column_dend = FALSE,
                  row_split = 2,
                  row_title = c('100pM RTA-408', 'Untreated'),
                  show_column_names = TRUE,
                  border_gp = gpar(col='black')
)

h2_mat = t(as.matrix(cyto_int_scaled))

h2 = Heatmap(h2_mat,
             # heatmap_legend_param = list(title = 'Scaled expression'),
             # column_labels = rep('', 17),
             column_title = 'Cytokine expression',
             col = col_fun,
             show_column_dend = FALSE,
             row_split = 2,
             row_title = c('100pM RTA-408', 'Untreated'),
             show_column_names = TRUE,
             show_heatmap_legend = FALSE,
             border_gp = gpar(col='black')
)

ht_list = h1 + h2

draw(ht_list)