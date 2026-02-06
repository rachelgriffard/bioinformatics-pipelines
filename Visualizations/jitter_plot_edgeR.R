# edgeR Jitter Plots
# Rachel Griffard-Smith
#
# Last updated: 071725

library(tidyverse)
library(ggsignif)

# load data
meta = read.csv('metadata.csv', row.names = 1)
cpm = read.csv('CPM_normalized.csv')
cpm = data.frame(cpm) %>%
  mutate(gene = sapply(str_split(X, "_"), function(x) ifelse(length(x) >= 2, x[2], NA))) %>%
  mutate(gene = ifelse(is.na(gene), remove, gene)) %>%
  mutate(gene = make.unique(gene)) %>%
  column_to_rownames('gene') %>%
  select(-X)

# load genes of interest
goi = readxl::read_xlsx('../../goi.xlsx', sheet = 2, col_names = FALSE)

# subset to genes in cpm
cpm_goi = cpm[rownames(cpm) %in% goi$...1,]

# set up significance from edgeR
## pull in DE frames
UC = read.csv('path/to/de_results.csv', row.names = 1)
CD = read.csv('path/to/de_results.csv', row.names = 1)

## pull out only genes of interest, rename based on comparison
UC = UC[UC$gene %in% goi$...1,] %>% select(gene, PValue) %>% rename(UC = PValue)
CD = CD[CD$gene %in% goi$...1,] %>% select(gene, PValue) %>% rename(CD = PValue)

## create long format data frame, provide comparisons column as formatted below that matches ggplot
pvals_edgeR = UC %>%
  full_join(CD, by = 'gene') %>%
  pivot_longer(cols = -gene) %>%
  rename(edgeR_gene = gene) %>%
  mutate(comparisons = case_when(
    name == 'CD' ~ list(c('CD', 'HC')),
    name == 'UC' ~ list(c('UC', 'HC')),
    TRUE ~ list(NA)
  ))

for (gene in rownames(cpm_goi)) {
  set.seed(42)
  
  cpm_gene = cpm_goi[rownames(cpm_goi) %in% gene,] %>%
    pivot_longer(cols = colnames(cpm_goi),
                 names_to = 'sample',
                 values_to = 'cpm') %>%
    full_join(meta, by = 'sample')
  
  pvals = pvals_edgeR %>% filter(edgeR_gene == !!gene)
  comparison_list = pvals$comparisons
  annotation_list = sprintf("p = %.3g", pvals$value)
  
  annotation_list = case_when(
    pvals$value < 0.001 ~ "***",
    pvals$value < 0.01  ~ "**",
    pvals$value < 0.05  ~ "*",
    TRUE                ~ "n.s."
  )
  
  max_y = max(cpm_gene$cpm, na.rm = TRUE)
  y_pos = seq(max_y * 1.05, length.out = nrow(pvals), by = max_y * 0.1)
  
  ggplot(cpm_gene, aes(group, cpm, shape = group)) +
    geom_jitter(width = 0.2, height = 0, size = 1.8) +
    stat_summary(
      fun = median, geom = "errorbar",
      color = "black", linewidth = 1, width = 0.3,
      fun.min = median, fun.max = median
    ) +
    labs(color = '',x = 'Group', y = 'Counts per\nmillion (CPM)', title = paste0(gene)) +
    guides(shape = 'none') +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = 'bold', size = 18, hjust = 0.5),
      plot.subtitle = element_text(face = 'italic', size = 18, hjust = 0.5),
      axis.title = element_text(face = 'bold', size = 14, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.background = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank()
    ) +
    geom_signif(
      comparisons = comparison_list,
      annotations = c(annotation_list),
      y_position = c(1.1 * max(cpm_gene$cpm), 1.3 * max(cpm_gene$cpm)),
      tip_length = 0.01,
      textsize = 5,
      vjust = -0.1,
      size = 0.6,
      fontface = 'bold'
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  fp = paste0('jitter_plots/', gene, '_jitter_edgeR_final.png')
  ggsave(fp, width = 3, height = 5, dpi = 'print', bg = 'white')
}