# Boxplot
# Rachel Griffard
# 062824

##### Boxplot with manual p-value
b1 = c(#genes of interest here as string)

b1_count = ncounts %>%
  subset(rownames(ncounts) %in% b1) %>%
  rownames_to_column('gene') %>%
  pivot_longer(!gene, names_to = 'patient')

b1_count$trt = ifelse(b1_count$patient == 'P1_100pM' | b1_count$patient ==  'P2_100pm' |
                        b1_count$patient == 'P3_100pm', '100pM RTA-408', 'Untreated')
b1_count$trt = factor(b1_count$trt, labels = c('Untreated', '100pM RTA-408'))
b1_count$trt = relevel(b1_count$trt, ref = 'Untreated')

de[de$Symbol == 'HMOX1',]

p1 = b1_count %>%
  filter(gene == 'HMOX1')
pl1 = ggboxplot(p1, x = 'trt', y = 'value', fill = 'trt', palette = c('white', 'blue')) +
  geom_bracket(xmin = 'Untreated', xmax = '100pM RTA-408', 
               y.position = 150,
               label = "p-value < 0.0001") +
  ggtitle('HMOX1') +
  ylab('CPM') + xlab('Treatment') +
  theme_bw() +
  guides(fill = FALSE)
ggsave('./plots/Figure3_HMOX1.png', pl1, height = 5, width = 5, dpi = 'print')