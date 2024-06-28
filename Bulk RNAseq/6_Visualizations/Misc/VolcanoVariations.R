# Volcano with labeled significant
# Rachel Griffard
# 062824

# import file de from differential expression

############## Lines pointing to significant/selected genes
# separate upregulated and downregulated and select top 10
up_de_labels = de %>%
  subset(log2FoldChange > 1) %>%
  arrange(pvalue) %>%
  select(Symbol)
up_de_labels = up_de_labels[1:10,]
down_de_labels = de %>%
  subset(log2FoldChange < -1) %>%
  arrange(pvalue) %>%
  select(Symbol)
down_de_labels = down_de_labels[1:10,]

all_labels = c(up_de_labels, down_de_labels)

# add label marker column to DE matrix
de$label =ifelse(de$Symbol %in% all_labels,
                 de$Symbol, '')

# volcano plot
ggplot(de, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, aes(color = ifelse(log2FoldChange > 1 & pvalue <= 0.05, "Upregulated",
                                          ifelse(log2FoldChange < -1 & pvalue <= 0.05, "Downregulated", "Non-significant"))),
             alpha = 0.4) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray"),
                     name = "") +
  geom_vline(xintercept=c(-1,1), alpha=0.5, color="black", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), alpha=0.5, color="black", linetype = 'dashed') +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  seed = 42.,
                  min.segment.length = 0,
                  segment.curvature = -0.1,
                  segment.ncp = 3, 
                  segment.angle = 20,
                  size = 3) +
  labs(title = "Untreated v 100pM RTA-408", y=expression(-log[10](pvalue)), x=expression(log[2]~Fold~Change)) +
  theme_bw() +
  theme(legend.position = "right") 

# Save the plot as a TIFF file with 300 DPI
ggsave(filename = 'plots/Figure4.png', plot = last_plot(), height = 6, width = 7.5, dpi = 'print')


########## Enlarged points and labels
# separate upregulated and downregulated and select top 10
up_de_labels = de %>%
  subset(log2FoldChange > 1) %>%
  arrange(pvalue) %>%
  select(Symbol)
up_de_labels = up_de_labels[1:10,]
down_de_labels = de %>%
  subset(log2FoldChange < -1) %>%
  arrange(pvalue) %>%
  select(Symbol)
down_de_labels = down_de_labels[1:10,]

other_labels = c('KEAP1')

all_labels = c(up_de_labels, down_de_labels, other_labels)

# add label marker column to DE matrix
de$label =ifelse(de$Symbol %in% all_labels,
                 de$Symbol, '')

# volcano plot
ggplot(de, aes(x = log2FoldChange, y = -log10(pvalue), size = ifelse(!Symbol %in% all_labels, 1, 5))) +
  geom_point(aes(alpha = 0.8, color = ifelse(log2FoldChange > 1 & pvalue <= 0.05, "Upregulated",
                                             ifelse(log2FoldChange < -1 & pvalue <= 0.05, "Downregulated", "Non-significant")))) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray"),
                     name = "") +
  geom_vline(xintercept=c(-1,1), alpha=0.5, color="black", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), alpha=0.5, color="black", linetype = 'dashed') +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf,
                  # box.padding = 0.5,
                  seed = 42.,
                  segment.color = 'transparent',
                  size = 4) +
  labs(title = "Untreated v 100pM RTA-408", y=expression(-log[10](pvalue)), x=expression(log[2]~Fold~Change)) +
  guides(size = FALSE, alpha = FALSE) +
  theme_bw() +
  theme(legend.position = "right") 