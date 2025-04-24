# volcano plot functions

library(tidyverse)

volcano_pretty = function(de_results, label_genes = NA,
                          p_threshold = 0.05, logFC_threshold = 2,
                          title,
                          gene_col = 'Gene',
                          logFC_col = 'logFC', pval_col = 'pval',
                          adj_pval_col = 'adj_pval',
                          pretty_type) {
  
  if (pretty_type == 'large_labeled') {
    
    de_results$label_gene_col = ifelse(de_results[[gene_col]] %in% label_genes,
                                       de_results[[gene_col]], '')
    de_results$size = !de_results[[gene_col]] %in% label_genes
    
    de_results = de_results %>% mutate(reg = ifelse(de_results[[logFC_col]] > logFC_threshold & de_results[[pval_col]] <= p_threshold, "Upregulated",
                                             ifelse(de_results[[logFC_col]] < -logFC_threshold & de_results[[pval_col]] <= p_threshold, "Downregulated", "Non-significant")))
    
    de_results$pval = de_results[[pval_col]]
    de_results$logFC = de_results[[logFC_col]]
    
    ggplot(de_results, aes(x = logFC, y = -log10(pval))) +
      geom_point(aes(size = size, color = reg),
                 alpha = 0.8) +
      scale_size_manual(values = c(3,1)) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray"),
                         name = "") +
      geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_hline(yintercept=-log10(p_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_text_repel(aes(label = label_gene_col),
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      seed = 42.,
                      min.segment.length = 0,
                      segment.curvature = -0.1,
                      segment.ncp = 3, 
                      segment.angle = 20,
                      segment.color = 'transparent',
                      size = 5,
                      fontface = 'bold') +
      labs(title = title, y=expression(-log[10](p~Value)), x=expression(log[2]~Fold~Change)) +
      guides(size = 'none', alpha = 'none') +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = 'bold', hjust = 0.5, size = 18),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 12))
  }
  
  if (pretty_type == 'large_significant') {
    
    de_results$label_gene_col = ifelse(de_results[[gene_col]] %in% label_genes,
                                       de_results[[gene_col]], '')

    
    de_results = de_results %>% mutate(reg = ifelse(de_results[[logFC_col]] > logFC_threshold & de_results[[pval_col]] <= p_threshold, "Upregulated",
                                                    ifelse(de_results[[logFC_col]] < -logFC_threshold & de_results[[pval_col]] <= p_threshold, "Downregulated", "Non-significant")))
    de_results$size = de_results$reg == 'Non-significant'
    
    de_results$pval = de_results[[pval_col]]
    de_results$logFC = de_results[[logFC_col]]
    
    ggplot(de_results, aes(x = logFC, y = -log10(pval))) +
      geom_point(aes(size = size, color = reg),
                 alpha = 0.8) +
      scale_size_manual(values = c(3,1)) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray"),
                         name = "") +
      geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_hline(yintercept=-log10(p_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_text_repel(aes(label = label_gene_col),
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      seed = 42.,
                      min.segment.length = 0,
                      segment.curvature = -0.1,
                      segment.ncp = 3, 
                      segment.angle = 20,
                      segment.color = 'transparent',
                      size = 5,
                      fontface = 'bold') +
      labs(title = title, y=expression(-log[10](p~Value)), x=expression(log[2]~Fold~Change)) +
      guides(size = 'none', alpha = 'none') +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = 'bold', hjust = 0.5, size = 18),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 12))
  }
  
  if (pretty_type == 'lots_o_labels') {
    
    de_results$label_gene_col = ifelse(de_results[[gene_col]] %in% label_genes,
                                       de_results[[gene_col]], '')
    
    
    de_results = de_results %>% mutate(reg = ifelse(de_results[[logFC_col]] > logFC_threshold & de_results[[pval_col]] <= p_threshold, "Upregulated",
                                                    ifelse(de_results[[logFC_col]] < -logFC_threshold & de_results[[pval_col]] <= p_threshold, "Downregulated", "Non-significant")))
    de_results$size = de_results$reg == 'Non-significant'
    
    de_results$pval = de_results[[pval_col]]
    de_results$logFC = de_results[[logFC_col]]
    
    ggplot(de_results, aes(x = logFC, y = -log10(pval))) +
      geom_point(aes(color = reg),
                 alpha = 0.4) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray"),
                         name = "") +
      geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_hline(yintercept=-log10(p_threshold), alpha=0.5, color="black", linetype = 'dashed') +
      geom_text_repel(aes(label = label_gene_col),
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      seed = 42.,
                      min.segment.length = 0,
                      segment.curvature = -0.1,
                      segment.ncp = 3, 
                      segment.angle = 20,
                      size = 3) +
    labs(title = title, y=expression(-log[10](p~Value)), x=expression(log[2]~Fold~Change)) +
      guides(size = 'none', alpha = 'none') +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = 'bold', hjust = 0.5, size = 18),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 12))
  }
}

volcano_pretty(d1_rna, label_genes, title = '24 Hour Post-Activation',
               logFC_col = 'RNA_logFC', pval_col = 'RNA_pvalue', pretty_type = 'lots_o_labels')
