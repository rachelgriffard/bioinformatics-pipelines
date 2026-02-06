# Qiagen IPA Analysis - Kalyani
# Rachel Griffard-Smith
#
# 071025

library(tidyverse)
library(ggrepel)

setwd() # set to where save plots

# canonical pathways

canon = read.delim('name_of_file.txt', skip = 2)

canon$z.score = as.numeric(canon$z.score)
canon$logp = as.numeric(canon$X.log.p.value.)
canon$ratio = as.numeric(canon$Ratio)

canon_toppz = canon %>%
  filter(logp > -log10(0.01), abs(z.score) > 4) %>%
  mutate(Pathway = factor(Ingenuity.Canonical.Pathways,
                          levels = rev(Ingenuity.Canonical.Pathways)))

ggplot(canon_toppz, aes(x = logp, y = Pathway)) +
  geom_col(aes(fill = z.score), width = 0.6, color = 'black', size = 0.7) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "gray80", name = "Z-score"
  ) +
  labs(
    x = expression(-log[10](p-value)),
    y = NULL,
    title = "Top COMPARISON Canonical Pathways",
    subtitle = 'abs(z-score) > 4 & p-value < 0.01'
  ) +
  geom_vline(xintercept = -log10(0.05), lwd = .7, linetype = 'dashed') +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = 'italic', hjust = 0.5),
    panel.border = element_rect(size = 1.5, color = 'black', fill = NA),
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_blank(),
    axis.ticks = element_line(size = 1)
  )
ggsave('comparison_signifpz_sortp.png', height = 8, width = 18, dpi = 'print', bg = 'white')

# upstream regulators

regs =  read.delim('regs_file.txt', skip = 2)

regs = regs %>%
  mutate(
    logp = -log10(p.value.of.overlap),
    zscore = as.numeric(Activation.z.score),
    Upstream.Regulator = fct_reorder(Upstream.Regulator, logp)
  )

top_by_p = regs %>% 
  arrange(desc(logp)) %>% 
  slice_head(n = 5)

top_by_z = regs %>% 
  arrange(desc(abs(zscore))) %>% 
  slice_head(n = 5)

top_inhibited = regs %>%
  filter(zscore < 0) %>%
  arrange(abs(zscore) %>% desc()) %>%
  slice_head(n = 5)

top_labels = bind_rows(top_by_p, top_by_z, top_inhibited) %>% 
  distinct()

ggplot(regs, aes(x = zscore, y = logp, label = Upstream.Regulator)) +
  geom_point(aes(color = Predicted.Activation.State), alpha = 0.8) +
  geom_text_repel(
    data = top_labels,
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("Activated" = "#D55E00", "Inhibited" = "#0072B2", "Not-significant" = "gray")) +
  labs(
    x = "Activation Z-score",
    y = expression(-log[10](p-value)),
    title = "Upstream Regulators"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = 'italic', hjust = 0.5),
    panel.border = element_rect(size = 1.5, color = 'black', fill = NA),
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_blank(),
    axis.ticks = element_line(size = 1)
  )
ggsave('upstreamRegs.png', height = 8, width = 10, dpi = 'print', bg = 'white')

# Diseases and Functions

dis = read.delim('diseasesFunctions.txt', skip = 2)

dis = dis %>%
  filter(!is.na(Activation.z.score)) %>%
  mutate(logp = -log10(p.value),
         zscore = Activation.z.score)

dis %>%
  filter(p.value < 0.05, abs(Activation.z.score) > 2) %>%
  slice_max(-log10(p.value), n = 25) %>%
  mutate(Diseases.or.Functions.Annotation = fct_reorder(Diseases.or.Functions.Annotation, -log10(p.value))) %>%
  ggplot(aes(x = -log10(p.value), y = Diseases.or.Functions.Annotation, fill = Activation.z.score)) +
  geom_col() +
  scale_fill_gradient2(high = "#D55E00", mid = "white", low = "#0072B2", midpoint = 0, name = "Z-score") +
  labs(
    title = "IPA Diseases and Functions",
    x = expression(-log[10](p-value)), y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(size = 1.5, color = 'black', fill = NA),
    panel.grid = element_blank()
  )
ggsave('diseaseFunctions.png', height = 8, width = 10, dpi = 'print', bg = 'white')