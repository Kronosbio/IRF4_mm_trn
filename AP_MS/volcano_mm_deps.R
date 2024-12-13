rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)
library(readxl)
library(ggrepel)
interactors_tot <- read.csv("~/ash_2024/interactors_apms_MM_deps.csv")

DE <- read.csv("~/ash_2024/20230519_KB450_MM1S_24h_24h_DE.csv")

DE$colors <- ""
DE$colors[DE$Gene %in% c("EP300","IRF4","IKZF3","IKZF1",interactors_tot$GENE)] <- "MM" 

DE$labels <- ""
DE$labels[DE$Gene %in% c("IRF4","IKZF1","IKZF3","ARID3B")] <- DE$Gene[DE$Gene %in% c("IRF4","IKZF1","IKZF3","ARID3B")] 
dat_select <- (DE %>% filter(colors == "MM"))

p <- DE %>%
  ggplot(aes(x = log2FC, y = logadjpval, color = colors, label = labels)) +
  geom_point(size = 3, alpha = 0.5) +
  xlim(-1, 1) +
  ylim(0, 4.5) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(max.overlaps = 10000) +
  labs(x = "log2FC", y = "-log10(adjpval)") +
  ggtitle("MM1S 24h vs 24h DE") +
  scale_color_manual(values = c("grey", "#54278F"))

p<- p + geom_point(data = dat_select,aes(x = log2FC, y = logadjpval, size = 3.5))

pdf("./volcano_mm_deps.pdf",height = 10, width = 10)
p
dev.off()