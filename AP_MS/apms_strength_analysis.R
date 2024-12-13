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
library(patchwork)

p300_sig <- read.csv("~/ash_2024/apms/p300_3R24M_noMito_average_score.csv")
IRF4_sig <- read.csv("~/ash_2024/apms/IRF4_3R24M_noMito_average_score.csv")

p300_all <- read.csv("~/ash_2024/apms/MM1S_3R24M_p300_noMito_FP_byMSstats_DE.csv")
IRF4_all <- read.csv("~/ash_2024/apms/MM1S_3R24M_IRF4_noMito_FP_byMSstats_DE.csv")
IKZF1_all <- read.csv("~/ash_2024/apms/MM1S_B5ran_IKZF_noMito_FP_byMSstats_DE_IKZF1.csv")
IKZF3_all <- read.csv("~/ash_2024/apms/MM1S_B5ran_IKZF_noMito_FP_byMSstats_DE_IKZF3.csv")

#IRF4 and p300 AP-MS are from three independent AP-MS experiments, IKZF1 and IKZF3 from a single experiment
#I use -Log10 adj.pval >2 for IRF4 (LFC > 1.2) and p300 (LFC >1) - Yupeng's suggestions
#use less stringent cutoff for IKZF1 and IKZF3 to account for the fact they are from single experiment, log10 adj.pval >1.3 + LFC >0.5

IKZF1_sig <- IKZF1_all %>% filter(logadjpval > 1.3, log2FC > 0.5)
IKZF3_sig <- IKZF3_all %>% filter(logadjpval > 1.3, log2FC > 0.5)

p300_genes <- p300_sig %>% pull(Prey) %>% unique()
p300_bg_gene <- p300_all %>% pull(Gene) %>% unique()

IRF4_genes <- IRF4_sig %>% pull(Prey) %>% unique()
IRF4_bg_gene <- IRF4_all %>% pull(Gene) %>% unique()

IKZF1_bg_genes <- IKZF1_all %>% pull(Gene) %>% unique()
IKZF3_bg_genes <- IKZF3_all %>% pull(Gene) %>% unique()

background_set <- c(p300_bg_gene,IRF4_bg_gene,IKZF1_bg_genes,IKZF3_bg_genes) %>% unique()

IKZF1_genes <- IKZF1_sig %>% pull(Gene) %>% unique() 
IKZF3_genes <- IKZF3_sig %>% pull(Gene) %>% unique() 

gene_list <- list("IRF4" = IRF4_genes,"IKZF1" = IKZF1_genes,"IKZF3" = IKZF3_genes,"p300" = p300_genes)

gene_nodes <- data.frame(gene = c("IRF4","IKZF1","IKZF3","p300"), sig_inter = c(length(IRF4_genes),length(IKZF1_genes),length(IKZF3_genes),length(p300_genes)), background = c(length(IRF4_bg_gene),length(IKZF1_bg_genes),length(IKZF3_bg_genes),length(p300_bg_gene)))

compute_chi_square <- function(gene_list, background_set) {
  results <- data.frame(
    gene1 = character(),
    gene2 = character(),
    chi_square_p = numeric(),
    stringsAsFactors = FALSE
  )
  gene_pairs <- combn(names(gene_list), 2, simplify = FALSE)
  for (pair in gene_pairs) {
    gene1 <- gene_list[[pair[1]]] %>% c()
    gene2 <- gene_list[[pair[2]]] %>% c()
    group1 <- gene1[gene1 %in% gene2]
    group4 <- background_set[!(background_set %in% gene1) & !(background_set %in% gene2)]
    group2 <- gene1[!(gene1 %in% gene2)]
    group3 <- gene2[!(gene2 %in% gene1)]
    matrix_2x2 <- matrix(c(length(group1), length(group2), 
                           length(group3), length(group4)),
                         nrow = 2, byrow = TRUE,
                         dimnames = list(c("sig both", "gene1 sig"),
                                         c("gene2 non sig", "gene2 sig")))
    chi_square_result <- chisq.test(matrix_2x2)
    results <- rbind(results, data.frame(
      gene1 = pair[1],
      gene2 = pair[2],
      chi_square_p = chi_square_result$p.value,
      shared = length(group1),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

result_table <- compute_chi_square(gene_list, background_set)
write.csv(result_table,"./apms_strength_network_edges.csv")
write.csv(gene_nodes,"./apms_strength_network_nodes.csv")

