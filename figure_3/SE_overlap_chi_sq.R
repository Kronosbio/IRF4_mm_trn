
##### This code is used for the chi-square test checking overlapping SE genes with the top 1000 genes chosen by enhancerPromoter p300 ChIP signal

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(org.Hs.eg.db)
library(IRanges)
library(readxl)
library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_3")

gene_list_boxplots <- read.csv("./data/gene_list_boxplots.csv")

### Using ROSE2 output superenhancer GENE_TO_ENHANCER to table. Essentially this identifies all nearby genes to SE loci.

DMSO_11_SuperEnhancers_GENE_TO_ENHANCER <- read.delim("~/irf4_mm_trn_code/ChIP_data/11_0EWX_0229Kronos_MM1S-No-Treatment-1_H3K27Ac_hs-dm_i20_peaks_SuperEnhancers_GENE_TO_ENHANCER.txt")
DMSO_12_SuperEnhancers_GENE_TO_ENHANCER <- read.delim("~/irf4_mm_trn_code/ChIP_data/12_0EWY_0229Kronos_MM1S-No-Treatment-2_H3K27Ac_hs-dm_i21_peaks_SuperEnhancers_GENE_TO_ENHANCER.txt")

### background gene list
fig3_heatmap_dat_supp <- read.csv("./data/fig3_heatmap_dat_supp.csv")
bg_gene_list<- fig3_heatmap_dat_supp$GENE %>% unique() 

'%!in%' <- function(x,y)!('%in%'(x,y))

genematch <- function(dat){
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL')
  genes_IDS <- mapIds(org.Hs.eg.db, entrez_IDS,  'SYMBOL','ENTREZID')
  
  genes_IDS <- genes_IDS %>% unlist()
  t <- as.data.frame(genes_IDS)
  t$entrez_IDS <- rownames(t)
  
  t2 <- as.data.frame(entrez_IDS)
  t2$GENE  <- rownames(t2)
  key_name <- merge(t,t2)
  
  key_name$entrez_IDS <- NULL
  dat <- merge(key_name,dat)
  
  rm(t,t2,key_name)
  dat$GENE <- NULL
  colnames(dat)[1] <- "GENE"
  return(dat)
}

### attempting to match gene names to entrez ids

DMSO_11_SuperEnhancers_GENE_TO_ENHANCER$GENE <- DMSO_11_SuperEnhancers_GENE_TO_ENHANCER$GENE_NAME
DMSO_11_SuperEnhancers_GENE_TO_ENHANCER <- genematch(DMSO_11_SuperEnhancers_GENE_TO_ENHANCER)

DMSO_12_SuperEnhancers_GENE_TO_ENHANCER$GENE <- DMSO_12_SuperEnhancers_GENE_TO_ENHANCER$GENE_NAME
DMSO_12_SuperEnhancers_GENE_TO_ENHANCER <- genematch(DMSO_12_SuperEnhancers_GENE_TO_ENHANCER)

data_expanded <- DMSO_11_SuperEnhancers_GENE_TO_ENHANCER  %>%
  separate_rows(PROXIMAL_ENHANCERS, sep = ",") %>%
  separate(PROXIMAL_ENHANCERS, into = c("chr", "positions"), sep = ":") %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(across(start:end, as.numeric))

# Define the target range for MYC, custom adding MYC
target_chr <- "chr22"
target_start <- 23182921
target_end <- 23344000

### enhancers near MYC-IGH loci are presentin both samples
data_expanded <- data_expanded %>% filter(chr == target_chr) %>%
  filter(start >= target_start & end <= target_end)

data_expanded <- DMSO_12_SuperEnhancers_GENE_TO_ENHANCER  %>%
  separate_rows(PROXIMAL_ENHANCERS, sep = ",") %>%
  separate(PROXIMAL_ENHANCERS, into = c("chr", "positions"), sep = ":") %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(across(start:end, as.numeric))

# Define the target range for MYC, custom adding MYC
target_chr <- "chr22"
target_start <- 23182921
target_end <- 23344000

### enhancers near MYC-IGH loci are present
data_expanded <- data_expanded %>% filter(chr == target_chr) %>%
  filter(start >= target_start & end <= target_end)

H3K27ac_SE_genes_1 <- DMSO_11_SuperEnhancers_GENE_TO_ENHANCER$GENE
H3K27ac_SE_genes_2 <- DMSO_12_SuperEnhancers_GENE_TO_ENHANCER$GENE

H3K27ac_SE_genes <- intersect(H3K27ac_SE_genes_1,H3K27ac_SE_genes_2)

###custom changing DUSP22 to IRF4
H3K27ac_SE_genes[H3K27ac_SE_genes == "DUSP22"] <- "IRF4"
H3K27ac_SE_genes <- c(H3K27ac_SE_genes,"MYC")

### selecting for genes chosen by enhancerPromoter
H3K27ac_SE_genes <- H3K27ac_SE_genes[H3K27ac_SE_genes %in% fig3_heatmap_dat_supp$GENE] 

top_genes <- gene_list_boxplots$top_1000

bg <- fig3_heatmap_dat_supp$GENE[fig3_heatmap_dat_supp$GENE %!in% c(top_genes,H3K27ac_SE_genes)]

a <- bg %>% length()
d <- length(intersect(top_genes,H3K27ac_SE_genes))

b <- length(top_genes) - d
c <- length(H3K27ac_SE_genes) - d

mat_dat <- matrix(c(a,b,c,d),nrow = 2)

chisq.test(mat_dat)

# > chisq.test(mat_dat)
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  mat_dat
# X-squared = 2325.2, df = 1, p-value < 2.2e-16

