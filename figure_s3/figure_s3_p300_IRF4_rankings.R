rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)
library(patchwork)

# > sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.2.0      lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
# [7] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0      org.Hs.eg.db_3.19.1 
# [13] AnnotationDbi_1.66.0 IRanges_2.38.0       S4Vectors_0.42.0     Biobase_2.64.0       BiocGenerics_0.50.0 
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1          
# [6] digest_0.6.35           timechange_0.3.0        lifecycle_1.0.4         KEGGREST_1.44.0         RSQLite_2.3.6          
# [11] magrittr_2.0.3          compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4             
# [16] data.table_1.15.4       ggsignif_0.6.4          bit_4.0.5               plyr_1.8.9              abind_1.4-5            
# [21] BiocParallel_1.38.0     withr_3.0.0             grid_4.4.0              fansi_1.0.6             GOSemSim_2.30.0        
# [26] ggpubr_0.6.0            colorspace_2.1-0        GO.db_3.19.1            scales_1.3.0            cli_3.6.2              
# [31] crayon_1.5.2            generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7              reshape2_1.4.4         
# [36] tzdb_0.4.0              DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1            
# [41] zlibbioc_1.50.0         splines_4.4.0           parallel_4.4.0          XVector_0.44.0          yulab.utils_0.1.4      
# [46] vctrs_0.6.5             Matrix_1.7-0            carData_3.0-5           jsonlite_1.8.8          car_3.1-2              
# [51] hms_1.1.3               rstatix_0.7.2           bit64_4.0.5             ggrepel_0.9.5           glue_1.7.0             
# [56] codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.5            GenomeInfoDb_1.40.0    
# [61] UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12
# [66] R6_2.5.1                lattice_0.22-6          backports_1.4.1         png_0.1-8               broom_1.0.5            
# [71] memoise_2.0.1           Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3  

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

setwd("~/irf4_mm_trn_code/figure_s3")

##############

### reading in RNA data to remove genes with DMSO < 1

dat <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_3_MM1S_exprs_matrix.txt")

dat <- dat %>% filter(DMSO_MM1S > 1)

### reading in enhancer promoter data, this follows the same pipeline for the omic heatmap in figure 3
### average the signal for replicates, rank the genes based on the signal, and then plot the rankings against IRF4 and p300 signal

enhancerPromoter_all <- read.csv("../ChIP_data/enhancerPromoter_DMSO_with_MYC.csv") %>% unique()

enhancerPromoter_all$comb_sig <- enhancerPromoter_all$DISTAL_signal + enhancerPromoter_all$TSS_signal

ep_dat <- enhancerPromoter_all %>% dplyr::select(GENE_ID,file,comb_sig) %>% spread(file,comb_sig)

ep_dat <- ep_dat[, c(1:2, (3+1):11, 3,12:ncol(ep_dat))]

ep_dat[is.na(ep_dat)] <- 0.0

mean_cols <- sapply(seq(2, ncol(ep_dat), by = 2), function(i) rowMeans(ep_dat[, i:(i+1)])) %>% as.data.frame()

colnames(mean_cols) <- colnames(ep_dat)[seq(3, ncol(ep_dat), by = 2)]
mean_cols$GENE <- ep_dat$GENE_ID

ep_dat <- mean_cols
ep_dat <- ep_dat %>% dplyr::select(GENE, everything())

ep_dat$rank_ep300 <- rank(-ep_dat$p300H3K27ac_34_DMSO_compare) 
ep_dat$rank_IRF4 <- rank(-ep_dat$IRF4H3K27ac_32_DMSO_compare) 

ep_dat$tot_rank <- ep_dat$rank_ep300 + ep_dat$rank_IRF4

ep_dat <- ep_dat %>% filter(GENE != "IRF4")
ep_dat$GENE[ep_dat$GENE == "DUSP22"] <- "IRF4"

ep_dat_p300_IRF4 <- ep_dat

rm(enhancerPromoter_all,mean_cols)

#######

### matching gene names

dat$GENE <- rownames(dat)
dat <- genematch(dat)
ep_dat <- genematch(ep_dat_p300_IRF4)

t1 <- ep_dat %>% dplyr::select(GENE,p300H3K27ac_34_DMSO_compare,IRF4H3K27ac_32_DMSO_compare,tot_rank,rank_ep300,rank_IRF4)

### tot_rank isn't re ranked originally because it is spearman correlations, it is re-ranked for the x-axis here. 

t1$re_tot_rank <- rank(t1$tot_rank)

### highlighting mm TRN tfs
picked_genes <- c("ATF4","PRDM10","IRF4","PRDM1","MEF2C","MYC","POU2F1","TCF3","ZNF217")

t1$top <- "No"
t1$top[t1$re_tot_rank < 1001] <- "Yes"
t1$label <- ""
t1$label[t1$GENE %in% picked_genes] <- t1$GENE[t1$GENE %in% picked_genes] 

p1 <- t1 %>% ggplot(aes(x = re_tot_rank, y = p300H3K27ac_34_DMSO_compare, color = top,label = label)) +
  scale_color_manual(values=c("#999999", "cyan")) +
  geom_point(size = 3, alpha = 0.5) + ylab("Total p300 AUC signal") + 
  xlab("Average ranking of IRF4 and p300 bound genes") + theme_classic(base_size = 15) +
  theme(legend.position = "none") +  ggrepel::geom_text_repel(color = "black",max.overlaps = Inf)

p2 <- t1 %>% ggplot(aes(x = re_tot_rank, y = IRF4H3K27ac_32_DMSO_compare, color = top,label = label)) +
  scale_color_manual(values=c("#999999", "magenta")) +
  geom_point(size = 3, alpha = 0.5) + ylab("Total IRF4AUC signal") + 
  xlab("Average ranking of IRF4 and p300 bound genes") + theme_classic(base_size = 15) +
  theme(legend.position = "none")+ ggrepel::geom_text_repel(color = "black",max.overlaps = Inf)


pdf("fig_S3_IRF4_p300_rankings.pdf", width = 10, height = 5)
p2+p1
dev.off()
