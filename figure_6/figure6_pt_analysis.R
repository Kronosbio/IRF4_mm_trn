###########

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(patchwork)
library(tidyverse)
library(preprocessCore)
library(ggpubr)

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.6.0          preprocessCore_1.61.0 lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
# [7] purrr_1.0.2           readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0      
# [13] patchwork_1.2.0      
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1          
# [6] digest_0.6.35           timechange_0.3.0        lifecycle_1.0.4         KEGGREST_1.44.0         RSQLite_2.3.6          
# [11] magrittr_2.0.3          compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4             
# [16] data.table_1.15.4       ggsignif_0.6.4          bit_4.0.5               plyr_1.8.9              abind_1.4-5            
# [21] BiocParallel_1.38.0     withr_3.0.0             BiocGenerics_0.50.0     grid_4.4.0              stats4_4.4.0           
# [26] fansi_1.0.6             GOSemSim_2.30.0         colorspace_2.1-0        GO.db_3.19.1            scales_1.3.0           
# [31] cli_3.6.2               crayon_1.5.2            generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7             
# [36] reshape2_1.4.4          tzdb_0.4.0              DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8           
# [41] DOSE_3.30.1             zlibbioc_1.50.0         splines_4.4.0           parallel_4.4.0          AnnotationDbi_1.66.0   
# [46] XVector_0.44.0          yulab.utils_0.1.4       vctrs_0.6.5             Matrix_1.7-0            carData_3.0-5          
# [51] jsonlite_1.8.8          car_3.1-2               IRanges_2.38.0          hms_1.1.3               S4Vectors_0.42.0       
# [56] bit64_4.0.5             rstatix_0.7.2           glue_1.7.0              codetools_0.2-20        cowplot_1.1.3          
# [61] stringi_1.8.4           gtable_0.3.5            GenomeInfoDb_1.40.0     UCSC.utils_1.0.0        munsell_0.5.1          
# [66] pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12 R6_2.5.1                lattice_0.22-6         
# [71] Biobase_2.64.0          backports_1.4.1         png_0.1-8               memoise_2.0.1           broom_1.0.5            
# [76] Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3        


'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
setwd("~/irf4_mm_trn_code/figure_6")
dat <- read.csv("./pt_KB528_TPM.csv")

###data contains a lot of zeroes, so quantile normalizing data to have fixed distributions

q_norm <-
  preprocessCore::normalize.quantiles(as.matrix(dat[-c(1, 2)]))
dat[-c(1, 2)] <- q_norm

colnames(dat)[-c(1, 2)]

###cleaning up column names

colnames(dat) <-
  stringr::str_replace(colnames(dat), "_L002_R1_001.fastq.gz", "")

colnames(dat)[-c(1, 2)] <-
  sub("_S[0-9]+$", "", colnames(dat)[-c(1, 2)])

print(dat[1:5, ])

dat$FEATURE_ID <- NULL

annotation_key <- read.csv("./table_supp_df_pt_annotations.csv")

annotation_key$group <-
  paste(annotation_key$Treatment, annotation_key$Time, sep = "_")
annotation_key$label <- annotation_key$Sample

dat$FEATURE_ID <- NULL
dat <- dat %>% gather(label, Exp, -GENE_SYMBOL)

###matching annotation key to data

dat <- merge(dat, annotation_key)
genes <- c("IKZF3", "IKZF1", "IRF4", "MYC", "CRBN", "UBC", "EP300")

###plotting bar plot for each gene of interest, p300 and UBC not used in manuscript

for (gene in genes) {
  t <- dat %>% filter(GENE_SYMBOL == gene, Time == 12)
  
  ### normalizing gene expression data to control sample means

  t$normalize_1 <- t$Exp
  ctrl_mean <-
    t %>% filter(Treatment == "ctrl") %>% pull(Exp) %>% mean()
  t$normalize_1 <- t$normalize_1 / ctrl_mean
  means <-
    t %>% group_by(Treatment) %>% summarize(normalize_1 = mean(normalize_1))
  se_error <-
    t %>% group_by(Treatment) %>% summarize(normalize_1_sd = sd(normalize_1))
  means <- merge(means, se_error)
  filename =  paste("~/paper/fig5/", gene, "_DF_pt_12_hours.pdf", sep = "")
  means$ymin = means$normalize_1 - means$normalize_1_sd
  
  ##setting ymin to 0 for error bar purposes.
  means$ymin[means$ymin < 0] <- 0
  
  p <- ggplot(t, aes(y = normalize_1, x = Treatment)) +
    geom_bar(
      data = means,
      stat = 'identity',
      width = 0.5,
      aes(fill = Treatment)
    ) +
    geom_errorbar(
      data = means,
      aes(ymin = ymin, ymax = normalize_1 + normalize_1_sd),
      width = .2,
      size = 1
    ) +
    geom_line(data = t, aes(group = rep), color = "gray") +
    geom_point(size = 5, alpha = 0.4) +
    xlab("Treatment") +
    ylab("Quantile Norm to CTRL") +
    theme(legend.position = "none") +
    theme_classic(base_size = 15) + # ylim(c(0,1.75)) +
    ggtitle(paste(gene, " 12 Hour")) +
    ylim(c(0, 1.75)) +
    stat_compare_means(method = "t.test", paired = T) +
    scale_fill_manual(values = c("#999999", "dodgerblue4")) + theme(legend.position = "none")
  pdf(filename, height = 5, width = 3)
  show(p)
  dev.off()
}
