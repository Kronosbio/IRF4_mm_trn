
### This code was used to generate the multi-omic plot for IKZF3 in MM1S cells treated with KB528. 
### The plot shows the normalized omic levels of RNA-Seq, Acetylomics, and Protein levels of IKZF3 at different time points. 

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

library(patchwork)
library(tidyverse)

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
#   [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
# [9] ggplot2_3.5.1   tidyverse_2.0.0 patchwork_1.2.0
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.44.0         fastmatch_1.1-4         gtable_0.3.5            GOSemSim_2.30.0         Biobase_2.64.0         
# [6] lattice_0.22-6          tzdb_0.4.0              vctrs_0.6.5             tools_4.4.0             generics_0.1.3         
# [11] yulab.utils_0.1.4       stats4_4.4.0            parallel_4.4.0          fansi_1.0.6             AnnotationDbi_1.66.0   
# [16] RSQLite_2.3.6           blob_1.2.4              pkgconfig_2.0.3         Matrix_1.7-0            data.table_1.15.4      
# [21] S4Vectors_0.42.0        lifecycle_1.0.4         GenomeInfoDbData_1.2.12 HDO.db_0.99.1           compiler_4.4.0         
# [26] Biostrings_2.72.0       munsell_0.5.1           fgsea_1.30.0            codetools_0.2-20        DOSE_3.30.1            
# [31] GenomeInfoDb_1.40.0     pillar_1.9.0            crayon_1.5.2            GO.db_3.19.1            BiocParallel_1.38.0    
# [36] cachem_1.0.8            tidyselect_1.2.1        digest_0.6.35           stringi_1.8.4           reshape2_1.4.4         
# [41] splines_4.4.0           cowplot_1.1.3           fastmap_1.1.1           grid_4.4.0              colorspace_2.1-0       
# [46] cli_3.6.2               magrittr_2.0.3          utf8_1.2.4              withr_3.0.0             scales_1.3.0           
# [51] UCSC.utils_1.0.0        bit64_4.0.5             timechange_0.3.0        XVector_0.44.0          httr_1.4.7             
# [56] bit_4.0.5               qvalue_2.36.0           hms_1.1.3               png_0.1-8               memoise_2.0.1          
# [61] IRanges_2.38.0          rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2              
# [66] BiocGenerics_0.50.0     rstudioapi_0.16.0       jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1               
# [71] fs_1.6.4                zlibbioc_1.50.0   

setwd("~/irf4_mm_trn_code/figure_3")

dat <- read.delim("~/irf4_mm_trn_code/RNA_data/MM1S_KB528_timeseries_TPM.csv")

colnames(dat) <- stringr::str_replace(colnames(dat), "_L003", "")
colnames(dat) <- stringr::str_replace(colnames(dat), "X", "S")
colnames(dat) <- stringr::str_replace(colnames(dat), "Copy.of.", "S")

annotation_MM1S <- read.csv("~/irf4_mm_trn_code/RNA_data/annotation_MM1S_timeseries.csv")

dat$FEATURE_ID <- NULL

dat <- dat %>% gather(label, TPM, -GENE_SYMBOL)

head(annotation_MM1S)
head(dat)

dat <- merge(dat, annotation_MM1S)

plot_dat <- dat %>% filter(GENE_SYMBOL == "IKZF3")

t <- plot_dat %>% group_by(condition) %>% summarize(med_tpm = median(TPM))

norm_val <- t %>% filter(condition == "DMSO") %>% pull(med_tpm)

plot_dat <- merge(plot_dat, t)

plot_dat$TPM_norm <- plot_dat$TPM / norm_val

t <- plot_dat %>% group_by(condition) %>% summarize(sd_tpm = sd(TPM_norm))

plot_dat <- merge(plot_dat, t)

plot_dat <- plot_dat %>% dplyr::select(condition, TPM_norm, sd_tpm) %>% unique()

colnames(plot_dat) <- c("condition", "norm", "norm_sd")

plot_dat$time <- c(0, 0, 0, 1, 1, 1, 24, 24, 24, 3, 3, 3, 6, 6, 6)
plot_dat$type <- "RNA-Seq"
plot_dat$condition <- NULL

IKZF3_omic_comparison <- read.csv("./data/IKZF3_omic_comparison.csv")
ace_dat <- IKZF3_omic_comparison %>% dplyr::select(Acetylomics, IKZF3.K369) %>% filter(!is.na(IKZF3.K369))

t <- ace_dat %>% group_by(Acetylomics) %>% summarize(med_tpm = median(IKZF3.K369))

norm_val <- t %>% filter(Acetylomics == "DMSO") %>% pull(med_tpm)

ace_dat <- merge(t, ace_dat)
ace_dat$Acetylomics_norm <- ace_dat$IKZF3.K369 / norm_val

t <- ace_dat %>% group_by(Acetylomics) %>% summarize(sd_tpm = sd(Acetylomics_norm))
ace_dat <- merge(t, ace_dat)

ace_dat$time <- c(0, 0, 0, 1, 1, 1, 2, 2, 2)
ace_dat$type <- "Acetylomics"

ace_dat <- ace_dat %>% dplyr::select(Acetylomics_norm, sd_tpm, time, type)

colnames(ace_dat) <- colnames(plot_dat)

plot_dat <- rbind(ace_dat, plot_dat)

pro_dat <- IKZF3_omic_comparison %>% dplyr::select(Protein.levels, IKZF3) %>% filter(!is.na(IKZF3))

t <- pro_dat  %>% group_by(Protein.levels) %>% summarize(sd_tpm = sd(IKZF3))
pro_dat <- merge(t, pro_dat)

pro_dat$type <- "Protein"
pro_dat$time <- c(0, 0, 0, 24, 24, 24, 24, 24, 6, 6, 6, 6)
pro_dat <- pro_dat %>% dplyr::select(IKZF3, sd_tpm, time, type)

colnames(pro_dat) <- colnames(plot_dat)

plot_dat <- rbind(pro_dat, plot_dat)

plot_dat_t <- plot_dat %>% group_by(time, type) %>% summarize(med_norm = median(norm))

plot_dat$norm <- NULL

plot_dat <- merge(plot_dat_t, plot_dat) %>% unique()


p <- ggplot(plot_dat, aes(
  y = med_norm,
  x = time,
  color = type,
  group = type
)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_line(size = 1, alpha = 0.4) +
  geom_errorbar(
    data = plot_dat,
    aes(ymin = med_norm - norm_sd, ymax = med_norm + norm_sd),
    width = 0.3,
    size = 1
  ) +
  theme_classic(base_size = 18) + xlab("Time") + ylab("Normalized Omic Level") + ggtitle("IKZF3") +
  scale_color_manual(values = c("#301C09", "#8B5E3C", "#C49A6C")) + ylim(c(0, 1.65))

p

pdf("./fig3_IKZF3_KB528_omics_combine.pdf",
    height = 4,
    width = 5)
p
dev.off()
