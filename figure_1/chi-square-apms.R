rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(readxl)
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
#   [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1  
# [10] tidyverse_2.0.0 readxl_1.4.3   
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.44.0         fastmatch_1.1-4         gtable_0.3.5            GOSemSim_2.30.0         Biobase_2.64.0          lattice_0.22-6         
# [7] tzdb_0.4.0              vctrs_0.6.5             tools_4.4.0             generics_0.1.3          yulab.utils_0.1.4       stats4_4.4.0           
# [13] parallel_4.4.0          fansi_1.0.6             AnnotationDbi_1.66.0    RSQLite_2.3.6           blob_1.2.4              pkgconfig_2.0.3        
# [19] Matrix_1.7-0            data.table_1.15.4       S4Vectors_0.42.0        lifecycle_1.0.4         GenomeInfoDbData_1.2.12 HDO.db_0.99.1          
# [25] compiler_4.4.0          Biostrings_2.72.0       munsell_0.5.1           fgsea_1.30.0            codetools_0.2-20        DOSE_3.30.1            
# [31] GenomeInfoDb_1.40.0     pillar_1.9.0            crayon_1.5.2            GO.db_3.19.1            BiocParallel_1.38.0     cachem_1.0.8           
# [37] tidyselect_1.2.1        digest_0.6.35           stringi_1.8.4           reshape2_1.4.4          splines_4.4.0           cowplot_1.1.3          
# [43] fastmap_1.1.1           grid_4.4.0              colorspace_2.1-0        cli_3.6.2               magrittr_2.0.3          utf8_1.2.4             
# [49] withr_3.0.0             scales_1.3.0            UCSC.utils_1.0.0        bit64_4.0.5             timechange_0.3.0        XVector_0.44.0         
# [55] httr_1.4.7              bit_4.0.5               qvalue_2.36.0           cellranger_1.1.0        hms_1.1.3               png_0.1-8              
# [61] memoise_2.0.1           IRanges_2.38.0          rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2              
# [67] BiocGenerics_0.50.0     rstudioapi_0.16.0       jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1                fs_1.6.4               
# [73] zlibbioc_1.50.0   

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x,y)!('%in%'(x,y))

###  significant APMS genes

IRF4_MS <- read_excel("./Table S2.xlsx",sheet = 1)
p300_MS <- read_excel("./Table S2.xlsx",sheet = 2)

### reading in all genes, background 

IRF4_MS_bg <- read_excel("./Table S2.xlsx",sheet = 3)
p300_MS_bg <- read_excel("./Table S2.xlsx",sheet = 4)

IRF4_sig_MS <- IRF4_MS$`Gene Symbol` %>% unique()
p300_sig_MS <- p300_MS$`Gene Symbol` %>% unique()

IRF4_MS_bg <- IRF4_MS_bg$`Gene Symbol` %>% unique()
p300_MS_bg <- p300_MS_bg$`Gene Symbol` %>% unique()

IRF4_MS_bg <- IRF4_MS_bg[ IRF4_MS_bg %!in% IRF4_sig_MS]
p300_MS_bg <- p300_MS_bg[ p300_MS_bg %!in% p300_sig_MS]

a <- length(intersect(IRF4_MS_bg,p300_MS_bg))
d <- length(intersect(IRF4_sig_MS,p300_sig_MS))

b <- length(IRF4_sig_MS) - d
c <- length(p300_sig_MS) - d

mat_dat <- matrix(c(a,b,c,d),nrow = 2)

chisq.test(mat_dat)

# > chisq.test(mat_dat)
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  mat_dat
# X-squared = 20.963, df = 1, p-value = 4.683e-06
