rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(tidyverse)
library(ggpubr)
library(cowplot)

sessionInfo()

setwd("~/irf4_mm_trn_code/figure_s1/")

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
#   [1] cowplot_1.1.3   ggpubr_0.6.0    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1    
# [10] tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1           digest_0.6.35          
# [7] timechange_0.3.0        lifecycle_1.0.4         KEGGREST_1.44.0         RSQLite_2.3.6           magrittr_2.0.3          compiler_4.4.0         
# [13] rlang_1.1.3             tools_4.4.0             utf8_1.2.4              data.table_1.15.4       ggsignif_0.6.4          bit_4.0.5              
# [19] plyr_1.8.9              abind_1.4-5             BiocParallel_1.38.0     withr_3.0.0             BiocGenerics_0.50.0     grid_4.4.0             
# [25] stats4_4.4.0            fansi_1.0.6             GOSemSim_2.30.0         colorspace_2.1-0        GO.db_3.19.1            scales_1.3.0           
# [31] cli_3.6.2               crayon_1.5.2            generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7              reshape2_1.4.4         
# [37] tzdb_0.4.0              DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1             zlibbioc_1.50.0        
# [43] splines_4.4.0           parallel_4.4.0          AnnotationDbi_1.66.0    XVector_0.44.0          yulab.utils_0.1.4       vctrs_0.6.5            
# [49] Matrix_1.7-0            carData_3.0-5           jsonlite_1.8.8          car_3.1-2               IRanges_2.38.0          hms_1.1.3              
# [55] S4Vectors_0.42.0        bit64_4.0.5             ggrepel_0.9.5           rstatix_0.7.2           glue_1.7.0              codetools_0.2-20       
# [61] stringi_1.8.4           gtable_0.3.5            GenomeInfoDb_1.40.0     UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0           
# [67] fgsea_1.30.0            GenomeInfoDbData_1.2.12 R6_2.5.1                lattice_0.22-6          Biobase_2.64.0          png_0.1-8              
# [73] backports_1.4.1         memoise_2.0.1           broom_1.0.5             Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4               
# [79] pkgconfig_2.0.3   

'%!in%' <- function(x,y)!('%in%'(x,y))

metadata <-  read.csv("../public_data/Model.csv")
Chronos <-  read.csv("../public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")

Chronos <- Chronos %>% gather(GENE,Dep,-X)

Chronos_chisq <- Chronos 

### binarizing dependency (0 = non dependent, -1 = dependent)

Chronos_chisq$Dep_bin <- 0
Chronos_chisq$Dep_bin[Chronos_chisq$Dep < -1] <- 1

### selecting 3 genes of interest

Chronos <- Chronos %>% filter(GENE %in% c("EP300","IRF4","CREBBP")) %>% spread(GENE,Dep) 

Chronos$Type <- "Other"
Chronos$Type[Chronos$X %in% metadata$ModelID] <- "MM"

Chronos_chisq$Dep <- NULL

Chronos_chisq_test <- Chronos_chisq %>% filter(GENE %in% c("EP300","IRF4")) %>% spread(GENE,Dep_bin) 

### calculating chi square of co-dependence

chi_sq_res <- chisq.test(table(Chronos_chisq_test$EP300,Chronos_chisq_test$IRF4))

### p300 vs IRF4

a <- Chronos %>% ggplot(aes(x = EP300,y = IRF4, color = Type)) + geom_point(size = 5,alpha = 0.7) +
  theme_classic(base_size = 20) + scale_color_manual(values=c( "#4c00a4","#999999")) +
  geom_hline(yintercept = -1,linetype = "dashed") +
  geom_vline(xintercept = -1,linetype = "dashed") +
  xlim(c(-2.5,1.25)) +
  ylim(c(-2.5,1.25))  + 
  xlab("p300 Chronos DepMap22Q4") +
  ylab("IRF4 Chronos DepMap22Q4") +
  theme(legend.position =  "none") +
  geom_text( aes( x=-2, y=1,
                  label=paste("Chi-Square p = ",format(chi_sq_res$p.value,digits = 3))), 
            color="black", 
            size=6 )

a

Chronos_chisq_test <- Chronos_chisq %>% filter(GENE %in% c("CREBBP","IRF4")) %>% spread(GENE,Dep_bin) 

chi_sq_res <- chisq.test(table(Chronos_chisq_test$CREBBP,Chronos_chisq_test$IRF4))

### CREBBP vs IRF4

b <- Chronos %>% ggplot(aes(x = CREBBP,y = IRF4, color = Type)) + geom_point(size = 5,alpha = 0.7) +
  theme_classic(base_size = 20) + scale_color_manual(values=c( "#4c00a4","#999999")) +
  geom_hline(yintercept = -1,linetype = "dashed") +
  geom_vline(xintercept = -1,linetype = "dashed") +
  xlim(c(-2.5,1.25))+ ylim(c(-2.5,1.25))  + 
  xlab("CBP Chronos DepMap22Q4") +
  ylab("IRF4 Chronos DepMap22Q4") +
  theme(legend.position =  "none") +
  geom_text( aes( x=-2, y=1,
                  label=paste("Chi-Square p = ",format(chi_sq_res$p.value,digits = 3))), 
             color="black", 
             size=6 )


### running these seperately to calculate wilcox test with stat_compare means (ggpubr)

Chronos %>% ggplot(aes(y = CREBBP, x = as.factor(Type)),outlier.alpha = 0) + geom_boxplot() + stat_compare_means()
#Wilcox Test p = 1E-4
dev.off()

cbp_fig <- Chronos %>% ggplot() +
  xlim(c(-2.5,1.25)) + 
  geom_boxplot(data = Chronos,aes(CREBBP, y = as.factor(Type),fill = Type,color = Type),outlier.alpha = 0,
               width = 0.3,alpha = 0.9,
               color =c("lightgray","black"),
               fill =c("#4c00a4","white")) +
  theme_classic(base_size = 20) + 
  xlab("") + scale_color_manual(values=c("#4c00a4","white")) + 
  geom_vline(xintercept = -1,color  = "red", linetype = "dashed") + 
  ylab("CBP Chronos Score") + 
  annotate("text", x=1.5, y=1, label= "Wilcox Test p = 1E-4")

### running these seperately to calculate wilcox test with stat_compare means (ggpubr)
Chronos %>% ggplot(aes(y = EP300, x = as.factor(Type)),outlier.alpha = 0) + geom_boxplot() + stat_compare_means()
#Wilcox Test p = 3.9E-10
dev.off()

p300_fig <-  Chronos %>% ggplot() +
  xlim(c(-2.5,1.25)) + 
  geom_boxplot(data = Chronos,aes(EP300, y = as.factor(Type),fill = Type,color = Type),outlier.alpha = 0,
               width = 0.3,alpha = 0.9,
               color =c("lightgray","black"),
               fill =c("#4c00a4","white")) +
  theme_classic(base_size = 20) + 
  xlab("") + scale_color_manual(values=c("#4c00a4","white")) + 
  geom_vline(xintercept = -1,color  = "red", linetype = "dashed") + 
  ylab("p300 Chronos Score") + 
  annotate("text", x=1.5, y=0, label= "Wilcox Test p = 3.9E-10")

p300_fig

Chronos %>% ggplot(aes(IRF4, x = as.factor(Type)),outlier.alpha = 0) + geom_boxplot() + stat_compare_means()
dev.off()
#Wilcox Test p = 5.4E-14

IRF4_fig <- Chronos %>% ggplot() +
  ylim(c(-2.5,1.25)) + 
  geom_boxplot(data = Chronos,aes(IRF4, x = as.factor(Type),fill = Type,color = Type),outlier.alpha = 0,
               width = 0.3,alpha = 0.9,
               color =c("lightgray","black"),
               fill =c("#4c00a4","white")) +
  theme_classic(base_size = 20) + 
  xlab("") + scale_color_manual(values=c("#4c00a4","white")) + 
  geom_hline(yintercept = -1,color  = "red", linetype = "dashed") + 
  ylab("IRF4 Chronos Score") + 
  annotate("text", x=1.5, y=1, label= "Wilcox Test p = 5.4E-14")

### an update causes warnings now about missing values. The non plotted points are geom_text, and were the wilcox tests. Stats added in illustrator
### This warning came up after an update that caused geom_text to break with the factor y-axis.

## This output was cleaned in illustrator, but data and proportions are the same

pdf("./fig_s1_IRF4_p300_Chronos_score.pdf", height = 10,width = 15)
plot_grid(p300_fig,cbp_fig,IRF4_fig,
          a,b,IRF4_fig,ncol = 3)
dev.off()


### Not used, but calculating correlation of dependency score between IRF4 and p300/CBP in all, MM, and non-MM cell lines

cor.test(Chronos$CREBBP,Chronos$IRF4)
# Pearson's product-moment correlation
# 
# data:  Chronos$CREBBP and Chronos$IRF4
# t = 3.4573, df = 1076, p-value = 0.000567
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.04539334 0.16350049
# sample estimates:
#       cor 
# 0.1048165 

cor.test(Chronos$EP300,Chronos$IRF4)
# Pearson's product-moment correlation
# 
# data:  Chronos$EP300 and Chronos$IRF4
# t = 12.803, df = 1076, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3106260 0.4143031
# sample estimates:
#       cor 
# 0.3635899 

t <- Chronos %>% filter(Type == "MM")
cor.test(t$CREBBP,t$IRF4)
# Pearson's product-moment correlation
# 
# data:  t$CREBBP and t$IRF4
# t = -0.36103, df = 18, p-value = 0.7223
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5082399  0.3716778
# sample estimates:
#         cor 
# -0.08478861 

cor.test(t$EP300,t$IRF4)
# Pearson's product-moment correlation
# 
# data:  t$EP300 and t$IRF4
# t = 1.311, df = 18, p-value = 0.2063
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1694264  0.6525036
# sample estimates:
#       cor 
# 0.2952291 


t <- Chronos %>% filter(Type == "Other")
cor.test(t$CREBBP,t$IRF4)

# Pearson's product-moment correlation
# 
# data:  t$CREBBP and t$IRF4
# t = 0.27933, df = 1056, p-value = 0.78
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.05170034  0.06882912
# sample estimates:
#        cor 
# 0.00859561 

cor.test(t$EP300,t$IRF4)

# Pearson's product-moment correlation
# 
# data:  t$EP300 and t$IRF4
# t = 9.2617, df = 1056, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2174173 0.3289302
# sample estimates:
#       cor 
# 0.2740948 

