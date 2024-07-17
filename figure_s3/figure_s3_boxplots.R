### this code is used for generating the boxplots used in figure 3 and S3
### essentially we are identifying the top and bottom 1,000 bound eligible genes and examining the logFC of these genes in RNA-seq, H3K18ac, and H3K27ac

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(org.Hs.eg.db)
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2          readr_2.1.5         
# [7] tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0      org.Hs.eg.db_3.19.1  AnnotationDbi_1.66.0
# [13] IRanges_2.38.0       S4Vectors_0.42.0     Biobase_2.64.0       BiocGenerics_0.50.0 
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1          
# [6] digest_0.6.35           timechange_0.3.0        lifecycle_1.0.4         KEGGREST_1.44.0         RSQLite_2.3.6          
# [11] magrittr_2.0.3          compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4             
# [16] data.table_1.15.4       ggsignif_0.6.4          bit_4.0.5               plyr_1.8.9              abind_1.4-5            
# [21] BiocParallel_1.38.0     withr_3.0.0             grid_4.4.0              fansi_1.0.6             GOSemSim_2.30.0        
# [26] ggpubr_0.6.0            colorspace_2.1-0        GO.db_3.19.1            scales_1.3.0            cli_3.6.2              
# [31] crayon_1.5.2            generics_0.1.3          rstudioapi_0.16.0       tzdb_0.4.0              httr_1.4.7             
# [36] reshape2_1.4.4          DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1            
# [41] zlibbioc_1.50.0         splines_4.4.0           parallel_4.4.0          XVector_0.44.0          yulab.utils_0.1.4      
# [46] vctrs_0.6.5             Matrix_1.7-0            jsonlite_1.8.8          carData_3.0-5           car_3.1-2              
# [51] hms_1.1.3               bit64_4.0.5             ggrepel_0.9.5           rstatix_0.7.2           glue_1.7.0             
# [56] codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.5            GenomeInfoDb_1.40.0    
# [61] UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12
# [66] R6_2.5.1                lattice_0.22-6          png_0.1-8               backports_1.4.1         memoise_2.0.1          
# [71] broom_1.0.5             Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3  

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


setwd("~/irf4_mm_trn_code/figure_s3/")

##########
### reading in gene lists from figure 3

gene_list_boxplots <- read.csv("../figure_3/data/gene_list_boxplots.csv")
### reading in RNA-seq

dat1 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_1_MM1S_exprs_matrix.txt")
dat3 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_3_MM1S_exprs_matrix.txt")
dat6 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_6_MM1S_exprs_matrix.txt")
dat24 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_24_MM1S_exprs_matrix.txt")

dat1$GENE <- rownames(dat1)
dat3$GENE <- rownames(dat3)
dat6$GENE <- rownames(dat6)
dat24$GENE <- rownames(dat24)

colnames(dat1)[3:4] <- c("LogFC_1Hr","p_1Hr")
colnames(dat3)[3:4] <- c("LogFC_3Hr","p_3Hr")
colnames(dat6)[3:4] <- c("LogFC_6Hr","p_6Hr")
colnames(dat24)[3:4] <- c("LogFC_24Hr","p_24Hr")

dat <- merge(dat1,dat3)
dat <- merge(dat,dat6)
dat <- merge(dat,dat24)

rm(dat1,dat3,dat6,dat24)

##########
### reading in spike in controls

MM_p300_internal_QC_report <- read.csv("../ChIP_data/MM_p300_internal_QC_report.csv") 
MM_p300_internal_QC_report$Scale_factor <- MM_p300_internal_QC_report$Total_Reads_dm3/1000000
MM_p300_internal_QC_report$Scale_factor <- 1/MM_p300_internal_QC_report$Scale_factor

MM_p300_internal_QC_report$Scale_factor[21:34] <- 1

MM_p300_internal_QC_report$X <- NULL
MM_p300_internal_QC_report$X.1 <- NULL
MM_p300_internal_QC_report$X.2 <- NULL

numbers <- gsub("_.*", "", MM_p300_internal_QC_report$Sample) %>% as.numeric()

MM_p300_internal_QC_report$sample_num <- numbers

##########
### reading in ChIP Enhancer promoter data for each metric
enhancerPromoter_DMSO <- read.csv("../ChIP_data/enhancerPromoter_DMSO_with_MYC.csv") %>% unique()

enhancerPromoter_DMSO$comb_sig <- enhancerPromoter_DMSO$DISTAL_signal + enhancerPromoter_DMSO$TSS_signal

enhancerPromoter_DMSO$sample_num <- enhancerPromoter_DMSO$file

enhancerPromoter_DMSO$sample_num <- gsub("_DMSO_compare","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k27ac_","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k18ac_","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k27me3_","",enhancerPromoter_DMSO$sample_num)

enhancerPromoter_DMSO$sample_num <- as.numeric(enhancerPromoter_DMSO$sample_num)

MM_p300_internal_QC_report <- MM_p300_internal_QC_report %>% dplyr::select(sample_num,Scale_factor)

enhancerPromoter_DMSO <- merge(enhancerPromoter_DMSO,MM_p300_internal_QC_report)

enhancerPromoter_DMSO$comb_sig <- enhancerPromoter_DMSO$comb_sig*enhancerPromoter_DMSO$Scale_factor

ep_dat <- enhancerPromoter_DMSO %>% dplyr::select(GENE_ID,file,comb_sig) %>% spread(file,comb_sig)

ep_dat <- ep_dat[, c(1:2, (3+1):11, 3,12:ncol(ep_dat))]

ep_dat[is.na(ep_dat)] <- 0.0

### averaging chip auc signal per replicate.

mean_cols <- sapply(seq(2, ncol(ep_dat), by = 2), function(i) rowMeans(ep_dat[, i:(i+1)])) %>% as.data.frame()

colnames(mean_cols) <- colnames(ep_dat)[seq(3, ncol(ep_dat), by = 2)]
mean_cols$GENE <- ep_dat$GENE_ID

ep_dat <- mean_cols
ep_dat_DMSO <- ep_dat %>% dplyr::select(GENE, everything())

##########
### calculating logFC for each time point, data with NAs are removed from boxplots, as some metrics do not have data for all time points

ep_dat_DMSO$H3K18ac_1Hr <- log(ep_dat_DMSO$H3k18ac_4_DMSO_compare/ep_dat_DMSO$H3k18ac_2_DMSO_compare)
ep_dat_DMSO$H3K18ac_3Hr <- log(ep_dat_DMSO$H3k18ac_6_DMSO_compare/ep_dat_DMSO$H3k18ac_2_DMSO_compare)
ep_dat_DMSO$H3K18ac_6Hr <- log(ep_dat_DMSO$H3k18ac_8_DMSO_compare/ep_dat_DMSO$H3k18ac_2_DMSO_compare)
ep_dat_DMSO$H3K18ac_24Hr <- log(ep_dat_DMSO$H3k18ac_10_DMSO_compare/ep_dat_DMSO$H3k18ac_2_DMSO_compare)

ep_dat_DMSO$H3K27ac_1Hr <- log(ep_dat_DMSO$H3k27ac_14_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_3Hr <- log(ep_dat_DMSO$H3k27ac_16_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_6Hr <- log(ep_dat_DMSO$H3k27ac_18_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_24Hr <- log(ep_dat_DMSO$H3k27ac_20_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)

ep_dat_DMSO$H3K27me3_1Hr <- log(ep_dat_DMSO$H3k27me3_24_DMSO_compare/ep_dat_DMSO$H3k27me3_22_DMSO_compare)
ep_dat_DMSO$H3K27me3_3Hr <- log(ep_dat_DMSO$H3k27me3_26_DMSO_compare/ep_dat_DMSO$H3k27me3_22_DMSO_compare)
ep_dat_DMSO$H3K27me3_6Hr <- log(ep_dat_DMSO$H3k27me3_28_DMSO_compare/ep_dat_DMSO$H3k27me3_22_DMSO_compare)
ep_dat_DMSO$H3K27me3_24Hr <- log(ep_dat_DMSO$H3k27me3_30_DMSO_compare/ep_dat_DMSO$H3k27me3_22_DMSO_compare)


ep_dat_DMSO <- ep_dat_DMSO %>% dplyr::select(GENE,H3K27ac_1Hr,H3K27ac_3Hr,H3K27ac_6Hr,H3K27ac_24Hr,
                                             H3K18ac_1Hr,H3K18ac_3Hr,H3K18ac_6Hr,H3K18ac_24Hr,
                                             H3K27me3_1Hr,H3K27me3_3Hr,H3K27me3_6Hr,H3K27me3_24Hr)

rm(enhancerPromoter_DMSO,mean_cols,ep_dat,MM_p300_internal_QC_report)

########

dat <- genematch(dat)
ep_dat_DMSO <- genematch(ep_dat_DMSO)

dat$type <- "Other"
dat$type[dat$GENE %in% gene_list_boxplots$top_1000] <- "top1000"
dat$type[dat$GENE %in% gene_list_boxplots$bot_1000] <- "bot1000"

ep_dat_DMSO$type <- "Other"
ep_dat_DMSO$type[ep_dat_DMSO$GENE %in% gene_list_boxplots$top_1000] <- "top1000"
ep_dat_DMSO$type[ep_dat_DMSO$GENE %in% gene_list_boxplots$bot_1000] <- "bot1000"

boxplot_dat <- ep_dat_DMSO %>% filter(type != "Other") %>% dplyr::select(type,H3K27ac_1Hr,H3K27ac_3Hr,H3K27ac_6Hr,H3K27ac_24Hr) %>% gather(time,logFC,-type)

boxplot_dat$time <- factor(boxplot_dat$time, levels = c("H3K27ac_1Hr","H3K27ac_3Hr","H3K27ac_6Hr","H3K27ac_24Hr"))
boxplot_dat$type <- factor(boxplot_dat$type, levels = c("top1000","bot1000") )

##########
### generating box/violin plots per each metric type

pdf("H3K27ac_grouped_time.pdf",height = 6,width = 12)
boxplot_dat %>% ggplot(aes(x = time, y = logFC,fill = type)) +
  geom_violin( position = position_dodge(0.9), alpha = 0.5) + facet_wrap(.~time,scales = "free_x",nrow = 1)+
  ggpubr::stat_compare_means() + xlab("") + theme_classic(base_size =  15) + geom_boxplot(alpha = 0,width=0.6, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("red","#999999")) + geom_hline(yintercept = 0,linetype="dashed") + ggtitle("H3K27ac LogFC dm3-normalized")
dev.off()

boxplot_dat <- ep_dat_DMSO %>% filter(type != "Other") %>% dplyr::select(type,H3K18ac_1Hr,H3K18ac_3Hr,H3K18ac_6Hr,H3K18ac_24Hr) %>% gather(time,logFC,-type)

boxplot_dat$time <- factor(boxplot_dat$time, levels = c("H3K18ac_1Hr","H3K18ac_3Hr","H3K18ac_6Hr","H3K18ac_24Hr"))
boxplot_dat$type <- factor(boxplot_dat$type, levels = c("top1000","bot1000") )

pdf("H3K18ac_grouped_time.pdf",height = 6,width = 12)
boxplot_dat %>% ggplot(aes(x = time, y = logFC,fill = type)) +
  geom_violin( position = position_dodge(0.9), alpha = 0.5) + facet_wrap(.~time,scales = "free_x",nrow = 1)+
  ggpubr::stat_compare_means() + xlab("") + theme_classic(base_size =  15) + geom_boxplot(alpha = 0,width=0.6, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("red","#999999")) + geom_hline(yintercept = 0,linetype="dashed") + ggtitle("H3K18ac LogFC dm3-normalized")
dev.off() 

boxplot_dat <- ep_dat_DMSO %>% filter(type != "Other") %>% dplyr::select(type,H3K27me3_1Hr,H3K27me3_3Hr,H3K27me3_6Hr,H3K27me3_24Hr) %>% gather(time,logFC,-type)

boxplot_dat$time <- factor(boxplot_dat$time, levels = c("H3K27me3_1Hr","H3K27me3_3Hr","H3K27me3_6Hr","H3K27me3_24Hr"))
boxplot_dat$type <- factor(boxplot_dat$type, levels = c("top1000","bot1000") )

### H3K27me3 not used, many top and bottom genes had little to zero signal, and thus had NA logFCs. This is not surprising, as prior in figure 1 ChIP heatmap we show low overlap of H3K ac marks and H3K27me3

pdf("H3K27me3_grouped_time.pdf",height = 6,width = 12)
boxplot_dat %>% ggplot(aes(x = time, y = logFC,fill = type)) +
  geom_violin( position = position_dodge(0.9), alpha = 0.5) + facet_wrap(.~time,scales = "free_x",nrow = 1)+
  ggpubr::stat_compare_means() + xlab("") + theme_classic(base_size =  15) + geom_boxplot(alpha = 0,width=0.6, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("red","#999999")) + geom_hline(yintercept = 0,linetype="dashed") + ggtitle("H3K27me3 LogFC dm3-normalized")
dev.off()

boxplot_dat <- dat %>% filter(type != "Other") %>% dplyr::select(type,LogFC_1Hr,LogFC_3Hr,LogFC_6Hr,LogFC_24Hr) %>% gather(time,logFC,-type)

boxplot_dat$time <- factor(boxplot_dat$time, levels = c("LogFC_1Hr","LogFC_3Hr","LogFC_6Hr","LogFC_24Hr"))
boxplot_dat$type <- factor(boxplot_dat$type, levels = c("top1000","bot1000") )

pdf("RNA_Seq_grouped_time.pdf",height = 6,width = 12)
boxplot_dat %>% ggplot(aes(x = time, y = logFC,fill = type)) +
  geom_violin( position = position_dodge(0.9), alpha = 0.5) + facet_wrap(.~time,scales = "free_x",nrow = 1)+
  ggpubr::stat_compare_means() + xlab("") + theme_classic(base_size =  15) + geom_boxplot(alpha = 0,width=0.6, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("red","#999999")) + geom_hline(yintercept = 0,linetype="dashed") + ggtitle("RNA-Seq")

dev.off()
