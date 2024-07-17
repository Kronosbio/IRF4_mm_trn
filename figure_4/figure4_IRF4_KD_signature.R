rm(list = ls())
if (is.integer(dev.list())) {
  dev.off()
}

cat("\014")
set.seed(1)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(patchwork)
library(ggpubr)

# >  sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# time zone: UTC
# tzcode source: system (glibc)
#
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggpubr_0.6.0          patchwork_1.2.0       circlize_0.4.16       lubridate_1.9.3       forcats_1.0.0
# [6] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.5           tidyr_1.3.1
# [11] tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0       ComplexHeatmap_2.20.0
#
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               rlang_1.1.3             magrittr_2.0.3          clue_0.3-65
# [5] GetoptLong_1.0.5        DOSE_3.30.1             matrixStats_1.3.0       compiler_4.4.0
# [9] RSQLite_2.3.6           png_0.1-8               vctrs_0.6.5             reshape2_1.4.4
# [13] pkgconfig_2.0.3         shape_1.4.6.1           crayon_1.5.2            fastmap_1.1.1
# [17] backports_1.4.1         magick_2.8.3            XVector_0.44.0          utf8_1.2.4
# [21] HDO.db_0.99.1           tzdb_0.4.0              UCSC.utils_1.0.0        bit_4.0.5
# [25] zlibbioc_1.50.0         cachem_1.0.8            GenomeInfoDb_1.40.0     jsonlite_1.8.8
# [29] blob_1.2.4              BiocParallel_1.38.0     broom_1.0.5             parallel_4.4.0
# [33] cluster_2.1.6           R6_2.5.1                stringi_1.8.4           RColorBrewer_1.1-3
# [37] car_3.1-2               GOSemSim_2.30.0         Rcpp_1.0.12             iterators_1.0.14
# [41] IRanges_2.38.0          Matrix_1.7-0            splines_4.4.0           timechange_0.3.0
# [45] tidyselect_1.2.1        abind_1.4-5             qvalue_2.36.0           rstudioapi_0.16.0
# [49] doParallel_1.0.17       codetools_0.2-20        lattice_0.22-6          plyr_1.8.9
# [53] Biobase_2.64.0          withr_3.0.0             KEGGREST_1.44.0         Biostrings_2.72.0
# [57] pillar_1.9.0            carData_3.0-5           foreach_1.5.2           stats4_4.4.0
# [61] generics_0.1.3          S4Vectors_0.42.0        hms_1.1.3               munsell_0.5.1
# [65] scales_1.3.0            glue_1.7.0              tools_4.4.0             data.table_1.15.4
# [69] fgsea_1.30.0            ggsignif_0.6.4          fs_1.6.4                fastmatch_1.1-4
# [73] cowplot_1.1.3           Cairo_1.6-2             AnnotationDbi_1.66.0    colorspace_2.1-0
# [77] GenomeInfoDbData_1.2.12 cli_3.6.2               fansi_1.0.6             gtable_0.3.5
# [81] rstatix_0.7.2           yulab.utils_0.1.4       digest_0.6.35           BiocGenerics_0.50.0
# [85] rjson_0.2.21            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7
# [89] GlobalOptions_0.1.2     GO.db_3.19.1            bit64_4.0.5

setwd("~/irf4_mm_trn_code/figure_4/")
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

MM1S_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_siIRF4 MM1S_vs_CTRL MM1S_exprs_matrix.txt")
MM1R_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_CTRL MM1R_vs_siIRF4 MM1R_exprs_matrix.txt")
KMS12BM_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_siIRF4 KMS12BM_vs_CTRL KMS12BM_exprs_matrix.txt")
KMS12PE_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_CTRL KMS12PE_vs_siIRF4 KMS12PE_exprs_matrix.txt")

#### these two need to be flipped in terms of direction. Their Log2FC had CTRL and siIRF4 flipped.
MM1S_siIRF4$LOG2_FOLD_CHANGE <- -MM1S_siIRF4$LOG2_FOLD_CHANGE
KMS12BM_siIRF4$LOG2_FOLD_CHANGE <- -KMS12BM_siIRF4$LOG2_FOLD_CHANGE

MM1S_siIRF4$GENE <- row.names(MM1S_siIRF4)
MM1R_siIRF4$GENE <- row.names(MM1R_siIRF4)
KMS12BM_siIRF4$GENE <- row.names(KMS12BM_siIRF4)
KMS12PE_siIRF4$GENE <- row.names(KMS12PE_siIRF4)

colnames(MM1S_siIRF4)[3:4] <- c("Log2FC_MM1S_siIRF4", "p_MM1S")
colnames(MM1R_siIRF4)[3:4] <- c("Log2FC_MM1R_siIRF4", "p_MM1R")
colnames(KMS12BM_siIRF4)[3:4] <- c("Log2FC_KMS12BM_siIRF4", "p_KMS12BM")
colnames(KMS12PE_siIRF4)[3:4] <- c("Log2FC_KMS12PE_siIRF4", "p_KMS12PE")

dat <- merge(MM1R_siIRF4, MM1S_siIRF4)
dat <- merge(dat, KMS12BM_siIRF4)
dat <- merge(dat, KMS12PE_siIRF4)

### Fedele et al. IRF4 KD signature
MM1S_sgIRF4_DE <- read.csv("../public_data/MM1S_sgIRF4_DE_2.csv")

### Filtering for genes that are differential expressed in at least 3 conditions
across_genes <- table(MM1S_sgIRF4_DE$Symbol)[table(MM1S_sgIRF4_DE$Symbol) > 3]

MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %in% names(across_genes))
MM1S_sgIRF4_DE$comb <- paste(MM1S_sgIRF4_DE$Symbol, MM1S_sgIRF4_DE$Type)

### Removing 3 genes that were not consistent in direction.
MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %!in% c("ARHGEF3", "CDKN1A", "RHOB")) %>% dplyr::select(Symbol, Type) %>% unique()

MM1S_KB528 <- read.delim(
  "../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_24_MM1S_exprs_matrix.txt"
)

RPMI8226_KB528 <- read.delim(
  "../RNA_data/KB528_time_series/KB528_RPMI8226_analysis_DMSO_RPMI8226_vs_KB528_24Hr_RPMI8226_exprs_matrix.txt"
)

MM1S_KB528$GENE <- row.names(MM1S_KB528)
RPMI8226_KB528$GENE <- row.names(RPMI8226_KB528)

colnames(MM1S_KB528)[3:4] <- c("Log2FC_KB528_24", "p_KB528_24")
colnames(RPMI8226_KB528)[3:4] <- c("Log2FC_KB528_24_RPMI8226", "p_KB528_24_RPMI8226")

dat <- merge(dat, MM1S_KB528)
dat <- merge(dat, RPMI8226_KB528)

dat$type <- "Other"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- "Down"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- "Up"
dat$type[dat$GENE == "IRF4"] <- "IRF4"

dat$alpha_ <- 0
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- 0.6
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- 0.6
dat$alpha_[dat$GENE == "IRF4"] <- 0.6

heatmap_dat <- dat %>% filter(type != "Other") %>%
  dplyr::select(GENE, type, colnames(dat)[grep("Log2FC", colnames(dat))])
rownames(heatmap_dat) <- heatmap_dat$GENE

column_ha = columnAnnotation(type = heatmap_dat$type, col = list(type = c(
  "Up" = "#e66a78",
  "IRF4"  = "magenta",
  "Down" = "#2171B5"
)))
Heatmap(t(heatmap_dat[, -c(1, 2)]),
        clustering_distance_rows =  "pearson",
        top_annotation = column_ha)

pdf("./fig4_IRF4_signature_heatmap.pdf",
    height = 4,
    width = 10)
Heatmap(
  t(heatmap_dat[, -c(1, 2)]),
  clustering_distance_rows =  "pearson",
  clustering_distance_columns =  "manhattan",
  top_annotation = column_ha
)
dev.off()

##################
rm(list = ls())
if (is.integer(dev.list())) {
  dev.off()
}

cat("\014")
set.seed(1)

setwd("~/irf4_mm_trn_code/figure_4/")
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

MM1S_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_siIRF4 MM1S_vs_CTRL MM1S_exprs_matrix.txt")

MM1S_siIRF4$LOG2_FOLD_CHANGE <- -MM1S_siIRF4$LOG2_FOLD_CHANGE

MM1S_siIRF4$GENE <- row.names(MM1S_siIRF4)

colnames(MM1S_siIRF4)[3:4] <- c("Log2FC_MM1S_siIRF4", "p_MM1S")

MM1S_sgIRF4_DE <- read.csv("../public_data/MM1S_sgIRF4_DE_2.csv")
across_genes <- table(MM1S_sgIRF4_DE$Symbol)[table(MM1S_sgIRF4_DE$Symbol) > 3]

MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %in% names(across_genes))
MM1S_sgIRF4_DE$comb <- paste(MM1S_sgIRF4_DE$Symbol, MM1S_sgIRF4_DE$Type)

MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %!in% c("ARHGEF3", "CDKN1A", "RHOB")) %>% dplyr::select(Symbol, Type) %>% unique()

###write_csv(MM1S_sgIRF4_DE,"./Fedele_et_al_2021_MM1S_sgIRF4_DE.csv")

MM1S_KB528 <- read.delim(
  "../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_24_MM1S_exprs_matrix.txt"
)
MM1S_KB528$GENE <- row.names(MM1S_KB528)

colnames(MM1S_KB528)[3:4] <- c("Log2FC_KB528_24", "p_KB528_24")

dat <- merge(MM1S_siIRF4, MM1S_KB528)

dat$type <- "Other"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- "Down"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- "Up"

dat$alpha_ <- 0
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- 0.6
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- 0.6

dat_t <- dat %>% filter(type != "Other")
cor.test(dat_t$Log2FC_MM1S_siIRF4, dat_t$Log2FC_KB528_24)

# Pearson's product-moment correlation
# 
# data:  dat_t$Log2FC_MM1S_siIRF4 and dat_t$Log2FC_KB528_24
# t = 15.775, df = 67, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8240605 0.9291110
# sample estimates:
#       cor 
# 0.8876204 

cor.test(dat$Log2FC_MM1S_siIRF4, dat$Log2FC_KB528_24)

# Pearson's product-moment correlation
# 
# data:  dat$Log2FC_MM1S_siIRF4 and dat$Log2FC_KB528_24
# t = 79.159, df = 11402, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5835577 0.6072495
# sample estimates:
#       cor 
# 0.5955331 

plot1 <- dat %>% ggplot(aes(y = Log2FC_KB528_24, x = Log2FC_MM1S_siIRF4)) +
  geom_point(alpha = 0.6, aes(color = type, size = 2)) +
  scale_color_manual(values = c("#0072B250", "gray", "#e66a78")) +
  theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", color = "darkgray") +
  theme(legend.position = "none") + xlab("Log2FC siIRF4") + xlim(c(-7, 7)) + ylim(c(-7.11, 7.11))

plot1 <- plot1 + geom_point(data = dat_t,
                            aes(
                              y = Log2FC_KB528_24,
                              x = Log2FC_MM1S_siIRF4,
                              color = type,
                              size = 2
                            )) +
  scale_color_manual(values = c("#2171B5", "gray", "#e66a78"))
plot1
dens1 <- ggplot(dat_t, aes(x = Log2FC_MM1S_siIRF4, fill = type)) +
  geom_density(alpha = 0.4) + xlim(c(-7.11, 7.11)) +
  theme_void() +
  theme(legend.position = "none") + scale_fill_manual(values = c("#2171B5", "#e66a78"))

dens2 <- ggplot(dat_t , aes(x =  Log2FC_KB528_24, fill = type)) +
  geom_density(alpha = 0.4) + xlim(c(-7.11, 7.11)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip() + scale_fill_manual(values = c("#2171B5", "#e66a78"))

c <- dens1 + plot_spacer() + plot1 + dens2 +
  plot_layout(
    ncol = 2,
    nrow = 2,
    widths = c(4, 1),
    heights = c(1, 4)
  )

pdf("./fig4_IRF4_signature_point.pdf",
    height = 5,
    width = 5)
c
dev.off()
