#### this code was used for the RNA-seq related to KB528 and other standard of care compounds.

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)
setwd("~/irf4_mm_trn_code/figure_4/")

library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# time zone: UTC
# tzcode source: system (glibc)
#
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] patchwork_1.2.0       circlize_0.4.16       lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1
# [6] dplyr_1.1.4           purrr_1.0.2           readr_2.1.5           tidyr_1.3.1           tibble_3.2.1
# [11] ggplot2_3.5.1         tidyverse_2.0.0       ComplexHeatmap_2.20.0 org.Hs.eg.db_3.19.1   AnnotationDbi_1.66.0
# [16] IRanges_2.38.0        S4Vectors_0.42.0      Biobase_2.64.0        BiocGenerics_0.50.0
#
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        viridisLite_0.4.2       HDO.db_0.99.1           blob_1.2.4
# [5] viridis_0.6.5           Biostrings_2.72.0       fastmap_1.1.1           digest_0.6.35
# [9] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.6           KEGGREST_1.44.0
# [13] RSQLite_2.3.6           magrittr_2.0.3          compiler_4.4.0          rlang_1.1.3
# [17] tools_4.4.0             utf8_1.2.4              data.table_1.15.4       bit_4.0.5
# [21] plyr_1.8.9              RColorBrewer_1.1-3      BiocParallel_1.38.0     withr_3.0.0
# [25] fansi_1.0.6             GOSemSim_2.30.0         colorspace_2.1-0        GO.db_3.19.1
# [29] scales_1.3.0            iterators_1.0.14        cli_3.6.2               crayon_1.5.2
# [33] generics_0.1.3          rstudioapi_0.16.0       tzdb_0.4.0              httr_1.4.7
# [37] reshape2_1.4.4          rjson_0.2.21            DBI_1.2.2               qvalue_2.36.0
# [41] cachem_1.0.8            DOSE_3.30.1             zlibbioc_1.50.0         splines_4.4.0
# [45] parallel_4.4.0          XVector_0.44.0          matrixStats_1.3.0       yulab.utils_0.1.4
# [49] vctrs_0.6.5             Matrix_1.7-0            jsonlite_1.8.8          hms_1.1.3
# [53] GetoptLong_1.0.5        bit64_4.0.5             ggrepel_0.9.5           clue_0.3-65
# [57] foreach_1.5.2           glue_1.7.0              codetools_0.2-20        cowplot_1.1.3
# [61] shape_1.4.6.1           stringi_1.8.4           gtable_0.3.5            GenomeInfoDb_1.40.0
# [65] UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0
# [69] GenomeInfoDbData_1.2.12 R6_2.5.1                doParallel_1.0.17       lattice_0.22-6
# [73] png_0.1-8               memoise_2.0.1           Rcpp_1.0.12             fastmatch_1.1-4
# [77] gridExtra_2.3           fs_1.6.4                pkgconfig_2.0.3         GlobalOptions_0.1.2

genematch <- function(dat) {
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL')
  genes_IDS <- mapIds(org.Hs.eg.db, entrez_IDS, 'SYMBOL', 'ENTREZID')
  
  genes_IDS <- genes_IDS %>% unlist()
  t <- as.data.frame(genes_IDS)
  t$entrez_IDS <- rownames(t)
  
  t2 <- as.data.frame(entrez_IDS)
  t2$GENE  <- rownames(t2)
  key_name <- merge(t, t2)
  
  key_name$entrez_IDS <- NULL
  dat <- merge(key_name, dat)
  
  rm(t, t2, key_name)
  dat$GENE <- NULL
  colnames(dat)[1] <- "GENE"
  return(dat)
}

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

#### Reading in data, 6 hours was not used in the manuscript, as most effects were seen at 24 hours.

MM1S_KB528_6 <- read.delim("../RNA_data/MM1S_SOC/MM1S_SOC_analysis_KB528_0_vs_KB528_6_exprs_matrix.txt")
MM1S_KB528_24 <- read.delim("../RNA_data/MM1S_SOC/MM1S_SOC_analysis_KB528_0_vs_KB528_24_exprs_matrix.txt")

MM1S_KB528_24$GENE <- row.names(MM1S_KB528_24)
MM1S_KB528_6$GENE <- row.names(MM1S_KB528_6)

colnames(MM1S_KB528_6)[3:4] <- c("Log2FC_KB528_6Hr", "p_KB528_6Hr")
colnames(MM1S_KB528_24)[3:4] <- c("Log2FC_KB528_24Hr", "p_KB528_24Hr")

MM1S_CCS1477_24 <- read.delim(
  "../RNA_data/MM1S_SOC/MM1S_SOC_analysis_CCS1477_0_vs_CCS1477_24_exprs_matrix.txt"
)
MM1S_CCS1477_6 <- read.delim(
  "../RNA_data/MM1S_SOC/MM1S_SOC_analysis_CCS1477_0_vs_CCS1477_6_exprs_matrix.txt"
)

MM1S_CCS1477_24$GENE <- row.names(MM1S_CCS1477_24)
MM1S_CCS1477_6$GENE <- row.names(MM1S_CCS1477_6)

colnames(MM1S_CCS1477_6)[3:4] <- c("Log2FC_CCS1477_6Hr", "p_CCS1477_6Hr")
colnames(MM1S_CCS1477_24)[3:4] <- c("Log2FC_CCS1477_24Hr", "p_CCS1477_24Hr")

MM1S_Lenalidomide_6 <- read.delim(
  "../RNA_data/MM1S_SOC/MM1S_SOC_analysis_Lenalidomide_0_vs_Lenalidomide_6_exprs_matrix.txt"
)
MM1S_Lenalidomide_24 <- read.delim(
  "../RNA_data/MM1S_SOC/MM1S_SOC_analysis_Lenalidomide_0_vs_Lenalidomide_24_exprs_matrix.txt"
)


MM1S_Lenalidomide_24$GENE <- row.names(MM1S_Lenalidomide_24)
MM1S_Lenalidomide_6$GENE <- row.names(MM1S_Lenalidomide_6)

colnames(MM1S_Lenalidomide_6)[3:4] <- c("Log2FC_Lenalidomide_6Hr", "p_Lenalidomide_6Hr")
colnames(MM1S_Lenalidomide_24)[3:4] <- c("Log2FC_Lenalidomide_24Hr", "p_Lenalidomide_24Hr")

MM1S_JQ1_6 <- read.delim("../RNA_data/MM1S_SOC/MM1S_SOC_analysis_JQ1_0_vs_JQ1_6_exprs_matrix.txt")
MM1S_JQ1_24 <- read.delim("../RNA_data/MM1S_SOC/MM1S_SOC_analysis_JQ1_0_vs_JQ1_24_exprs_matrix.txt")

MM1S_JQ1_24$GENE <- row.names(MM1S_JQ1_24)
MM1S_JQ1_6$GENE <- row.names(MM1S_JQ1_6)

colnames(MM1S_JQ1_6)[3:4] <- c("Log2FC_JQ1_6Hr", "p_JQ1_6Hr")
colnames(MM1S_JQ1_24)[3:4] <- c("Log2FC_JQ1_24Hr", "p_JQ1_24Hr")

MM1S_GNE781_6 <- read.delim("../RNA_data/MM1S_SOC/MM1S_SOC_analysis_GNE781_0_vs_GNE781_6_exprs_matrix.txt")
MM1S_GNE781_24 <- read.delim(
  "../RNA_data/MM1S_SOC/MM1S_SOC_analysis_GNE781_0_vs_GNE781_24_exprs_matrix.txt"
)

MM1S_GNE781_24$GENE <- row.names(MM1S_GNE781_24)
MM1S_GNE781_6$GENE <- row.names(MM1S_GNE781_6)

colnames(MM1S_GNE781_6)[3:4] <- c("Log2FC_GNE781_6Hr", "p_GNE781_6Hr")
colnames(MM1S_GNE781_24)[3:4] <- c("Log2FC_GNE781_24Hr", "p_GNE781_24Hr")

dfs <- list(
  MM1S_KB528_6,
  MM1S_KB528_24,
  MM1S_GNE781_6,
  MM1S_GNE781_24,
  MM1S_JQ1_6,
  MM1S_JQ1_24,
  MM1S_Lenalidomide_6,
  MM1S_Lenalidomide_24,
  MM1S_CCS1477_6,
  MM1S_CCS1477_24
)

dat <- dfs %>% purrr::reduce(full_join, by = "GENE") %>% as.data.frame()

#### used for supplemental tables
# write_csv(dat,"./MM1S_SOC_logFC.csv")

rm(
  dfs,
  MM1S_KB528_6,
  MM1S_KB528_24,
  MM1S_GNE781_6,
  MM1S_GNE781_24,
  MM1S_JQ1_6,
  MM1S_JQ1_24,
  MM1S_Lenalidomide_6,
  MM1S_Lenalidomide_24,
  MM1S_CCS1477_6,
  MM1S_CCS1477_24
)

### Fedele et al. IRF4 KD signature
MM1S_sgIRF4_DE <- read.csv("../public_data/MM1S_sgIRF4_DE_2.csv")

### Filtering for genes that are differential expressed in at least 3 conditions
across_genes <- table(MM1S_sgIRF4_DE$Symbol)[table(MM1S_sgIRF4_DE$Symbol) > 3]

MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %in% names(across_genes))
MM1S_sgIRF4_DE$comb <- paste(MM1S_sgIRF4_DE$Symbol, MM1S_sgIRF4_DE$Type)

### Removing 3 genes that were not consistent in direction.
MM1S_sgIRF4_DE <- MM1S_sgIRF4_DE %>% filter(Symbol %!in% c("ARHGEF3", "CDKN1A", "RHOB")) %>% dplyr::select(Symbol, Type) %>% unique()

MM1S_siIRF4 <- read.delim("../RNA_data/IRF4_KD/_siIRF4 MM1S_vs_CTRL MM1S_exprs_matrix.txt")

MM1S_siIRF4$LOG2_FOLD_CHANGE <- -MM1S_siIRF4$LOG2_FOLD_CHANGE
MM1S_siIRF4$GENE <- row.names(MM1S_siIRF4)

colnames(MM1S_siIRF4)[3:4] <- c("Log2FC_si", "p_si")

dat <- merge(MM1S_siIRF4, dat)
dat_IRF4 <- dat %>% filter(GENE == "IRF4")

dat <- dat %>% filter(GENE %in% c("IRF4", MM1S_sgIRF4_DE$Symbol))
dat$type <- "Other"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- "Down"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- "Up"
dat$type[dat$GENE == "IRF4"] <- "IRF4"

dat$alpha_ <- 0
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "down") %>% pull(Symbol))] <- 0.6
dat$alpha_[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- 0.6

heatmap_dat <- dat %>% filter(type != "Other") %>%
  dplyr::select(GENE, type, colnames(dat)[grep("Log2FC", colnames(dat))])

heatmap_dat <- heatmap_dat %>%  dplyr::select(GENE, type, Log2FC_si, colnames(heatmap_dat)[grep("24Hr", colnames(heatmap_dat))])
rownames(heatmap_dat) <- heatmap_dat$GENE

column_ha = columnAnnotation(type = heatmap_dat$type, col = list(type = c(
  "IRF4" = "magenta",
  "Up" = "#e66a78",
  "Down" = "#2171B5"
)))

###colors changed in illustrator to match aesthetics with other figures

row_ha = rowAnnotation(type = colnames(heatmap_dat)[-c(1, 2)], col = list(
  type = c(
    "Log2FC_si" = "magenta",
    "Log2FC_KB528_24Hr" = "#3A53A4",
    "Log2FC_GNE781_24Hr" = "#F58020",
    "Log2FC_JQ1_24Hr" = "#32B44A",
    "Log2FC_Lenalidomide_24Hr" = "black",
    "Log2FC_CCS1477_24Hr" = "gray"
  )
))

pdf(
  "./figS4_IRF4_signature_heatmap_soc.pdf",
  height = 4,
  width = 15
)
Heatmap(
  t(heatmap_dat[, -c(1, 2)]),
  clustering_distance_rows =  "spearman",
  clustering_distance_columns =  "manhattan",
  right_annotation = row_ha,
  top_annotation = column_ha
)
dev.off()

#########
### plotting the the IRF4 Log2FC vs correlation with siIRF4 signature

dat <- dat %>% filter(GENE %in% MM1S_sgIRF4_DE$Symbol)

dat$type <- "Down"
dat$type[dat$GENE %in% (MM1S_sgIRF4_DE %>% filter(Type == "up") %>% pull(Symbol))] <- "Up"

temp <- dat %>% dplyr::select(
  GENE,
  Log2FC_si,
  Log2FC_KB528_24Hr,
  Log2FC_GNE781_24Hr,
  Log2FC_JQ1_24Hr,
  Log2FC_Lenalidomide_24Hr,
  Log2FC_CCS1477_24Hr
)

#### correlation df with siIRF4 signature

cor_mat_IRF4_sig <- cor(temp[-1]) %>% as.data.frame()

cor_mat_IRF4_sig$Sample <- rownames(cor_mat_IRF4_sig)

dat_IRF4 <- t(dat_IRF4[-1]) %>% as.data.frame()
dat_IRF4$Sample <- rownames(dat_IRF4)
cor_mat_IRF4_sig$Sample <- rownames(cor_mat_IRF4_sig)

plot_dat <- merge(cor_mat_IRF4_sig, dat_IRF4)

plot_dat$Sample
plot_dat$Time <- c("24_Hr", "24_Hr", "24_Hr", "24_Hr", "24_Hr", "48Hour")

plot_dat <- plot_dat %>% filter(Time != "48Hour")

p <- plot_dat %>% ggplot(aes(
  x = V1,
  y = Log2FC_si,
  shape = Time,
  label = Sample,
  fill = Sample,
  color = Sample
)) +
  ###colors changed in illustrator to match aesthetics with other figures
  geom_point(size = 5, color = "black") +
  ggrepel::geom_label_repel() + theme_classic(base_size = 20) +
  scale_fill_manual(values = c("gray", "#F58020", "#32B44A", "#3A53A4", "black")) +
  scale_color_manual(values = c("black", "white", "white", "white", "white")) + theme(legend.position = "none") +
  ylab("Correlation (r) LogFC siIRF4 Signature") +
  xlab("Log FC IRF4")

pdf("fig4_IRF4_signature_IRF4.pdf",
    height = 6,
    width = 6)
p
dev.off()
