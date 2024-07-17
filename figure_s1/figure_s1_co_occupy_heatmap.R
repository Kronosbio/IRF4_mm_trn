###### First defining common peaks between replicates of each ChIP-seq data
###### essentially if a peak had overlaps within both replicates, it was considered a 'common peak'
###### common peaks found for p300, IRF4, and DMSO H3K27ac/H3K18ac/H3K27me3 ChIP-seq data
###### These commons peaks were used for enhancerpromoter 

rm(list = ls())

if (is.integer(dev.list())) {
  #dev.off()
}
cat("\014")
set.seed(1)
library(tidyverse)
library(GenomicRanges)

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.20.0 circlize_0.4.16       GenomicRanges_1.56.0  GenomeInfoDb_1.40.0   IRanges_2.38.0        S4Vectors_0.42.0     
# [7] BiocGenerics_0.50.0   lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
# [13] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0      
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1           digest_0.6.35          
# [7] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.6           KEGGREST_1.44.0         RSQLite_2.3.6           magrittr_2.0.3         
# [13] compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4              data.table_1.15.4       bit_4.0.5              
# [19] RColorBrewer_1.1-3      plyr_1.8.9              BiocParallel_1.38.0     withr_3.0.0             fansi_1.0.6             GOSemSim_2.30.0        
# [25] colorspace_2.1-0        GO.db_3.19.1            iterators_1.0.14        scales_1.3.0            cli_3.6.2               crayon_1.5.2           
# [31] generics_0.1.3          rstudioapi_0.16.0       rjson_0.2.21            httr_1.4.7              reshape2_1.4.4          tzdb_0.4.0             
# [37] DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1             zlibbioc_1.50.0         splines_4.4.0          
# [43] parallel_4.4.0          AnnotationDbi_1.66.0    XVector_0.44.0          matrixStats_1.3.0       yulab.utils_0.1.4       vctrs_0.6.5            
# [49] Matrix_1.7-0            jsonlite_1.8.8          GetoptLong_1.0.5        hms_1.1.3               bit64_4.0.5             clue_0.3-65            
# [55] foreach_1.5.2           glue_1.7.0              codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.5           
# [61] shape_1.4.6.1           UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12
# [67] R6_2.5.1                doParallel_1.0.17       lattice_0.22-6          Biobase_2.64.0          png_0.1-8               memoise_2.0.1          
# [73] Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3         GlobalOptions_0.1.2  

'%!in%' <- function(x,y)!('%in%'(x,y))

setwd("~/irf4_mm_trn_code/figure_s1/")

### bed files coming from MACS1.4 peak calling

IRF4_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/31_0EX7_0229Kronos_MM1S-No-Treatment-1_IRF4_hs_i45_peaks.bed", header=FALSE)
IRF4_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/32_0EX8_0229Kronos_MM1S-No-Treatment-2_IRF4_hs_i48_peaks.bed", header=FALSE)
p300_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/33_0EX9_0229Kronos_MM1S-No-Treatment-1_p300_hs_i49_peaks.bed", header=FALSE)
p300_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/34_0EXA_0229Kronos_MM1S-No-Treatment-2_p300_hs_i50_peaks.bed", header=FALSE)

### creating GRanges objects for each ChIP-seq data set
gr1 <- GRanges(seqnames = IRF4_chip_1$V1, ranges = IRanges(start = IRF4_chip_1$V2, end = IRF4_chip_1$V3))
gr2 <- GRanges(seqnames = IRF4_chip_2$V1, ranges = IRanges(start = IRF4_chip_2$V2, end = IRF4_chip_2$V3))

### identifying overlaps and then creating a peak dataset for the common peaks. 
overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

IRF4_chip <- IRF4_chip_1[overlapping_indices, ]
rm(IRF4_chip_1,IRF4_chip_2)

gr1 <- GRanges(seqnames = p300_chip_1$V1, ranges = IRanges(start = p300_chip_1$V2, end = p300_chip_1$V3))
gr2 <- GRanges(seqnames = p300_chip_2$V1, ranges = IRanges(start = p300_chip_2$V2, end = p300_chip_2$V3))

overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

p300_chip <- p300_chip_1[overlapping_indices, ]

rm(p300_chip_1,p300_chip_2,gr1,gr2)

write_tsv(p300_chip,"../ChIP_data/p300_common_peaks.bed")
write_tsv(IRF4_chip,"../ChIP_data/IRF4_common_peaks.bed")

H3K18ac_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/01_0EWN_0229Kronos_MM1S-No-Treatment-1_H3K18Ac_hs-dm_i02_peaks.bed", header=FALSE)
H3K18ac_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/02_0EWO_0229Kronos_MM1S-No-Treatment-2_H3K18Ac_hs-dm_i04_peaks.bed", header=FALSE)

gr1 <- GRanges(seqnames = H3K18ac_chip_1$V1, ranges = IRanges(start = H3K18ac_chip_1$V2, end = H3K18ac_chip_1$V3))
gr2 <- GRanges(seqnames = H3K18ac_chip_2$V1, ranges = IRanges(start = H3K18ac_chip_2$V2, end = H3K18ac_chip_2$V3))

overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

H3K18ac_chip <- H3K18ac_chip_1[overlapping_indices, ]
rm(H3K18ac_chip_1,H3K18ac_chip_2,gr1,gr2)

write_tsv(p300_chip,"../ChIP_data/H3K18ac_DMSO_common_peaks.bed")

H3K18ac_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/01_0EWN_0229Kronos_MM1S-No-Treatment-1_H3K18Ac_hs-dm_i02_peaks.bed", header=FALSE)
H3K18ac_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/02_0EWO_0229Kronos_MM1S-No-Treatment-2_H3K18Ac_hs-dm_i04_peaks.bed", header=FALSE)

gr1 <- GRanges(seqnames = H3K18ac_chip_1$V1, ranges = IRanges(start = H3K18ac_chip_1$V2, end = H3K18ac_chip_1$V3))
gr2 <- GRanges(seqnames = H3K18ac_chip_2$V1, ranges = IRanges(start = H3K18ac_chip_2$V2, end = H3K18ac_chip_2$V3))

overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

H3K18ac_chip <- H3K18ac_chip_1[overlapping_indices, ]
rm(H3K18ac_chip_1,H3K18ac_chip_2,gr1,gr2)

write_tsv(H3K18ac_chip,"../ChIP_data/H3K18ac_DMSO_common_peaks.bed")

H3K27ac_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/11_0EWX_0229Kronos_MM1S-No-Treatment-1_H3K27Ac_hs-dm_i20_peaks.bed", header=FALSE)
H3K27ac_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/12_0EWY_0229Kronos_MM1S-No-Treatment-2_H3K27Ac_hs-dm_i21_peaks.bed", header=FALSE)

gr1 <- GRanges(seqnames = H3K27ac_chip_1$V1, ranges = IRanges(start = H3K27ac_chip_1$V2, end = H3K27ac_chip_1$V3))
gr2 <- GRanges(seqnames = H3K27ac_chip_2$V1, ranges = IRanges(start = H3K27ac_chip_2$V2, end = H3K27ac_chip_2$V3))

overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

H3K27ac_chip <- H3K27ac_chip_1[overlapping_indices, ]
rm(H3K27ac_chip_1,H3K27ac_chip_2,gr1,gr2)

write_tsv(H3K27ac_chip,"../ChIP_data/H3K27ac_DMSO_common_peaks.bed")

H3K27me3_chip_1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/21_0F27_0229Kronos_MM1S-No-Treatment-1_H3K27me3_hs-dm_i02_peaks.bed", header=FALSE)
H3K27me3_chip_2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/22_0F28_0229Kronos_MM1S-No-Treatment-2_H3K27me3_hs-dm_i04_peaks.bed", header=FALSE)

gr1 <- GRanges(seqnames = H3K27me3_chip_1$V1, ranges = IRanges(start = H3K27me3_chip_1$V2, end = H3K27me3_chip_1$V3))
gr2 <- GRanges(seqnames = H3K27me3_chip_2$V1, ranges = IRanges(start = H3K27me3_chip_2$V2, end = H3K27me3_chip_2$V3))

overlaps <- findOverlaps(gr1, gr2)
overlapping_indices <- queryHits(overlaps)

H3K27me3_chip <- H3K27me3_chip_1[overlapping_indices, ]
rm(H3K27me3_chip_1,H3K27me3_chip_2,gr1,gr2)

write_tsv(H3K27me3_chip,"../ChIP_data/H3K27me3_DMSO_common_peaks.bed")


#####################

#### Calculating the overlaps used within the ChIP overlap heatmap for figure S1

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(tidyverse)
library(circlize)
library(GenomicRanges)
library(ComplexHeatmap)

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
 #   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 # [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
 # [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
 # 
 # time zone: UTC
 # tzcode source: system (glibc)
 # 
 # attached base packages:
 #   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
 # 
 # other attached packages:
 #   [1] ComplexHeatmap_2.20.0 circlize_0.4.16       GenomicRanges_1.56.0  GenomeInfoDb_1.40.0   IRanges_2.38.0        S4Vectors_0.42.0     
 # [7] BiocGenerics_0.50.0   lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
 # [13] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0      
 # 
 # loaded via a namespace (and not attached):
 #   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1           digest_0.6.35          
 # [7] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.6           KEGGREST_1.44.0         RSQLite_2.3.6           magrittr_2.0.3         
 # [13] compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4              data.table_1.15.4       bit_4.0.5              
 # [19] RColorBrewer_1.1-3      plyr_1.8.9              BiocParallel_1.38.0     withr_3.0.0             fansi_1.0.6             GOSemSim_2.30.0        
 # [25] colorspace_2.1-0        GO.db_3.19.1            iterators_1.0.14        scales_1.3.0            cli_3.6.2               crayon_1.5.2           
 # [31] generics_0.1.3          rstudioapi_0.16.0       rjson_0.2.21            httr_1.4.7              reshape2_1.4.4          tzdb_0.4.0             
 # [37] DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1             zlibbioc_1.50.0         splines_4.4.0          
 # [43] parallel_4.4.0          AnnotationDbi_1.66.0    XVector_0.44.0          matrixStats_1.3.0       yulab.utils_0.1.4       vctrs_0.6.5            
 # [49] Matrix_1.7-0            jsonlite_1.8.8          GetoptLong_1.0.5        hms_1.1.3               bit64_4.0.5             clue_0.3-65            
 # [55] foreach_1.5.2           glue_1.7.0              codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.5           
 # [61] shape_1.4.6.1           UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12
 # [67] R6_2.5.1                doParallel_1.0.17       lattice_0.22-6          Biobase_2.64.0          png_0.1-8               memoise_2.0.1          
 # [73] Rcpp_1.0.12             fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3         GlobalOptions_0.1.2  

 
 setwd("~/irf4_mm_trn_code/figure_s1")

MM1S_H3K27ac <- read.delim("../ChIP_data/H3K27ac_DMSO_common_peaks.bed", header=TRUE)
MM1S_H3K27me3 <- read.delim("../ChIP_data/H3K27me3_DMSO_common_peaks.bed", header=TRUE)
MM1S_H3K18ac <- read.delim("../ChIP_data/H3K18ac_DMSO_common_peaks.bed", header=TRUE)

IRF4_chip <- read.delim("../ChIP_data/IRF4_common_peaks.bed", header=TRUE)
p300_chip<- read.delim("../ChIP_data/p300_common_peaks.bed", header=TRUE)

p300_gr <- GRanges(seqnames = p300_chip$V1, ranges = IRanges(start = p300_chip$V2, end = p300_chip$V3))
IRF4_gr <- GRanges(seqnames = IRF4_chip$V1, ranges = IRanges(start = IRF4_chip$V2, end = IRF4_chip$V3))
H3K27ac_gr <- GRanges(seqnames = MM1S_H3K27ac$V1, ranges = IRanges(start = MM1S_H3K27ac$V2, end = MM1S_H3K27ac$V3))
H3K18ac_gr <- GRanges(seqnames = MM1S_H3K18ac$V1, ranges = IRanges(start = MM1S_H3K18ac$V2, end = MM1S_H3K18ac$V3))
H3K27me3_gr <- GRanges(seqnames = MM1S_H3K27me3$V1, ranges = IRanges(start = MM1S_H3K27me3$V2, end = MM1S_H3K27me3$V3))

type <- c("p300","IRF4","H3K27ac","H3K18ac","H3K27me3")

mat_data <- matrix(ncol = length(type),nrow = length(type))

colnames(mat_data) <- type
rownames(mat_data) <- type

### Tracking overlap for each individual ChIP type, and calculating the overlap

gr_ranges <- list(p300_gr,IRF4_gr,H3K27ac_gr,H3K18ac_gr,H3K27me3_gr)
for (i in 1:length(type)) {
  print(type[i])
  print(length(gr_ranges[[i]]))
  
        for (j in 1:length(type)) {
    overlaps <- findOverlaps(gr_ranges[[i]],gr_ranges[[j]])
    overlapping_indices <- queryHits(overlaps) %>% unique()
    mat_data[i,j] <- length(overlapping_indices)/length(gr_ranges[[i]])
  }
}


### creating the heatmap with 3 different heatmap colors and combining them in illustrator 

col_fun = colorRamp2(c(0,1), c("white","cyan"))
pdf("./figure_s1_co_occupy_heatmap_p300.pdf",height = 4,width = 5)
Heatmap(mat_data, cluster_rows = FALSE,cluster_columns = FALSE, col = col_fun)
dev.off()

col_fun = colorRamp2(c(0,1), c("white","magenta"))
pdf("./figure_s1_co_occupy_heatmap_IRF4.pdf",height = 4,width = 5)
Heatmap(mat_data, cluster_rows = FALSE,cluster_columns = FALSE, col = col_fun)
dev.off()

col_fun = colorRamp2(c(0,1), c("white","black"))
pdf("./figure_s1_co_occupy_heatmap_H3K.pdf",height = 4,width = 5)
Heatmap(mat_data, cluster_rows = FALSE,cluster_columns = FALSE, col = col_fun)
dev.off()
