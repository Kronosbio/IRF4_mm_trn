rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

gc()

setwd("~/irf4_mm_trn_code/figure_s1/")

library(tidyverse)
library(scales)
library(patchwork)
library(gridExtra)

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
#   [1] gridExtra_2.3   patchwork_1.2.0 scales_1.3.0    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
# [10] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
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
# [49] withr_3.0.0             UCSC.utils_1.0.0        bit64_4.0.5             timechange_0.3.0        XVector_0.44.0          httr_1.4.7             
# [55] bit_4.0.5               qvalue_2.36.0           hms_1.1.3               png_0.1-8               memoise_2.0.1           IRanges_2.38.0         
# [61] rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2               BiocGenerics_0.50.0     rstudioapi_0.16.0      
# [67] jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1                fs_1.6.4                zlibbioc_1.50.0   

'%!in%' <- function(x,y)!('%in%'(x,y))

### TF list

TF_list <- read.delim("../public_data/TFlist_NMid_hg19.txt", header=FALSE)

crc_degree <- read.csv("../ChIP_data/crc_degree.csv")

crc_degree$...1 <- NULL
crc_degree$X <- NULL

### identifying patient samples 

patients <- c("MM11.H3K27ac_DEGREE_TABLE.txt","MM24.H3K27ac_DEGREE_TABLE.txt","MM1.H3K27ac_DEGREE_TABLE.txt","MM2.H3K27ac_DEGREE_TABLE.txt",
              "MM3.H3K27ac_DEGREE_TABLE.txt","MM1_jia_DEGREE_TABLE.txt","MM4_jia_DEGREE_TABLE.txt","MM8.H3K27ac_DEGREE_TABLE.txt",
              "MM5.H3K27ac_DEGREE_TABLE.txt","MM7_jia_DEGREE_TABLE.txt","MM2_jia_DEGREE_TABLE.txt","MM6_jia_DEGREE_TABLE.txt",
              "MM3_jia_DEGREE_TABLE.txt","MM10_jia_DEGREE_TABLE.txt","MM5_jia_DEGREE_TABLE.txt","MM12.H3K27ac_DEGREE_TABLE.txt",
              "MM8_jia_DEGREE_TABLE.txt","MM25.H3K27ac_DEGREE_TABLE.txt","MM14.H3K27ac_DEGREE_TABLE.txt","MM9_jia_DEGREE_TABLE.txt"
)

crc_complete_pats <- crc_degree %>% filter(file %in% patients) 

t <- crc_degree %>% filter(file %!in% patients) 


#### aggregating cell line replicates together, TF must be observed in both replicates to be included in analysis. 

U266 <- t %>% filter(file %in% c("U266.H3K27ac_DEGREE_TABLE.txt","U266_H3K27ac_DEGREE_TABLE.txt")) 

U266 <- table(U266$Tf)[table(U266$Tf) > 1]
U266_comb  <- t %>% filter(file %in% c("U266_H3K27ac_DEGREE_TABLE.txt"))  %>% filter(Tf %in% names(U266))

t <- t %>% filter(file %!in% c("U266.H3K27ac_DEGREE_TABLE.txt","U266_H3K27ac_DEGREE_TABLE.txt")) 
t <- rbind(t,U266_comb)

rm(U266_comb,U266)

RPMI8226 <- t %>% filter(file %in% c("RPMI8226.H3K27ac_DEGREE_TABLE.txt","RPMI8226_H3K27ac_DEGREE_TABLE.txt")) 

RPMI8226 <- table(RPMI8226$Tf)[table(RPMI8226$Tf) > 1]
RPMI8226_comb  <- t %>% filter(file %in% c("RPMI8226_H3K27ac_DEGREE_TABLE.txt"))  %>% filter(Tf %in% names(RPMI8226))

t <- t %>% filter(file %!in% c("RPMI8226.H3K27ac_DEGREE_TABLE.txt","RPMI8226_H3K27ac_DEGREE_TABLE.txt")) 
t <- rbind(t,RPMI8226_comb)
rm(RPMI8226_comb,RPMI8226)

KMS12 <- t %>% filter(file %in% c("KMS-12-BM.H3K27ac_DEGREE_TABLE.txt","KMS12_H3K27ac_DEGREE_TABLE.txt")) 

KMS12 <- table(KMS12$Tf)[table(KMS12$Tf) > 1]
KMS12_comb  <- t %>% filter(file %in% c("KMS12_H3K27ac_DEGREE_TABLE.txt"))  %>% filter(Tf %in% names(KMS12))

t <- t %>% filter(file %!in% c("KMS-12-BM.H3K27ac_DEGREE_TABLE.txt","KMS12_H3K27ac_DEGREE_TABLE.txt")) 
t <- rbind(t,KMS12_comb)
rm(KMS12_comb,KMS12)

H929 <- t %>% filter(file %in% c("NCI-H929.H3K27ac_DEGREE_TABLE.txt","H929_H3K27ac_DEGREE_TABLE.txt")) 

H929 <- table(H929$Tf)[table(H929$Tf) > 1]
H929_comb  <- t %>% filter(file %in% c("H929_H3K27ac_DEGREE_TABLE.txt"))  %>% filter(Tf %in% names(H929))

t <- t %>% filter(file %!in% c("NCI-H929.H3K27ac_DEGREE_TABLE.txt","H929_H3K27ac_DEGREE_TABLE.txt")) 
t <- rbind(t,H929_comb)
rm(H929_comb,H929)

crc_complete_cells <- t

### running saturation analysis on patient and cell line samples individiually 

colnames(crc_complete_cells)[1] <- "SOURCE"
colnames(crc_complete_cells)[5] <- "ORIGIN"
colnames(crc_complete_pats)[1] <- "SOURCE"
colnames(crc_complete_pats)[5] <- "ORIGIN"

####################

temp <- crc_complete_pats %>% dplyr::select(SOURCE, ORIGIN) %>% unique()

temp <- table(temp$SOURCE) %>% as.data.frame()

tab_temp <- table(temp$Freq) %>% as.data.frame(stringsAsFactors = FALSE)

obs_tfs_pats <- tab_temp 

obs_tfs_pats$type <- "original"

obs_tfs_pats$Var1 <- as.numeric(obs_tfs_pats$Var1)
obs_tfs_pats$Freq <- as.numeric(obs_tfs_pats$Freq)
obs_tfs_pats <- obs_tfs_pats[order(obs_tfs_pats$Var1),]

obs_tfs_pats <- obs_tfs_pats%>% map_df(rev)
obs_tfs_pats$cumulative_observations <- cumsum(obs_tfs_pats$Freq)

### Running 1000 iterations for the patient samples, to see how often we would expect to see a TF given the number of observed number of TFs in a random sample.


for (i in 1:1000) {
  set.seed(i)
  
  temp <- crc_complete_pats %>% dplyr::select(SOURCE, ORIGIN) %>% unique() %>% dplyr::select(ORIGIN)
  
  temp$SOURCE <- sample(size = length(temp$ORIGIN),x = unique(TF_list$V2),replace = TRUE)

  temp <- table(temp$SOURCE) %>% as.data.frame()
  
  tab_temp <- table(temp$Freq) %>% as.data.frame(stringsAsFactors = FALSE)
  
  tab_temp$type <- "null"
  tab_temp$cumulative_observations <- 0
  
  obs_tfs_pats  <- rbind(obs_tfs_pats,tab_temp)

}

null_dist <- obs_tfs_pats %>% filter(type == "null") %>% group_by(Var1) %>% summarize(freq_mean = sum(Freq)/1000,
                                                                                      sample_SD = sd(Freq),
                                                                                      sample_se = sample_SD/sqrt(1000))

alpha = 0.05
degrees.freedom = 999
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)

### Calculating the margin of error for   for the null distribution

null_dist$margin.error <- t.score * null_dist$sample_se
null_dist$lower.bound <- null_dist$freq_mean - null_dist$margin.error
null_dist$upper.bound <- null_dist$freq_mean + null_dist$margin.error

null_dist$Var1 <- as.numeric(null_dist$Var1)
null_dist$freq_mean <- as.numeric(null_dist$freq_mean)
null_dist$type <- "null"

null_dist <- null_dist[null_dist$Var1,]
rownames(null_dist) <- null_dist$Var1
null_dist <- null_dist[order(null_dist$Var1),]

null_dist <- null_dist %>% map_df(rev)
null_dist$cumulative_observations <- cumsum(null_dist$freq_mean)

null_dist[is.na(null_dist)] <- 0

obs_tfs_pats <- obs_tfs_pats %>% filter(type != "null") 

null_dist$Var1 <- as.numeric(null_dist$Var1)/max(as.numeric(obs_tfs_pats$Var1))
obs_tfs_pats$Var1 <- as.numeric(obs_tfs_pats$Var1)/max(as.numeric(obs_tfs_pats$Var1))

a1 <- ggplot() + 
  geom_point(data = obs_tfs_pats,aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "black") +
  geom_line(data = obs_tfs_pats, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "black") +
  geom_point(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "#BFC9D0") +
  geom_line(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "#BFC9D0") +
  geom_errorbar(data = null_dist, aes(x = as.numeric(Var1),ymin=cumulative_observations-lower.bound,ymax=cumulative_observations+upper.bound),color = "#BFC9D0",width=0.05) + 
  ylab("Number of TFs") + xlab("Percent of Samples") +
  theme_classic(base_size = 16) + geom_vline(xintercept = 0.25,color = "red",linetype = "dashed") +
  ggtitle("Patients H3K27ac") + ylim(c(0,1500))

a2 <- ggplot() + 
  geom_point(data = obs_tfs_pats,aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "black") +
  geom_line(data = obs_tfs_pats, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "black") +
  geom_point(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "#BFC9D0") +
  geom_line(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "#BFC9D0") +
  geom_errorbar(data = null_dist, aes(x = as.numeric(Var1),ymin=cumulative_observations-lower.bound,ymax=cumulative_observations+upper.bound),color = "#BFC9D0",width=0.05) + 
  ylab("Number of TFs") + xlab("Percent of Samples") +
  theme_classic(base_size = 16) + geom_vline(xintercept = 0.25,color = "red",linetype = "dashed") +
  ggtitle("Patients H3K27ac - Zoomed In") + ylim(c(0,1500)) + coord_cartesian(ylim=c(0, 200)) 

a <- grid.arrange(a1, a2, nrow = 1)
###################

### doing same analysis now for cell line samples

temp <- crc_complete_cells %>% dplyr::select(SOURCE, ORIGIN) %>% unique()

temp <- table(temp$SOURCE) %>% as.data.frame()

tab_temp <- table(temp$Freq) %>% as.data.frame(stringsAsFactors = FALSE)

obs_tfs_cells <- tab_temp 

obs_tfs_cells$type <- "original"

obs_tfs_cells$Var1 <- as.numeric(obs_tfs_cells$Var1)
obs_tfs_cells$Freq <- as.numeric(obs_tfs_cells$Freq)
obs_tfs_cells <- obs_tfs_cells[order(obs_tfs_cells$Var1),]

obs_tfs_cells <- obs_tfs_cells %>% map_df(rev)
obs_tfs_cells$cumulative_observations <- cumsum(obs_tfs_cells$Freq)

for (i in 1:1000) {
  set.seed(i)
  
  temp <- crc_complete_cells %>% dplyr::select(SOURCE, ORIGIN) %>% unique() %>% dplyr::select(ORIGIN)
  
  temp$SOURCE <- sample(size = length(temp$ORIGIN),x = unique(TF_list$V2),replace = TRUE)
  
  temp <- table(temp$SOURCE) %>% as.data.frame()
  
  tab_temp <- table(temp$Freq) %>% as.data.frame(stringsAsFactors = FALSE)
  
  tab_temp$type <- "null"
  tab_temp$cumulative_observations <- 0
  
  obs_tfs_cells  <- rbind(obs_tfs_cells,tab_temp)
  
}

null_dist <- obs_tfs_cells %>% filter(type == "null") %>% group_by(Var1) %>% summarize(freq_mean = sum(Freq)/1000,
                                                                                      sample_SD = sd(Freq),
                                                                                      sample_se = sample_SD/sqrt(1000))

alpha = 0.05
degrees.freedom = 999
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)

null_dist$margin.error <- t.score * null_dist$sample_se
null_dist$lower.bound <- null_dist$freq_mean - null_dist$margin.error
null_dist$upper.bound <- null_dist$freq_mean + null_dist$margin.error

null_dist$Var1 <- as.numeric(null_dist$Var1)
null_dist$freq_mean <- as.numeric(null_dist$freq_mean)
null_dist$type <- "null"

null_dist <- null_dist[null_dist$Var1,]
rownames(null_dist) <- null_dist$Var1
null_dist <- null_dist[order(null_dist$Var1),]

null_dist <- null_dist %>% map_df(rev)
null_dist$cumulative_observations <- cumsum(null_dist$freq_mean)

null_dist[is.na(null_dist)] <- 0

obs_tfs_cells <- obs_tfs_cells %>% filter(type != "null") 

null_dist$Var1 <- as.numeric(null_dist$Var1)/max(as.numeric(obs_tfs_cells$Var1))
obs_tfs_cells$Var1 <- as.numeric(obs_tfs_cells$Var1)/max(as.numeric(obs_tfs_cells$Var1))

b1 <- ggplot() + 
  geom_point(data = obs_tfs_cells,aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "black") +
  geom_line(data = obs_tfs_cells, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "black") +
  geom_point(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "#BFC9D0") +
  geom_line(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "#BFC9D0") +
  geom_errorbar(data = null_dist, aes(x = as.numeric(Var1),ymin=cumulative_observations-lower.bound,ymax=cumulative_observations+upper.bound),color = "#BFC9D0",width=0.05) + 
  ylab("Number of TFs") +  xlab("Percent of Samples") +
  theme_classic(base_size = 16) + geom_vline(xintercept = 0.25,color = "red",linetype = "dashed") +
  ggtitle("Cell H3K27ac") + ylim(c(0,1500))

b2 <- ggplot() + 
  geom_point(data = obs_tfs_cells,aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "black") +
  geom_line(data = obs_tfs_cells, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "black") +
  geom_point(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 3, color = "#BFC9D0") +
  geom_line(data = null_dist, aes(x = as.numeric(Var1),y = cumulative_observations, color = type ), size = 2,color = "#BFC9D0") +
  geom_errorbar(data = null_dist, aes(x = as.numeric(Var1),ymin=cumulative_observations-lower.bound,ymax=cumulative_observations+upper.bound),color = "#BFC9D0",width=0.05) + 
  ylab("Number of TFs") + xlab("Percent of Samples") +
  theme_classic(base_size = 16) + geom_vline(xintercept = 0.25,color = "red",linetype = "dashed") +
  ggtitle("Cell H3K27ac - Zoomed In")+ ylim(c(0,1500)) + coord_cartesian(ylim=c(0, 200)) 


pdf("./fig_s1_sat_a.pdf",height = 5,width = 10)
grid.arrange(a1, a2, nrow = 1)
dev.off() 

pdf("./fig_s1_sat_b.pdf",height = 5,width = 10)
grid.arrange(b1, b2, nrow = 1)
dev.off() 

