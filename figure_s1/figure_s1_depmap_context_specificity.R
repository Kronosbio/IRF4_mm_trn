rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

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
# [10] tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.44.0         fastmatch_1.1-4         gtable_0.3.5            GOSemSim_2.30.0         ggrepel_0.9.5           Biobase_2.64.0         
# [7] lattice_0.22-6          tzdb_0.4.0              vctrs_0.6.5             tools_4.4.0             generics_0.1.3          yulab.utils_0.1.4      
# [13] stats4_4.4.0            parallel_4.4.0          fansi_1.0.6             AnnotationDbi_1.66.0    RSQLite_2.3.6           blob_1.2.4             
# [19] pkgconfig_2.0.3         Matrix_1.7-0            data.table_1.15.4       S4Vectors_0.42.0        lifecycle_1.0.4         GenomeInfoDbData_1.2.12
# [25] HDO.db_0.99.1           compiler_4.4.0          Biostrings_2.72.0       munsell_0.5.1           fgsea_1.30.0            codetools_0.2-20       
# [31] DOSE_3.30.1             GenomeInfoDb_1.40.0     pillar_1.9.0            crayon_1.5.2            GO.db_3.19.1            BiocParallel_1.38.0    
# [37] cachem_1.0.8            tidyselect_1.2.1        digest_0.6.35           stringi_1.8.4           reshape2_1.4.4          splines_4.4.0          
# [43] cowplot_1.1.3           fastmap_1.1.1           grid_4.4.0              colorspace_2.1-0        cli_3.6.2               magrittr_2.0.3         
# [49] utf8_1.2.4              withr_3.0.0             scales_1.3.0            UCSC.utils_1.0.0        bit64_4.0.5             timechange_0.3.0       
# [55] XVector_0.44.0          httr_1.4.7              bit_4.0.5               qvalue_2.36.0           hms_1.1.3               png_0.1-8              
# [61] memoise_2.0.1           IRanges_2.38.0          rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2              
# [67] BiocGenerics_0.50.0     rstudioapi_0.16.0       jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1                fs_1.6.4               
# [73] zlibbioc_1.50.0

setwd("~/irf4_mm_trn_code/figure_1/")

'%!in%' <- function(x,y)!('%in%'(x,y))

metadata <-  read.csv("../public_data/Model.csv")
Chronos <-  read.csv("../public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")

tf_list <- read.delim("~/data/TFlist_NMid_hg19.txt", header=FALSE)

Chronos <- Chronos %>% gather(GENE,Dep,-X)

Chronos$MM <- 0 
Chronos$MM[Chronos$X %in% metadata$ModelID] <- 1

genes <- Chronos %>% filter(Dep < -1,MM == 1)
Chronos <- Chronos %>% filter(GENE %in% genes$GENE)

temp_table <- Chronos %>% filter(Dep < -1 ) 
temp_table <- table(temp_table$GENE)

total_cells <- unique(Chronos$X) %>% length()
temp_table <- temp_table %>% as.data.frame()
temp_table$percent_dep <- temp_table$Freq/total_cells

### Identifying all potential context dependency genes (not pan essential, and essential in > 1 MM cell line)

temp_table <- temp_table %>% filter(percent_dep < 0.90) %>% filter(Freq > 1)

Chronos <- Chronos %>% filter(GENE %in% temp_table$Var1)

genes <- Chronos %>% filter(GENE %in% temp_table$Var1) %>% dplyr::select(GENE) %>% unique()

rm(temp_table)

#### creating a list of genes to loop through, and calculating specific statistics 

gene_list <- genes %>%
  distinct(GENE) %>%
  mutate(wilcox_p = 1, delta = 0, percent_dep = 0,median_chronos = 0)

for (gene in gene_list$GENE){
  x <- Chronos %>% filter(GENE == gene, MM == 1) %>% pull(Dep)
  median_chro <- median(x)
  y <- Chronos %>% filter(GENE == gene, MM == 0) %>% pull(Dep)
  
  ### Wilcox test of MM dependencies scores vs every other cell line
  wilcox_p <- wilcox.test(x, y, "less")$p.value
  
  ### Delta of median dependency scores between condition
  delta <- median(y) - median(x)
  
  ### Percent of cell lines that are dependencies
  percent_n <- Chronos %>% filter(GENE == gene, Dep < -1) %>% nrow()
  percent_d <- Chronos %>% filter(GENE == gene) %>% nrow()
  percent_dep <- percent_n / percent_d
  
  gene_list[gene_list$GENE == gene,1:5] <- c(gene,wilcox_p, delta, percent_dep,median_chro)
}

gene_list$wilcox_p<- as.numeric(gene_list$wilcox_p)
gene_list$delta<- as.numeric(gene_list$delta)
gene_list$percent_dep<- as.numeric(gene_list$percent_dep)
gene_list$median_chronos<- as.numeric(gene_list$median_chronos)

### converting to log p values for plotting

gene_list$log_p <- -log((gene_list$wilcox_p))

### generating a line for relative ideal target ID, has to have sufficient impact on fitness (median chronos score < -0.5) and be context specific

dat <- data.frame(median_chronos = gene_list$median_chronos,
                  y=(max(-log(gene_list$wilcox_p))/(pi^-(gene_list$median_chronos-0.5))) + 1)

gene_list <- merge(gene_list,dat)

###original color scheme, removed
gene_list$color <- "dodgerblue4"
gene_list$color[gene_list$log_p > gene_list$y] <- "seagreen3"
keep <- gene_list %>% filter(color == "seagreen3")
keep$Type <- "MM"

gene_list$Ideal_Target <- "No"
gene_list$Ideal_Target[gene_list$log_p > gene_list$y] <- "Yes"

gene_list$median_chronos <- round(gene_list$median_chronos,3)
gene_list$percent_dep <- round(gene_list$percent_dep,3)
gene_list$log_p <- round(gene_list$log_p,3)
gene_list$delta <- round(gene_list$delta,3)

dat$median_chronos <- round(dat$median_chronos,3)
dat$y <- round(dat$y,3)

gene_list$delta_size <- abs(gene_list$delta)
gene_list$delta_size[gene_list$delta_size > 0.3] <- 0.3

### MM TRN TFs + CREBBP

SE_TFs <- c("ATF4","PRDM10","IRF4","PRDM1","MEF2C","MYC","POU2F1","TCF3","ZNF217","CREBBP")

gene_list$shape_style <- 22
gene_list$shape_style[gene_list$GENE %in% tf_list$V2] <- 23
gene_list$shape_style[gene_list$GENE %in% SE_TFs ] <- 21

gene_list$Gene_Type <- "Other"
gene_list$Gene_Type[gene_list$GENE %in% tf_list$V2] <- "TF"
gene_list$Gene_Type[gene_list$GENE %in% SE_TFs] <- "SE TF"

t <- gene_list
gene_list <- gene_list %>% filter(percent_dep < 0.9)

### arbitrary label, attempting to get top left of plot of TRN genes

gene_list$label <- NULL
gene_list$label[(gene_list$log_p > 20 & gene_list$median_chronos < -0.75)] <- gene_list$GENE[(gene_list$log_p > 20 & gene_list$median_chronos < -0.75)]

node_table <- read.csv("../figure_1/mm_trn/trn_node_table_step4.csv")

gene_list$TRN_network <- "No"
gene_list$TRN_network[gene_list$GENE %in% node_table$GENE] <- "Yes"

gene_list$label <- NULL
gene_list$label[(gene_list$log_p > 20 & gene_list$median_chronos < -0.75)] <- gene_list$GENE[(gene_list$log_p > 20 & gene_list$median_chronos < -0.75)]
gene_list$label[gene_list$TRN_network == "No"] <- ""

gene_list$label[gene_list$GENE == "CREBBP"] <- "CREBBP"

gene_list$TRN_network_alpha<- 0.3
gene_list$TRN_network_alpha[gene_list$GENE %in% node_table$GENE] <- 0.8

p <- ggplot() + 
  geom_line(aes(dat$median_chronos,dat$y), color = "black") +
  geom_point(data = gene_list, aes(
    x = median_chronos ,
    y = -log(wilcox_p),
    shape = Gene_Type,
    color = TRN_network,
    fill = Ideal_Target,
  ),
  alpha =  gene_list$TRN_network_alpha,
  size = 10) + 
  ggrepel::geom_label_repel(data = gene_list, aes(label = label,
                                                  x = median_chronos ,
                                                  y = -log(wilcox_p)),position = "dodge",max.overlaps = 1000) +
  ylim(c(0,max(-log((gene_list$wilcox_p)))+2)) +
  xlim(c(min(gene_list$median_chronos) - 0.2, 0.5)) +
  xlab("Median Chronos") +  
  scale_shape_manual(values = c(22,21,23)) +
  theme_classic(base_size = 20) +
  ylab("Context Specificity") + theme(legend.position = "bottom") +
  guides(guide_legend(nrow = 2))  +
  guides(color = guide_legend(nrow = 2,override.aes = list(size=8)), 
         shape = guide_legend(nrow = 2,override.aes = list(size=8))) +
  scale_color_manual(values = c("gray","black")) +
  scale_fill_manual(values = c("gray","purple")) 
p

pdf('./fig_s1_context_dep_nodes.pdf',height = 9,width = 9)
p
dev.off()


