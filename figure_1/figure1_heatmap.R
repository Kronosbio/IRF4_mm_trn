rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)
library(readxl)

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
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] org.Hs.eg.db_3.19.1   AnnotationDbi_1.66.0  IRanges_2.38.0        S4Vectors_0.42.0      Biobase_2.64.0        BiocGenerics_0.50.0  
# [7] lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.5          
# [13] tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0       ComplexHeatmap_2.20.0 circlize_0.4.16      
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        HDO.db_0.99.1           blob_1.2.4              Biostrings_2.72.0       fastmap_1.1.1           digest_0.6.35          
# [7] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.6           KEGGREST_1.44.0         RSQLite_2.3.6           magrittr_2.0.3         
# [13] compiler_4.4.0          rlang_1.1.3             tools_4.4.0             utf8_1.2.4              data.table_1.15.4       bit_4.0.5              
# [19] plyr_1.8.9              RColorBrewer_1.1-3      BiocParallel_1.38.0     withr_3.0.0             fansi_1.0.6             GOSemSim_2.30.0        
# [25] colorspace_2.1-0        GO.db_3.19.1            scales_1.3.0            iterators_1.0.14        cli_3.6.2               crayon_1.5.2           
# [31] generics_0.1.3          rstudioapi_0.16.0       tzdb_0.4.0              httr_1.4.7              reshape2_1.4.4          rjson_0.2.21           
# [37] DBI_1.2.2               qvalue_2.36.0           cachem_1.0.8            DOSE_3.30.1             zlibbioc_1.50.0         splines_4.4.0          
# [43] parallel_4.4.0          XVector_0.44.0          matrixStats_1.3.0       yulab.utils_0.1.4       vctrs_0.6.5             Matrix_1.7-0           
# [49] jsonlite_1.8.8          hms_1.1.3               GetoptLong_1.0.5        bit64_4.0.5             clue_0.3-65             foreach_1.5.2          
# [55] glue_1.7.0              codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.5            shape_1.4.6.1          
# [61] GenomeInfoDb_1.40.0     UCSC.utils_1.0.0        munsell_0.5.1           pillar_1.9.0            fgsea_1.30.0            GenomeInfoDbData_1.2.12
# [67] R6_2.5.1                doParallel_1.0.17       lattice_0.22-6          png_0.1-8               memoise_2.0.1           Rcpp_1.0.12            
# [73] fastmatch_1.1-4         fs_1.6.4                pkgconfig_2.0.3         GlobalOptions_0.1.2 

genematch <- function(dat){
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL') %>% unlist()
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

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x,y)!('%in%'(x,y))

### Reading in APMS data

IRF4_MS <- read_excel("./Table S2.xlsx",sheet = 1)
p300_MS <- read_excel("./Table S2.xlsx",sheet = 2)

p300_MS <- p300_MS %>% dplyr::select(`Gene Symbol`) %>% filter(!is.na(`Gene Symbol`)) %>% unique()
p300_MS$GENE <- p300_MS$`Gene Symbol`

p300_MS <- genematch(p300_MS)

IRF4_MS <- IRF4_MS %>% dplyr::select(`Gene Symbol`) %>% filter(!is.na(`Gene Symbol`)) %>% unique()
IRF4_MS$GENE <- IRF4_MS$`Gene Symbol`

IRF4_MS <- genematch(IRF4_MS)

### Reading in ChIP peaks - this is IRF4/p300 post ROSE2, so the ChIP peaks used are 'SE' peaks within IRF4/p300 ChIP datasets.

p300_chip1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/33_0EX9_0229Kronos_MM1S-No-Treatment-1_p300_hs_i49_peaks_SuperEnhancers_ENHANCER_TO_GENE_with_MYC.txt")
p300_chip2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/34_0EXA_0229Kronos_MM1S-No-Treatment-2_p300_hs_i50_peaks_SuperEnhancers_ENHANCER_TO_GENE_with_MYC.txt")

p300_chip <- c(p300_chip1$CLOSEST_GENE,p300_chip2$CLOSEST_GENE) 

p300_chip  <- table(p300_chip) %>% as.data.frame()

p300_chip <- p300_chip %>% filter(Freq > 1)
p300_chip$GENE <- p300_chip$p300_chip

p300_chip <- genematch(p300_chip)
p300_chip$GENE[p300_chip$GENE == "DUSP22"] <- "IRF4"

irf4_chip1 <- read.delim("~/irf4_mm_trn_code/ChIP_data/31_0EX7_0229Kronos_MM1S-No-Treatment-1_IRF4_hs_i45_peaks_SuperEnhancers_ENHANCER_TO_GENE_with_MYC.txt")
irf4_chip2 <- read.delim("~/irf4_mm_trn_code/ChIP_data/32_0EX8_0229Kronos_MM1S-No-Treatment-2_IRF4_hs_i48_peaks_SuperEnhancers_ENHANCER_TO_GENE_with_MYC.txt")

irf4_chip <- c(irf4_chip1$CLOSEST_GENE,irf4_chip2$CLOSEST_GENE) 

irf4_chip  <- table(irf4_chip) %>% as.data.frame()

irf4_chip <- irf4_chip %>% filter(Freq > 1)
irf4_chip$GENE <- irf4_chip$irf4_chip

irf4_chip <- genematch(irf4_chip)
irf4_chip$GENE[irf4_chip$GENE == "DUSP22"] <- "IRF4"

### reading in TRN table 4

node_table <- read.csv("~/irf4_mm_trn_code/figure_1/mm_trn/trn_node_table_step4.csv")
edge_table <- read.csv("~/irf4_mm_trn_code/figure_1/mm_trn/trn_crc_ppi_dep_step4.csv")

crc_nodes <- node_table %>% filter(type == "SE_TF") %>% pull(GENE)

node_table <- node_table %>% filter(type == "SE_TF") %>% dplyr::select(GENE)

### other TFs and interactors were chosen by hand for representation

other_tfs <- c("IKZF1", "IKZF3", "MEF2D", "RUNX1","ARID1A","ARID1B","MEF2A", "NFKB1")
interactors <- c("DPF2","SMARCC1","SMARCA4","SMARCD1","SMARCD2")

t <- data_frame("GENE" = c(other_tfs,interactors))

node_table <- rbind(node_table,t)

rm(t)


### wilcox test used for context specificity of MM, method = "less" for heatmap itself. 

wilcox_test_function <- function(df) {
  result <- wilcox.test(Exp ~ MM, data = df, method = "less")
  return(data.frame(GENE = df$GENE[1], p_value = result$p.value))
}

### chi-square test used for co-depdendency relationships 

perform_chi_square_test <- function(gene1, gene2, data) {
  data <- data %>% dplyr::select(X,GENE,Dep) %>% filter(GENE %in% c(gene1,gene2)) %>% spread(GENE,Dep)
  data$X <- NULL
  contingency_table <- as.data.frame.matrix(table(data))
  if(dim(contingency_table)[1] == 1 & dim(contingency_table)[2] == 1){
    return(1)
  }
  else{
    return(chisq.test(contingency_table, simulate.p.value=TRUE)$p.value)
  }
}

### reading in STRING PPI data

PPI_info <- read_table("~/irf4_mm_trn_code/public_data/9606.protein.info.v11.5.txt")
PPI_links <- read_table("~/irf4_mm_trn_code/public_data/9606.protein.physical.links.v11.5.txt")

colnames(PPI_info)[2] <- "GENE"

PPI_info <- genematch(PPI_info)

colnames(PPI_links)[1] <- "string_protein_id"
PPI_links <- merge(PPI_links,PPI_info)
PPI_links$string_protein_id <- NULL
colnames(PPI_links)[1] <- "string_protein_id"
colnames(PPI_links)[3] <- "GENE1"
PPI_links <- merge(PPI_links,PPI_info)
PPI_links$string_protein_id <- NULL
colnames(PPI_links)[3] <- "GENE2"

String_ppi <- PPI_links %>% filter(GENE1 %in% c("IRF4","EP300"))
String_ppi$GENE <- String_ppi$GENE2
String_ppi <- genematch(String_ppi)


### DepMap Chronos scores used for context specificity metric, and co-dependency rows. 

metadata <-  read.csv("~/irf4_mm_trn_code/public_data/Model.csv")
Chronos <-  read.csv("~/irf4_mm_trn_code/public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")
Chronos <- Chronos %>% gather(GENE,Exp,-X)

Chronos$MM <- 0 
Chronos$MM[Chronos$X %in% metadata$ModelID] <- 1

Chronos <- Chronos %>% filter(GENE %in% c("EP300",node_table$GENE))

result_df <- Chronos %>%
  group_by(GENE) %>%
  do(wilcox_test_function(.))

node_table <- merge(node_table,result_df)

node_table$p300_PPI_edge <- 0
node_table$IRF4_PPI_edge <- 0
node_table$p300_coDep_edge <- 0
node_table$IRF4_coDep_edge <- 0
node_table$p300_chip <- 0
node_table$IRF4_chip <- 0
node_table$p300_MS <- 0
node_table$IRF4_MS <- 0

edge_table <- edge_table %>% filter(To %in% c("EP300",crc_nodes,other_tfs,interactors)) %>% filter(From %in% c("EP300",crc_nodes,other_tfs,interactors))

genes_list <- c(node_table$GENE,"EP300")
gene1s <- c()
gene2s <- c()
p_vals <- c()

Chronos$Dep <- 0 
Chronos$Dep[Chronos$Exp < -1] <- 1
data <- Chronos %>% filter(GENE %in% c("EP300",node_table$GENE)) %>% group_by(GENE) %>% mutate(sum_dep = sum(Dep))

# Loop through all pairs of genes and perform chi-square test
for (i in 1:(length(genes_list) - 1)) {
  for (j in (i + 1):length(genes_list)) {
    gene1 <- genes_list[i]
    
    gene2 <- genes_list[j]
    
    p_value <- perform_chi_square_test(gene1, gene2, data)
    
    gene1s <- c(gene1s,gene1)
    gene2s <- c(gene2s,gene2)
    p_vals <- c(p_vals,p_value)
    
  }
}

codep <- data.frame(gene1s,gene2s,p_vals)

### putting all the data together to generate heatmap in figure 1

for(gene in node_table$GENE){
  
  t_string_ppi <- String_ppi %>% filter(GENE == gene)
  t_MS_p300 <- p300_MS %>% filter(GENE == gene)
  t_MS_IRF4 <- IRF4_MS %>% filter(GENE == gene)
  
  t_ChIP_p300 <- p300_chip %>% filter(GENE == gene)
  t_ChIP_IRF4 <- irf4_chip %>% filter(GENE == gene)
  
  pvals  <- codep %>% filter(gene1s == gene | gene2s  == gene)
  col <- grep("EP300",pvals)
  pval_heatmap <- pvals$p_vals[which(pvals[col] == "EP300")]
  node_table$p300_coDep_edge[node_table$GENE == gene] <- pval_heatmap
  
  ##code breaks with IRF4, but IRF4 co-dependence is omitted in heatmap
  
  col <- grep("IRF4",pvals)
  pval_heatmap <- pvals$p_vals[which(pvals[col] == "IRF4")]
  node_table$IRF4_coDep_edge[node_table$GENE == gene] <- pval_heatmap

  if(nrow(t_string_ppi) > 0){
    if("IRF4" %in% t_string_ppi$GENE1){
      node_table$IRF4_PPI_edge[node_table$GENE == gene] <-1
    }
    if("EP300" %in% t_string_ppi$GENE1){
      node_table$p300_PPI_edge[node_table$GENE == gene] <-1
    }
  }
  
  if(nrow(t_MS_p300) > 0){
    node_table$p300_MS[node_table$GENE == gene] <-1
  }
  if(nrow(t_ChIP_p300) > 0){
    node_table$p300_chip[node_table$GENE == gene] <-1
  }
  if(nrow(t_MS_IRF4) > 0){
    node_table$IRF4_MS[node_table$GENE == gene] <-1
  }
  if(nrow(t_ChIP_IRF4) > 0){
    node_table$IRF4_chip[node_table$GENE == gene] <-1
  }
}

### setting IRF4-IRF4 co-dep to null
node_table$IRF4_coDep_edge[node_table$GENE == "IRF4"] <- NA

node_table <- node_table %>% column_to_rownames("GENE")
crc_table <- node_table[crc_nodes,]
other_tfs_table <- node_table[other_tfs,]
interactors_table <- node_table[interactors,]

crc_table <- crc_table[order(crc_table$p_value,decreasing = TRUE), ]
other_tfs_table <- other_tfs_table[order(other_tfs_table$p_value,decreasing = TRUE), ]
interactors_table <- interactors_table[order(interactors_table$p_value,decreasing = TRUE), ]

node_table <- rbind(crc_table,other_tfs_table)
node_table <- rbind(node_table,interactors_table)

col_fun_c = colorRamp2(c(0, 30), c( "white", "purple"))
ha = HeatmapAnnotation(context_dep = -log(node_table$p_value), col = list(context_dep = col_fun_c))



### The co-depend color scheme needed to be unique to the rest of the heatmap. 
### The first heatmap is the one used in the figure, the second was just for the co-depend heatmap rows which was used to replace the co-dependency heatmap in the figure.

node_table$p_value <- NULL
col_fun = colorRamp2(c(0, 1), c( "white", "black"))

pdf("./fig1b_tf_heatmap.pdf",height = 4,width = 8)
Heatmap(t(node_table),cluster_rows =  FALSE,width = ncol(t(node_table))*unit(5, "mm"), 
        height = nrow(t(node_table))*unit(5, "mm"),col = col_fun,cluster_columns = FALSE,top_annotation = ha)
dev.off()


col_fun = colorRamp2(c(0,0.01, 1), c( "black","gray", "white"))

pdf("./fig1b_tf_heatmap_copend_color.pdf",height = 4,width = 8)
Heatmap(t(node_table),cluster_rows =  FALSE,width = ncol(t(node_table))*unit(5, "mm"), 
        height = nrow(t(node_table))*unit(5, "mm"),col = col_fun,cluster_columns = FALSE,top_annotation = ha)
dev.off()
