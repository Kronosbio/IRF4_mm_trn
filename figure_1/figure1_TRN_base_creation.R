##################

### Step 1
rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)

# sessionInfo()
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
#   [1] KEGGREST_1.44.0         fastmatch_1.1-4         gtable_0.3.5            GOSemSim_2.30.0         lattice_0.22-6          tzdb_0.4.0             
# [7] vctrs_0.6.5             tools_4.4.0             generics_0.1.3          yulab.utils_0.1.4       parallel_4.4.0          fansi_1.0.6            
# [13] RSQLite_2.3.6           blob_1.2.4              pkgconfig_2.0.3         Matrix_1.7-0            data.table_1.15.4       lifecycle_1.0.4        
# [19] GenomeInfoDbData_1.2.12 HDO.db_0.99.1           compiler_4.4.0          Biostrings_2.72.0       munsell_0.5.1           fgsea_1.30.0           
# [25] codetools_0.2-20        DOSE_3.30.1             GenomeInfoDb_1.40.0     pillar_1.9.0            crayon_1.5.2            GO.db_3.19.1           
# [31] BiocParallel_1.38.0     cachem_1.0.8            tidyselect_1.2.1        digest_0.6.35           stringi_1.8.4           reshape2_1.4.4         
# [37] splines_4.4.0           cowplot_1.1.3           fastmap_1.1.1           grid_4.4.0              colorspace_2.1-0        cli_3.6.2              
# [43] magrittr_2.0.3          utf8_1.2.4              withr_3.0.0             scales_1.3.0            UCSC.utils_1.0.0        bit64_4.0.5            
# [49] timechange_0.3.0        XVector_0.44.0          httr_1.4.7              bit_4.0.5               qvalue_2.36.0           hms_1.1.3              
# [55] png_0.1-8               memoise_2.0.1           rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2              
# [61] vroom_1.6.5             rstudioapi_0.16.0       jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1                fs_1.6.4               
# [67] zlibbioc_1.50.0  

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

genematch <- function(dat) {
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL') %>% unlist()
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

### STEP 1 CRC Pruning created in CRC_edge_pruning.

crc_total <- read.csv("./mm_trn/crc_total_public_chip.csv")

crc_total$X <- NULL

patients <- c(
  "MM11.H3K27ac_EDGE_LIST.txt",
  "MM24.H3K27ac_EDGE_LIST.txt",
  "MM1.H3K27ac_EDGE_LIST.txt",
  "MM2.H3K27ac_EDGE_LIST.txt",
  "MM3.H3K27ac_EDGE_LIST.txt",
  "MM1_jia_EDGE_LIST.txt",
  "MM4_jia_EDGE_LIST.txt",
  "MM8.H3K27ac_EDGE_LIST.txt",
  "MM5.H3K27ac_EDGE_LIST.txt",
  "MM7_jia_EDGE_LIST.txt",
  "MM2_jia_EDGE_LIST.txt",
  "MM6_jia_EDGE_LIST.txt",
  "MM3_jia_EDGE_LIST.txt",
  "MM10_jia_EDGE_LIST.txt",
  "MM5_jia_EDGE_LIST.txt",
  "MM12.H3K27ac_EDGE_LIST.txt",
  "MM8_jia_EDGE_LIST.txt",
  "MM25.H3K27ac_EDGE_LIST.txt",
  "MM14.H3K27ac_EDGE_LIST.txt",
  "MM9_jia_EDGE_LIST.txt"
)

t <- crc_total %>% filter(file %in% patients)

t$comb <- paste(t$From, t$To)
crc_tab <- table(t$comb)
crc_tab <- as.data.frame(crc_tab)

### filtering for edges observed in at least 5 patient samples

crc_tab <- crc_tab %>% filter(Freq > 5)

t_patients <- t %>% filter(comb %in% crc_tab$Var1) %>% dplyr::select(-file) %>% unique()

### aggregating cell line samples. Edges for repeated samples have to be observed in both sample replicates. 

t <- crc_total %>% filter(file %!in% patients)

t$comb <- paste(t$From, t$To)

U266 <- t %>% filter(file %in% c("U266.H3K27ac_EDGE_LIST.txt", "U266_H3K27ac_EDGE_LIST.txt"))

U266 <- table(U266$comb)[table(U266$comb) > 1]
U266_comb  <- t %>% filter(file %in% c("U266_H3K27ac_EDGE_LIST.txt"))  %>% filter(comb %in% names(U266))

t <- t %>% filter(file %!in% c("U266.H3K27ac_EDGE_LIST.txt", "U266_H3K27ac_EDGE_LIST.txt"))
t <- rbind(t, U266_comb)

rm(U266_comb, U266)

RPMI8226 <- t %>% filter(file %in% c(
  "RPMI8226.H3K27ac_EDGE_LIST.txt",
  "RPMI8226_H3K27ac_EDGE_LIST.txt"
))

RPMI8226 <- table(RPMI8226$comb)[table(RPMI8226$comb) > 1]
RPMI8226_comb  <- t %>% filter(file %in% c("RPMI8226_H3K27ac_EDGE_LIST.txt"))  %>% filter(comb %in% names(RPMI8226))

t <- t %>% filter(file %!in% c(
  "RPMI8226.H3K27ac_EDGE_LIST.txt",
  "RPMI8226_H3K27ac_EDGE_LIST.txt"
))
t <- rbind(t, RPMI8226_comb)
rm(RPMI8226_comb, RPMI8226)

KMS12 <- t %>% filter(file %in% c(
  "KMS-12-BM.H3K27ac_EDGE_LIST.txt",
  "KMS12_H3K27ac_EDGE_LIST.txt"
))

KMS12 <- table(KMS12$comb)[table(KMS12$comb) > 1]
KMS12_comb  <- t %>% filter(file %in% c("KMS12_H3K27ac_EDGE_LIST.txt"))  %>% filter(comb %in% names(KMS12))

t <- t %>% filter(file %!in% c(
  "KMS-12-BM.H3K27ac_EDGE_LIST.txt",
  "KMS12_H3K27ac_EDGE_LIST.txt"
))
t <- rbind(t, KMS12_comb)
rm(KMS12_comb, KMS12)

H929 <- t %>% filter(file %in% c(
  "NCI-H929.H3K27ac_EDGE_LIST.txt",
  "H929_H3K27ac_EDGE_LIST.txt"
))

H929 <- table(H929$comb)[table(H929$comb) > 1]
H929_comb  <- t %>% filter(file %in% c("H929_H3K27ac_EDGE_LIST.txt"))  %>% filter(comb %in% names(H929))

t <- t %>% filter(file %!in% c(
  "NCI-H929.H3K27ac_EDGE_LIST.txt",
  "H929_H3K27ac_EDGE_LIST.txt"
))
t <- rbind(t, H929_comb)
rm(H929_comb, H929)

t$comb <- paste(t$From, t$To)
crc_tab <- table(t$comb)
crc_tab <- as.data.frame(crc_tab)
crc_tab <- crc_tab %>% filter(Freq > 3)
t_cells <- t %>%
  filter(comb %in% crc_tab$Var1) %>%
  dplyr::select(-file) %>%
  unique()

t_cells <- t_cells %>% filter(comb %in% t_patients$comb) %>% filter(From != To)
network_data <- t_cells
network_data$comb <- NULL

### Sort the gene names for each edge
sorted_network_data <- t(apply(network_data, 1, function(row)
  sort(row)))

### Identify unique edges
unique_rows <- !duplicated(sorted_network_data)

network_data$type = ifelse(unique_rows, "Unique", "Duplicated")
network_data_uni <- network_data %>% filter(type == "Unique")
network_data_dup <- network_data %>% filter(type == "Duplicated")

network_data_uni$comb <- paste(network_data_uni$From, network_data_uni$To)
network_data_dup$comb <- paste(network_data_dup$To, network_data_dup$From)

network_data_uni$type <- "CRC_Unidirectional"
network_data_uni$type[network_data_uni$comb %in% network_data_dup$comb] <- "CRC_Bidirectional"

network_data_uni$comb <- NULL

write_csv(network_data_uni, "./mm_trn/trn_crc_step1.csv")


##################
### Step 2 adding PPI edges to core MM TFs
rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

genematch <- function(dat) {
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL') %>% unlist()
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

trn_MM_high_confidence <- read.csv("./mm_trn/trn_crc_step1.csv")

trn_MM_high_confidence$GENE <- trn_MM_high_confidence$From
trn_MM_high_confidence <- genematch(trn_MM_high_confidence)
trn_MM_high_confidence$From <- trn_MM_high_confidence$GENE
trn_MM_high_confidence$GENE <- trn_MM_high_confidence$To
trn_MM_high_confidence <- genematch(trn_MM_high_confidence)
trn_MM_high_confidence$To <- trn_MM_high_confidence$GENE
trn_MM_high_confidence$GENE <- NULL

### reading in PPI data from STRING v11.5

PPI_info <- read_table("~/irf4_mm_trn_code/public_data/9606.protein.info.v11.5.txt")
PPI_links <- read_table("~/irf4_mm_trn_code/public_data/9606.protein.physical.links.v11.5.txt")

colnames(PPI_info)[2] <- "GENE"

PPI_info <- genematch(PPI_info)

colnames(PPI_links)[1] <- "string_protein_id"
PPI_links <- merge(PPI_links, PPI_info)
PPI_links$string_protein_id <- NULL
colnames(PPI_links)[1] <- "string_protein_id"
colnames(PPI_links)[3] <- "GENE1"
PPI_links <- merge(PPI_links, PPI_info)
PPI_links$string_protein_id <- NULL
colnames(PPI_links)[3] <- "GENE2"

PPI_links_1 <- PPI_links %>% filter(GENE1 %in% trn_MM_high_confidence$From)
PPI_links_2 <- PPI_links %>% filter(GENE1 %in% trn_MM_high_confidence$To)

node_list <- c(
  PPI_links_1$GENE1,
  PPI_links_1$GENE2,
  PPI_links_2$GENE1,
  PPI_links_2$GENE2,
  trn_MM_high_confidence$From,
  trn_MM_high_confidence$To
) %>% unique()

PPI_links <- PPI_links %>% filter(GENE1 %in% node_list) %>% filter(GENE2 %in% node_list)

### Concatenate the columns and identify unique rows

unique_rows <- !duplicated(apply(PPI_links, 1, function(row)
  paste(sort(row), collapse = ',')))

### Filter the data frame for unique PPI
PPI_links <- PPI_links[unique_rows, ]

rm(node_list, PPI_links_1, PPI_links_2)

PPI_links$combined_score <- NULL

colnames(PPI_links) <- c("To", "From")

edge_table <- PPI_links
PPI_links$type <- "PPI"
PPI_links$combined_score <- NULL

TRN_links <- rbind(PPI_links, trn_MM_high_confidence)

write_csv(TRN_links, "./mm_trn/trn_crc_ppi_step2.csv")


##################
### Step 3 adding in DepMap data  

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_1")
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

genematch <- function(dat) {
  genes <- c(dat$GENE) %>% unique()
  V1 <- as.vector(genes) ### Adding entrez id to fold change
  entrez_IDS <- mapIds(org.Hs.eg.db, V1, 'ENTREZID', 'SYMBOL') %>% unlist()
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


### the wilcox test function is for relative context specificity of MM dependencies, 
### wilcox test can be one-sided, but we are using two-sided for simplicity, as the direction of the effect is selected for when selecting for potential context dependencies
### alternative = 'less' is used for subsequent figures
wilcox_test_function <- function(df) {
  result <- wilcox.test(Exp ~ MM, data = df)
  return(data.frame(GENE = df$GENE[1], p_value = result$p.value))
}

TRN_links <- read.csv("./mm_trn/trn_crc_ppi_step2.csv")

metadata <-  read.csv("~/irf4_mm_trn_code/public_data/Model.csv")
Chronos <-  read.csv("~/irf4_mm_trn_code/public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")
Chronos <- Chronos %>% gather(GENE, Exp, -X)

Chronos$MM <- 0
Chronos$MM[Chronos$X %in% metadata$ModelID] <- 1

### Identifying all MM dependencies in DepMap 22Q4

genes <- Chronos %>% filter(Exp < -1, MM == 1)

Chronos <- Chronos %>% filter(GENE %in% genes$GENE)

temp_table <- Chronos %>% filter(Exp < -1)
temp_table <- table(temp_table$GENE)

total_cells <- unique(Chronos$X) %>% length()
temp_table <- temp_table %>% as.data.frame()
temp_table$percent_dep <- temp_table$Freq/total_cells

### Identifying all potential MM context dependencies (can't be common essential and has to be depedenent in > 1 cell line)

temp_table <- temp_table %>% filter(percent_dep < 0.90) %>% filter(Freq > 1)

Chronos <- Chronos %>% filter(GENE %in% temp_table$Var1)

MM_genes <- temp_table

rm(temp_table, genes)

### calculating Wilcox test p-values. this is for relative context specificity of MM dependencies

result_df <- Chronos %>%
  group_by(GENE) %>%
  do(wilcox_test_function(.))

Chronos$Dep <- 0
Chronos$Dep[Chronos$Exp < -1] <- 1
Chronos_dat <- Chronos

result_df$p_value <- -log(result_df$p_value)

Chronos <- Chronos %>% filter(X %in% metadata$ModelID)

Chronos <- genematch(Chronos)
context_dep <- genematch(result_df)

### function to perform chi-square test for independence, this is used for co-dependency edge identification and calculation. 
perform_chi_square_test <- function(gene1, gene2, data) {
  data <- data %>% dplyr::select(X, GENE, Dep) %>% filter(GENE %in% c(gene1, gene2)) %>% spread(GENE, Dep)
  data$X <- NULL
  contingency_table <- as.data.frame.matrix(table(data))
  if (dim(contingency_table)[1] == 1 &
      dim(contingency_table)[2] == 1) {
    return(1)
  }
  else{
    return(chisq.test(contingency_table, simulate.p.value = TRUE)$p.value)
  }
}

node_table <- context_dep
node_table <- node_table %>% filter(GENE %in% MM_genes$Var1)
node_table <- node_table %>% filter(GENE %in% c(TRN_links$To, TRN_links$From))
TRN_links <- TRN_links %>% filter(To %in% node_table$GENE) %>%
  filter(From %in% node_table$GENE)

genes_list <- unique(node_table$GENE)

data <- Chronos_dat %>% filter(GENE %in% node_table$GENE) %>% group_by(GENE) %>% mutate(sum_dep = sum(Dep))

gene1s <- c()
gene2s <- c()
p_vals <- c()

### Loop through all pairs of genes and perform chi-square test to calculate co-dependency edges 
for (i in 1:(length(genes_list) - 1)) {
  for (j in (i + 1):length(genes_list)) {
    gene1 <- genes_list[i]
    
    gene2 <- genes_list[j]
    
    p_value <- perform_chi_square_test(gene1, gene2, data)
    
    gene1s <- c(gene1s, gene1)
    gene2s <- c(gene2s, gene2)
    p_vals <- c(p_vals, p_value)
    
  }
}

codep <- data.frame(gene1s, gene2s, p_vals)

codep$adjust_pval <- p.adjust(as.numeric(codep$p_vals), "fdr")
codep_fil <- codep %>% filter(adjust_pval < 0.01)

codep_fil <- codep_fil %>% dplyr::select(gene1s, gene2s)
colnames(codep_fil) <- c("To", "From")

codep_fil$type <- "codep"

TRN_links <- rbind(TRN_links, codep_fil)

write_csv(node_table, "./mm_trn/trn_node_table_step3.csv")
write_csv(TRN_links, "./mm_trn/trn_crc_ppi_dep_step3.csv")

##################
### Step 4 pruning for edges that are only connected to remaining TFs

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

TRN_links <- read.csv("./mm_trn/trn_crc_ppi_dep_step3.csv")

TFlist_NMid_hg19 <- read.delim("~/data/TFlist_NMid_hg19.txt", header = FALSE)

TRN_links_CRC_core <- TRN_links %>% filter(type %in% c("CRC_Unidirectional", "CRC_Bidirectional")) %>% filter(To %in% TFlist_NMid_hg19$V2) %>% filter(From %in% TFlist_NMid_hg19$V2)

coreTFs <- unique(c(TRN_links_CRC_core$From, TRN_links_CRC_core$To))

TRN_links1 <- TRN_links %>% filter(To %in% coreTFs, type == "PPI")
TRN_links2 <- TRN_links %>% filter(From %in% coreTFs, type == "PPI")

TRN_links_nodes <- c(TRN_links1$From,
                     TRN_links1$To,
                     TRN_links2$From,
                     TRN_links2$To) %>% unique()

TRN_links <- TRN_links %>% filter(To %in% TRN_links_nodes) %>% filter(From %in% TRN_links_nodes)

###step 4 is then taken to cytoscape and made for figure observed in Figure 1

write_csv(TRN_links, "./mm_trn/trn_crc_ppi_dep_step4.csv")

###########

### Network node stats

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_1")

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))

step1 <- read.csv("./mm_trn/trn_crc_step1.csv")
step2 <- read.csv("./mm_trn/trn_crc_ppi_step2.csv")
step3 <- read.csv("./mm_trn/trn_crc_ppi_dep_step3.csv")
step4 <- read.csv("./mm_trn/trn_crc_ppi_dep_step4.csv")

TFlist_NMid_hg19 <- read.delim("~/data/TFlist_NMid_hg19.txt", header = FALSE)

unique(c(step1$From, step1$To)) %in% TFlist_NMid_hg19$V2 %>% table()
### 58 SE TFs, 7 Other

node_table <- unique(c(step1$From, step1$To)) %>% as.data.frame()
colnames(node_table) <- "GENE"

node_table$type <- "Other"
node_table$type[node_table$GENE %in% TFlist_NMid_hg19$V2] <- "SE_TF"

SE_TFS <- node_table %>% filter(type == "SE_TF") %>% pull(GENE)

write_csv(node_table, "./mm_trn/trn_node_table_step1.csv")

unique(c(step2$From, step2$To)) %in% TFlist_NMid_hg19$V2 %>% table()
### 58 SE TFs, 703 total TFs, 4325 nodes

node_table <- unique(c(step2$From, step2$To)) %>% as.data.frame()
colnames(node_table) <- "GENE"

node_table$type <- "Other"
node_table$type[node_table$GENE %in% TFlist_NMid_hg19$V2] <- "TF"
node_table$type[node_table$GENE %in% SE_TFS] <- "SE_TF"

write_csv(node_table, "./mm_trn/trn_node_table_step2.csv")

node_table <- read.csv("./mm_trn/trn_node_table_step3.csv")

unique(c(step3$From, step3$To)) %in% TFlist_NMid_hg19$V2 %>% table()
### 9 SE TFs, 36 total TFs, 503 nodes

node_table$type <- "Other"
node_table$type[node_table$GENE %in% TFlist_NMid_hg19$V2] <- "TF"
node_table$type[node_table$GENE %in% SE_TFS] <- "SE_TF"

write_csv(node_table, "./mm_trn/trn_node_table_step3.csv")

unique(c(step4$From, step4$To)) %in% TFlist_NMid_hg19$V2 %>% table()
### 9 SE TFs, 28 total TFs, 220 nodes

node_table <- read.csv("./mm_trn/trn_node_table_step3.csv")

TRN_links_CRC_core <- step4 %>% filter(type %in% c("CRC_Unidirectional", "CRC_Bidirectional")) %>%
  filter(To %in% TFlist_NMid_hg19$V2) %>% filter(From %in% TFlist_NMid_hg19$V2)

coreTFs <- unique(c(TRN_links_CRC_core$From, TRN_links_CRC_core$To))

node_table <- node_table %>% filter(GENE %in% c(step4$From, step4$To))
write_csv(node_table, "./mm_trn/trn_node_table_step4.csv")
