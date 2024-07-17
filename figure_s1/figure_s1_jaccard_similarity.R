

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
gc()

set.seed(1)

library(tidyverse)
library(igraph)

setwd("~/irf4_mm_trn_code/figure_s1")

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] igraph_2.0.3    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
# [10] ggplot2_3.5.1   tidyverse_2.0.0
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
# [55] httr_1.4.7              bit_4.0.5               qvalue_2.36.0           hms_1.1.3               png_0.1-8               memoise_2.0.1          
# [61] IRanges_2.38.0          rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0              DBI_1.2.2               BiocGenerics_0.50.0    
# [67] rstudioapi_0.16.0       jsonlite_1.8.8          plyr_1.8.9              R6_2.5.1                fs_1.6.4                zlibbioc_1.50.0    

'%!in%' <- function(x,y)!('%in%'(x,y))

node_table <- read.csv("../figure_1/mm_trn/trn_node_table_step4.csv")
edge_table <- read.csv("../figure_1/mm_trn/trn_crc_ppi_dep_step4.csv")

g <- graph_from_data_frame(edge_table, directed = FALSE)

# Function to calculate Jaccard similarity between two nodes
jaccard_similarity <- function(graph, node1, node2) {
  neighbors1 <- neighbors(graph, node1)
  neighbors2 <- neighbors(graph, node2)
  
  intersection_size <- length(intersect(neighbors1, neighbors2))
  union_size <- length(union(neighbors1, neighbors2))
  
  if (union_size == 0) {
    return(0)  
  } else {
    return(intersection_size / union_size)
  }
}

# Calculate Jaccard similarity for all pairs of nodes in the graph
# Get node names
node_names <- V(g)$name

# Calculate Jaccard similarity for all pairs of nodes in the graph
num_nodes <- vcount(g)
jaccard_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes, dimnames = list(node_names, node_names))

for (i in 1:num_nodes) {
  for (j in 1:num_nodes) {
    jaccard_matrix[i, j] <- jaccard_similarity(g, node_names[i], node_names[j])
  }
}

# Print the Jaccard similarity matrix
print(jaccard_matrix)

crc_nodes <- node_table %>% filter(type == "SE_TF") %>% pull(GENE)

non_crc_nodes <- colnames(jaccard_matrix) [colnames(jaccard_matrix) %!in% crc_nodes]

IRF4_Jaccard <- jaccard_matrix[non_crc_nodes,"IRF4"] %>% as.data.frame()
colnames(IRF4_Jaccard) <- "jaccard"
IRF4_Jaccard$GENE <- rownames(IRF4_Jaccard)
a <- IRF4_Jaccard %>% #filter(jaccard > 0.15) %>%
  ggplot(aes(x = reorder(GENE,-jaccard),y = jaccard)) + ggtitle("Most similar nodes to IRF4") +
  geom_point(size = 4) + scale_x_discrete(guide = guide_axis(angle = 90)) + xlab("PPI Co-Factor Nodes") + theme_classic(base_size = 13) + ylab("Jaccard Similarity")

a

pdf("./fig_s1_jaccard_simarility_individual.pdf",height = 4,width = 6)
a
dev.off()
