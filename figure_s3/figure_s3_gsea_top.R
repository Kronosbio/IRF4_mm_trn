rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)

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
#   [1] clusterProfiler_4.12.0 msigdbr_7.5.1          lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1         
# [6] dplyr_1.1.4            purrr_1.0.2            readr_2.1.5            tidyr_1.3.1            tibble_3.2.1          
# [11] ggplot2_3.5.1          tidyverse_2.0.0        org.Hs.eg.db_3.19.1    AnnotationDbi_1.66.0   IRanges_2.38.0        
# [16] S4Vectors_0.42.0       Biobase_2.64.0         BiocGenerics_0.50.0   
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               gson_0.1.0              shadowtext_0.1.3        gridExtra_2.3           rlang_1.1.3            
# [6] magrittr_2.0.3          DOSE_3.30.1             compiler_4.4.0          RSQLite_2.3.6           png_0.1-8              
# [11] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         crayon_1.5.2            fastmap_1.1.1          
# [16] XVector_0.44.0          ggraph_2.2.1            utf8_1.2.4              HDO.db_0.99.1           tzdb_0.4.0             
# [21] enrichplot_1.24.0       UCSC.utils_1.0.0        bit_4.0.5               zlibbioc_1.50.0         cachem_1.0.8           
# [26] aplot_0.2.2             GenomeInfoDb_1.40.0     jsonlite_1.8.8          blob_1.2.4              BiocParallel_1.38.0    
# [31] tweenr_2.0.3            parallel_4.4.0          R6_2.5.1                stringi_1.8.4           RColorBrewer_1.1-3     
# [36] GOSemSim_2.30.0         Rcpp_1.0.12             Matrix_1.7-0            splines_4.4.0           igraph_2.0.3           
# [41] timechange_0.3.0        tidyselect_1.2.1        qvalue_2.36.0           rstudioapi_0.16.0       viridis_0.6.5          
# [46] codetools_0.2-20        lattice_0.22-6          plyr_1.8.9              treeio_1.28.0           withr_3.0.0            
# [51] KEGGREST_1.44.0         gridGraphics_0.5-1      scatterpie_0.2.3        polyclip_1.10-6         Biostrings_2.72.0      
# [56] ggtree_3.12.0           pillar_1.9.0            ggfun_0.1.4             generics_0.1.3          hms_1.1.3              
# [61] tidytree_0.4.6          munsell_0.5.1           scales_1.3.0            glue_1.7.0              lazyeval_0.2.2         
# [66] tools_4.4.0             data.table_1.15.4       fgsea_1.30.0            babelgene_22.9          fs_1.6.4               
# [71] graphlayouts_1.1.1      fastmatch_1.1-4         tidygraph_1.3.1         cowplot_1.1.3           grid_4.4.0             
# [76] ape_5.8                 colorspace_2.1-0        nlme_3.1-164            patchwork_1.2.0         GenomeInfoDbData_1.2.12
# [81] ggforce_0.4.2           cli_3.6.2               fansi_1.0.6             viridisLite_0.4.2       gtable_0.3.5           
# [86] yulab.utils_0.1.4       digest_0.6.35           ggrepel_0.9.5           ggplotify_0.1.2         farver_2.1.1           
# [91] memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7              GO.db_3.19.1            bit64_4.0.5            
# [96] MASS_7.3-60.2     

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

setwd("~/irf4_mm_trn_code/figure_s3")

##############

gene_list_boxplots <- read.csv("../figure_3/data/gene_list_boxplots.csv")

top_genes <- gene_list_boxplots %>% dplyr::select(top_1000)
bot_genes <- gene_list_boxplots %>% dplyr::select(bot_1000)

colnames(top_genes) <- "GENE"
colnames(bot_genes) <- "GENE"

enhancerPromoter_all <- read.csv("../ChIP_data/enhancerPromoter_DMSO_with_MYC.csv") %>% unique()

ep_dat <- enhancerPromoter_all %>% dplyr::select(GENE_ID) %>% unique()

colnames(ep_dat) <- "GENE"

### for over representation analysis, followed example found at :https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html

### matching gene names to ENTREZ ID

entrez_vector <- mapIds(
  org.Hs.eg.db,
  keys = ep_dat$GENE,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

gene_key_df <- data.frame(
  ensembl_id = names(entrez_vector),
  entrez_id = entrez_vector,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(entrez_id))


top_genes <- gene_key_df %>% filter(ensembl_id %in% top_genes$GENE)
hs_msigdb_df <- msigdbr(species = "Homo sapiens")

### pulling out wikipathways

hs_msigdb_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2",
    gs_subcat == "CP:WIKIPATHWAYS"
  )

background <- gene_key_df$entrez_id %>% unique()

ora_results <- enricher(TERM2GENE = dplyr::select(
  hs_msigdb_df ,
  gs_name,
  human_entrez_gene),
  gene = top_genes$entrez_id,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe = background)

result_df <- data.frame(ora_results@result)

write_csv(result_df,"./GSEA_results_figS3.csv")

### filtering for adjusted pvalue < 0.05

result_df <- result_df %>%
  dplyr::filter(p.adjust < 0.05)

### ordering by percentage of elgible pathway covered.

result_df$Bg_count <-sapply(strsplit(result_df$BgRatio, "/"), function(x) as.numeric(x[1]))
result_df$rep <- result_df$Count/result_df$Bg_count

pdf("./fig_S3_IRF4_p300_ora_gsea.pdf",height = 5,width = 15)
ggplot(result_df, aes(x = rep, y = reorder(ID,-(-log(rep))), size = -log(p.adjust))) + 
  geom_point(size = 3*-log(result_df$p.adjust)) + 
  xlab("Percentage of Pathway Represented") + ylab("") + theme_classic(base_size = 15) 
dev.off()
