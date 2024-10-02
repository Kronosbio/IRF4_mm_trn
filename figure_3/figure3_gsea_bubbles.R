rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(org.Hs.eg.db)
library(tidyverse)
library(fgsea)

#### This code is used to generate the GSEA bubbles in Figure 3.

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
#   [1] fgsea_1.30.0         lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
# [7] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0      org.Hs.eg.db_3.19.1 
# [13] AnnotationDbi_1.66.0 IRanges_2.38.0       S4Vectors_0.42.0     Biobase_2.64.0       BiocGenerics_0.50.0 
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               gridExtra_2.3           rlang_1.1.3             magrittr_2.0.3          clue_0.3-65            
# [6] GetoptLong_1.0.5        DOSE_3.30.1             matrixStats_1.3.0       compiler_4.4.0          RSQLite_2.3.6          
# [11] png_0.1-8               vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         shape_1.4.6.1          
# [16] crayon_1.5.2            fastmap_1.1.1           XVector_0.44.0          ggraph_2.2.1            utf8_1.2.4             
# [21] HDO.db_0.99.1           tzdb_0.4.0              UCSC.utils_1.0.0        bit_4.0.5               zlibbioc_1.50.0        
# [26] cachem_1.0.8            GenomeInfoDb_1.40.0     jsonlite_1.8.8          blob_1.2.4              BiocParallel_1.38.0    
# [31] tweenr_2.0.3            parallel_4.4.0          cluster_2.1.6           R6_2.5.1                stringi_1.8.4          
# [36] RColorBrewer_1.1-3      GOSemSim_2.30.0         Rcpp_1.0.12             iterators_1.0.14        timechange_0.3.0       
# [41] Matrix_1.7-0            splines_4.4.0           igraph_2.0.3            tidyselect_1.2.1        qvalue_2.36.0          
# [46] rstudioapi_0.16.0       viridis_0.6.5           doParallel_1.0.17       codetools_0.2-20        lattice_0.22-6         
# [51] plyr_1.8.9              withr_3.0.0             KEGGREST_1.44.0         polyclip_1.10-6         circlize_0.4.16        
# [56] Biostrings_2.72.0       pillar_1.9.0            foreach_1.5.2           generics_0.1.3          hms_1.1.3              
# [61] munsell_0.5.1           scales_1.3.0            glue_1.7.0              tools_4.4.0             data.table_1.15.4      
# [66] fs_1.6.4                graphlayouts_1.1.1      fastmatch_1.1-4         tidygraph_1.3.1         cowplot_1.1.3          
# [71] grid_4.4.0              colorspace_2.1-0        GenomeInfoDbData_1.2.12 ggforce_0.4.2           cli_3.6.2              
# [76] fansi_1.0.6             viridisLite_0.4.2       ComplexHeatmap_2.20.0   gtable_0.3.5            yulab.utils_0.1.4      
# [81] digest_0.6.35           ggrepel_0.9.5           rjson_0.2.21            farver_2.1.1            memoise_2.0.1          
# [86] lifecycle_1.0.4         httr_1.4.7              GlobalOptions_0.1.2     GO.db_3.19.1            bit64_4.0.5            
# [91] MASS_7.3-60.2    

#########

setwd("~/irf4_mm_trn_code/figure_3") 

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

wilcox_test_function <- function(df) {
  result <- wilcox.test(Exp ~ MM, data = df,alternative = "greater")
  return(data.frame(GENE = df$GENE[1], p_value = result$p.value))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

########

### first gathering the gene lists that would be used for common essentials, and context MM dependencies 

metadata <-  read.csv("../public_data/Model.csv")
Chronos <-  read.csv("../public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

Chronos <- Chronos %>% gather(GENE,Exp,-X)

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")

Chronos$MM <- 0 
Chronos$MM[Chronos$X %in% metadata$ModelID] <- 1

percent_dependent <- Chronos %>%
  group_by(GENE) %>%
  summarize(percentage = mean(Exp < -0.5) * 100)

result_df <- Chronos %>%
  group_by(GENE) %>%
  do(wilcox_test_function(.))

result_df$p_value_adj <- p.adjust(result_df$p_value)

result_df$p_value <- -log(result_df$p_value)

context_deps <- result_df %>% filter(p_value_adj < 0.01)
#### n = 113
core_ess <- percent_dependent %>% filter(percentage > 90)
#### n = 1048

#######

#### now reading in RNA-seq logFC that would be used for each bubble 

MM1S_KB528_1 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_1_MM1S_exprs_matrix.txt")
MM1S_KB528_3 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_3_MM1S_exprs_matrix.txt")
MM1S_KB528_6 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_6_MM1S_exprs_matrix.txt")
MM1S_KB528_24 <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_24_MM1S_exprs_matrix.txt")

MM1S_KB528_24$GENE <- row.names(MM1S_KB528_24)
MM1S_KB528_6$GENE <- row.names(MM1S_KB528_6)
MM1S_KB528_3$GENE <- row.names(MM1S_KB528_3)
MM1S_KB528_1$GENE <- row.names(MM1S_KB528_1)

colnames(MM1S_KB528_1)[3:4] <- c("Log2FC_KB528_1Hr","p_KB528_1Hr")
colnames(MM1S_KB528_3)[3:4] <- c("Log2FC_KB528_3Hr","p_KB528_3Hr")
colnames(MM1S_KB528_6)[3:4] <- c("Log2FC_KB528_6Hr","p_KB528_6Hr")
colnames(MM1S_KB528_24)[3:4] <- c("Log2FC_KB528_24Hr","p_KB528_24Hr")

dfs <- list(MM1S_KB528_1,MM1S_KB528_3,MM1S_KB528_6,MM1S_KB528_24)

dat <- dfs %>% purrr::reduce(full_join, by = "GENE") %>% as.data.frame()

rm(dfs,MM1S_KB528_1,MM1S_KB528_3,MM1S_KB528_6,MM1S_KB528_24)

# #### used for supplemental table.
# write_csv(dat,"./log2fc_MM1S_timeseries_KB528.csv")

##########
RPMI8226_KB528_1 <- read.delim("../RNA_data/KB528_time_series/KB528_RPMI8226_analysis_DMSO_RPMI8226_vs_KB528_1Hr_RPMI8226_exprs_matrix.txt")
RPMI8226_KB528_3 <- read.delim("../RNA_data/KB528_time_series/KB528_RPMI8226_analysis_DMSO_RPMI8226_vs_KB528_3Hr_RPMI8226_exprs_matrix.txt")
RPMI8226_KB528_6 <- read.delim("../RNA_data/KB528_time_series/KB528_RPMI8226_analysis_DMSO_RPMI8226_vs_KB528_6Hr_RPMI8226_exprs_matrix.txt")
RPMI8226_KB528_24 <- read.delim("../RNA_data/KB528_time_series/KB528_RPMI8226_analysis_DMSO_RPMI8226_vs_KB528_24Hr_RPMI8226_exprs_matrix.txt")

RPMI8226_KB528_24$GENE <- row.names(RPMI8226_KB528_24)
RPMI8226_KB528_6$GENE <- row.names(RPMI8226_KB528_6)
RPMI8226_KB528_3$GENE <- row.names(RPMI8226_KB528_3)
RPMI8226_KB528_1$GENE <- row.names(RPMI8226_KB528_1)

colnames(RPMI8226_KB528_1)[3:4] <- c("Log2FC_KB528_1Hr_RPMI8226","p_KB528_1Hr_RPMI8226")
colnames(RPMI8226_KB528_3)[3:4] <- c("Log2FC_KB528_3Hr_RPMI8226","p_KB528_3Hr_RPMI8226")
colnames(RPMI8226_KB528_6)[3:4] <- c("Log2FC_KB528_6Hr_RPMI8226","p_KB528_6Hr_RPMI8226")
colnames(RPMI8226_KB528_24)[3:4] <- c("Log2FC_KB528_24Hr_RPMI8226","p_KB528_24Hr_RPMI8226")

dfs <- list(RPMI8226_KB528_1,RPMI8226_KB528_3,RPMI8226_KB528_6,RPMI8226_KB528_24)

dat_RPMI8226 <- dfs %>% purrr::reduce(full_join, by = "GENE") %>% as.data.frame()

rm(dfs,RPMI8226_KB528_1,RPMI8226_KB528_3,RPMI8226_KB528_6,RPMI8226_KB528_24)
#write_csv(dat_RPMI8226,"./log2fc_RPMI8226_timeseries_KB528.csv")


##########

PBMC_KB528_1 <- read.delim("../RNA_data/KB528_time_series/PBMC_analysis_DMSO_0_vs_KB528_1_exprs_matrix.txt") %>% rownames_to_column("GENE")
PBMC_KB528_3 <- read.delim("../RNA_data/KB528_time_series/PBMC_analysis_DMSO_0_vs_KB528_3_exprs_matrix.txt") %>% rownames_to_column("GENE")
PBMC_KB528_6 <- read.delim("../RNA_data/KB528_time_series/PBMC_analysis_DMSO_0_vs_KB528_6_exprs_matrix.txt") %>% rownames_to_column("GENE")
PBMC_KB528_24 <- read.delim("../RNA_data/KB528_time_series/PBMC_analysis_DMSO_0_vs_KB528_24_exprs_matrix.txt") %>% rownames_to_column("GENE")

colnames(PBMC_KB528_1 )[4:5] <- c("Log2FC_PBMC_1Hr","p_PBMC_1Hr")
colnames(PBMC_KB528_3 )[4:5] <- c("Log2FC_PBMC_3Hr","p_PBMC_3Hr")
colnames(PBMC_KB528_6 )[4:5] <- c("Log2FC_PBMC_6Hr","p_PBMC_6Hr")
colnames(PBMC_KB528_24 )[4:5] <- c("Log2FC_PBMC_24Hr","p_PBMC_24Hr")

dfs <- list(PBMC_KB528_1,PBMC_KB528_3,PBMC_KB528_6,PBMC_KB528_24)

dat_pbmc <- dfs %>% purrr::reduce(full_join, by = "GENE") %>% as.data.frame()
dat_pbmc <- genematch(dat_pbmc)

rm(dfs,PBMC_KB528_1,PBMC_KB528_3,PBMC_KB528_6,PBMC_KB528_24)
#write_csv(dat_pbmc ,"./log2fc_PBMC_timeseries_KB528.csv")

dat <- merge(dat_pbmc,dat)
dat <- merge(dat_RPMI8226,dat)

########
### matching gene names

dat <- genematch(dat)
context_deps  <- genematch(context_deps )
core_ess <- genematch(core_ess)

dat_temp <- dat %>% dplyr::select(GENE,
                                  Log2FC_PBMC_1Hr,
                                  Log2FC_PBMC_3Hr,
                                  Log2FC_PBMC_6Hr,
                                  Log2FC_PBMC_24Hr,
                                  Log2FC_KB528_1Hr,
                                  Log2FC_KB528_3Hr,
                                  Log2FC_KB528_6Hr,
                                  Log2FC_KB528_24Hr,
                                  Log2FC_KB528_1Hr_RPMI8226,
                                  Log2FC_KB528_3Hr_RPMI8226,
                                  Log2FC_KB528_6Hr_RPMI8226,
                                  Log2FC_KB528_24Hr_RPMI8226
)

dat_temp <- dat_temp %>% gather(background, Log2FC,-GENE)

dat_temp_p <- dat %>% dplyr::select(GENE,
                                    p_PBMC_1Hr,
                                    p_PBMC_3Hr,
                                    p_PBMC_6Hr,
                                    p_PBMC_24Hr,
                                    p_KB528_1Hr,
                                    p_KB528_3Hr,
                                    p_KB528_6Hr,
                                    p_KB528_24Hr,
                                    p_KB528_1Hr_RPMI8226,
                                    p_KB528_3Hr_RPMI8226,
                                    p_KB528_6Hr_RPMI8226,
                                    p_KB528_24Hr_RPMI8226)

dat_temp_p <- dat_temp_p %>% gather(background_p, p,-GENE)

dat_temp$Time <- 0
dat_temp$Time[dat_temp$background %in% c("Log2FC_PBMC_1Hr","Log2FC_KB528_1Hr","Log2FC_KB528_1Hr_RPMI8226")] <- "1 Hour"
dat_temp$Time[dat_temp$background %in% c("Log2FC_PBMC_3Hr","Log2FC_KB528_3Hr","Log2FC_KB528_3Hr_RPMI8226")] <- "3 Hour"
dat_temp$Time[dat_temp$background %in% c("Log2FC_PBMC_6Hr","Log2FC_KB528_6Hr","Log2FC_KB528_6Hr_RPMI8226")] <- "6 Hour"
dat_temp$Time[dat_temp$background %in% c("Log2FC_PBMC_24Hr","Log2FC_KB528_24Hr","Log2FC_KB528_24Hr_RPMI8226")] <- "24 Hour"

dat_temp_p$Time <- 0
dat_temp_p$Time[dat_temp_p$background %in% c("p_PBMC_1Hr","p_KB528_1Hr","p_KB528_1Hr_RPMI8226")] <- "1 Hour"
dat_temp_p$Time[dat_temp_p$background %in% c("p_PBMC_3Hr","p_KB528_3Hr","p_KB528_3Hr_RPMI8226")] <- "3 Hour"
dat_temp_p$Time[dat_temp_p$background %in% c("p_PBMC_6Hr","p_KB528_6Hr","p_KB528_6Hr_RPMI8226")] <- "6 Hour"
dat_temp_p$Time[dat_temp_p$background %in% c("p_PBMC_24Hr","p_KB528_24Hr","p_KB528_24Hr_RPMI8226")] <- "24 Hour"

dat_temp$cell_type <- "Other"
dat_temp$cell_type[dat_temp$background %in% c("Log2FC_PBMC_1Hr","Log2FC_PBMC_3Hr","Log2FC_PBMC_6Hr","Log2FC_PBMC_24Hr")] <- "PBMC"
dat_temp$cell_type[dat_temp$background %in% c("Log2FC_KB528_1Hr","Log2FC_KB528_3Hr","Log2FC_KB528_6Hr","Log2FC_KB528_24Hr")] <- "MM1S"
dat_temp$cell_type[dat_temp$background %in% c("Log2FC_KB528_1Hr_RPMI8226","Log2FC_KB528_3Hr_RPMI8226","Log2FC_KB528_6Hr_RPMI8226","Log2FC_KB528_24Hr_RPMI8226")] <- "RPMI8226"

dat_temp_p$cell_type <- "Other"
dat_temp_p$cell_type[dat_temp_p$background %in% c("p_PBMC_1Hr","p_PBMC_3Hr","p_PBMC_6Hr","p_PBMC_24Hr")] <- "PBMC"
dat_temp_p$cell_type[dat_temp_p$background %in% c("p_KB528_1Hr","p_KB528_3Hr","p_KB528_6Hr","p_KB528_24Hr")] <- "MM1S"
dat_temp_p$cell_type[dat_temp_p$background %in% c("p_KB528_1Hr_RPMI8226","p_KB528_3Hr_RPMI8226","p_KB528_6Hr_RPMI8226","p_KB528_24Hr_RPMI8226")] <- "RPMI8226"

#### pathways are top 1000 irf4/p300 bound genes, common essential genes, and MM context dependencies

gene_list_boxplots <- read.csv("./data/gene_list_boxplots.csv")

dat_temp <- merge(dat_temp,dat_temp_p)
paths_h <- fgsea::gmtPathways("~/irf4_mm_trn_code/public_data/h.all.v7.4.symbols.gmt")
  
paths <- list(
              top = gene_list_boxplots$top_1000,
              mm_context = context_deps$GENE,
              common_ess = core_ess$GENE
              )

fgseaRes_final <- tibble( "pathway","pval","padj","log2err","ES","NES","size","leadingEdge","Time","type")
colnames(fgseaRes_final) <- c( "pathway","pval","padj","log2err","ES","NES","size","leadingEdge","Time","type")
fgseaRes_final<- fgseaRes_final %>% filter(type == "")

rm(paths_h)

### calculating gsea for specified pathways for each time point per cell line. 
### npermsimple = 10000 in order to calculate NES scores for one-sided GSEA npermsimple can be increased, but this increased compute time, and does not change interpretation of results, or direction of pathway
for(type in unique(dat_temp$cell_type)){
  for(hour in unique(dat_temp$Time)){
    t <- dat_temp %>% filter(Time == hour,cell_type == type)
    rankData <-  -t$Log2FC
    names(rankData) <-  t$GENE
    fgseaRes <- fgsea(paths,rankData,nPermSimple = 10000)
    fgseaRes$Time <- hour
    fgseaRes$type <- type
    fgseaRes_final <- rbind(fgseaRes,fgseaRes_final)
    write_csv(fgseaRes_final,"./data/fgsea_result_KB528.csv")
    
  }
}

########

rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(tidyverse)

setwd("~/irf4_mm_trn_code/figure_3") 

'%!in%' <- function(x,y)!('%in%'(x,y))

fgsea_result <- read.csv("./data/fgsea_result_KB528.csv")

### setting min pvals to 1e-50 gsea did not compute for 1 sample due to 1-sidedness of GSEA, and calculating inferred NES. 
top_scores <- fgsea_result %>% filter(pathway == "top")

### NES and ES are linearly related, since not all NES are calculated, they are inferred from ES.
top_scores %>% ggplot(aes(ES,NES)) + geom_point()

train = top_scores[!is.na(top_scores$NES),]
test = top_scores[is.na(top_scores$NES),]
test$NES <- NULL

model <- lm(NES~ES, data = train)
predicted_values <- predict(model, newdata = test)

top_scores$NES[is.na(top_scores$NES)] <- predicted_values

fgsea_result <- fgsea_result %>% filter(pathway != "top")

### NES and ES are linearly related, since not all NES are calculated, they are inferred from ES.
top_scores %>% ggplot(aes(ES,NES)) + geom_point()

fgsea_result <- rbind(fgsea_result,top_scores)

top_scores <- fgsea_result %>% filter(pathway == "common_ess")

top_scores %>% ggplot(aes(ES,NES)) + geom_point()

train = top_scores[!is.na(top_scores$NES),]
test = top_scores[is.na(top_scores$NES),]
test$NES <- NULL

model <- lm(NES~ES, data = train)
predicted_values <- predict(model, newdata = test)

top_scores$NES[is.na(top_scores$NES)] <- predicted_values

fgsea_result <- fgsea_result %>% filter(pathway != "common_ess")

fgsea_result <- rbind(fgsea_result,top_scores)

### setting min pvals to 1e-50 gsea did not compute for 1 sample due to 1-sidedness of GSEA. 
fgsea_result$pval[is.na(fgsea_result$pval)] <-  1.0e-50

fgsea_result$Time  <- factor(fgsea_result$Time  , levels = c("1 Hour","3 Hour","6 Hour","24 Hour"))

fgsea_result$pathway <- factor(fgsea_result$pathway ,
                                             levels = c("common_ess","top","mm_context"))

fgsea_result$log_pval <- -log(fgsea_result$pval)

fgsea_result$log_pval[fgsea_result$log_pval > 20] <- 20

write_csv(fgsea_result,"./data/fgsea_result_KB528_computed.csv")

### bubble plot in figure 3 

pdf("./fig3_gsea_dynamics.pdf",height = 5,width = 15)
fgsea_result %>%
  ggplot((aes(x = Time,y = pathway, size = log_pval, color = NES))) +
  geom_point() + facet_wrap(.~type) +
  scale_colour_gradient2() + theme_classic(base_size = 20) + scale_size(range = c(1, 15))
dev.off()

