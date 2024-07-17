### This code is used to generate the omic heatmap with arcs seen in figure 3. 
### The boxplots within figure 3E are generated in the R script figureS3_boxplots.R in the figure S3 folder.


rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)
library(fgsea)
library(circlize)
library(tidygraph)
library(ggraph)

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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggraph_2.2.1          tidygraph_1.3.1       circlize_0.4.16       fgsea_1.30.0          ggpubr_0.6.0          lubridate_1.9.3      
# [7] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.5           tidyr_1.3.1          
# [13] tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0       ComplexHeatmap_2.20.0 org.Hs.eg.db_3.19.1   AnnotationDbi_1.66.0 
# [19] IRanges_2.38.0        S4Vectors_0.42.0      Biobase_2.64.0        BiocGenerics_0.50.0  
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               gridExtra_2.3           rlang_1.1.3             magrittr_2.0.3          clue_0.3-65            
# [6] GetoptLong_1.0.5        DOSE_3.30.1             matrixStats_1.3.0       compiler_4.4.0          RSQLite_2.3.6          
# [11] png_0.1-8               vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         shape_1.4.6.1          
# [16] crayon_1.5.2            fastmap_1.1.1           backports_1.4.1         XVector_0.44.0          utf8_1.2.4             
# [21] HDO.db_0.99.1           tzdb_0.4.0              UCSC.utils_1.0.0        bit_4.0.5               zlibbioc_1.50.0        
# [26] cachem_1.0.8            GenomeInfoDb_1.40.0     jsonlite_1.8.8          blob_1.2.4              tweenr_2.0.3           
# [31] BiocParallel_1.38.0     broom_1.0.5             parallel_4.4.0          cluster_2.1.6           R6_2.5.1               
# [36] stringi_1.8.4           RColorBrewer_1.1-3      car_3.1-2               GOSemSim_2.30.0         Rcpp_1.0.12            
# [41] iterators_1.0.14        igraph_2.0.3            Matrix_1.7-0            splines_4.4.0           timechange_0.3.0       
# [46] tidyselect_1.2.1        viridis_0.6.5           qvalue_2.36.0           rstudioapi_0.16.0       abind_1.4-5            
# [51] doParallel_1.0.17       codetools_0.2-20        lattice_0.22-6          plyr_1.8.9              withr_3.0.0            
# [56] KEGGREST_1.44.0         polyclip_1.10-6         Biostrings_2.72.0       pillar_1.9.0            carData_3.0-5          
# [61] foreach_1.5.2           generics_0.1.3          hms_1.1.3               munsell_0.5.1           scales_1.3.0           
# [66] glue_1.7.0              tools_4.4.0             data.table_1.15.4       ggsignif_0.6.4          graphlayouts_1.1.1     
# [71] fs_1.6.4                fastmatch_1.1-4         cowplot_1.1.3           colorspace_2.1-0        GenomeInfoDbData_1.2.12
# [76] ggforce_0.4.2           cli_3.6.2               fansi_1.0.6             viridisLite_0.4.2       gtable_0.3.5           
# [81] rstatix_0.7.2           yulab.utils_0.1.4       digest_0.6.35           ggrepel_0.9.5           farver_2.1.1           
# [86] rjson_0.2.21            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7              GlobalOptions_0.1.2    
# [91] GO.db_3.19.1            MASS_7.3-60.2           bit64_4.0.5    

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

wilcox_test_function <- function(df) {
  result <- wilcox.test(Exp ~ MM, data = df)
  return(data.frame(GENE = df$GENE[1], p_value = result$p.value))
}

setwd("~/irf4_mm_trn_code/figure_3")

metadata <-  read.csv("../public_data/Model.csv")
Chronos <-  read.csv("../public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

Chronos <- Chronos %>% gather(GENE,Exp,-X)

metadata <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")

percent_dependent <- Chronos %>% filter(X %in% metadata$ModelID) %>%
  group_by(GENE) %>%
  summarize(percentage = mean(Exp < -0.5) * 100)

Chronos$MM <- 0 
Chronos$MM[Chronos$X %in% metadata$ModelID] <- 1

### Context specificity calculation for all genes. 

result_df <- Chronos %>%
  group_by(GENE) %>%
  do(wilcox_test_function(.))
result_df$p_value_adj <- p.adjust(result_df$p_value)

result_df$p_value <- -log(result_df$p_value)

##############

### This is the logFC, and p-values for RNA-Seq data at 3 hours
dat <- read.delim("../RNA_data/KB528_time_series/KB528_MM1S_analysis_DMSO_MM1S_vs_KB528_3_MM1S_exprs_matrix.txt")

### Filtering genes with TPM DMSO > 1
dat <- dat %>% filter(DMSO_MM1S > 1)

### enhancerPromoter output for ChIP data at H3K27ac DMSO peak loci 
enhancerPromoter_all <- read.csv("../ChIP_data/enhancerPromoter_DMSO_with_MYC.csv") %>% unique()

### combining signal to aggregate per gene
enhancerPromoter_all$comb_sig <- enhancerPromoter_all$DISTAL_signal + enhancerPromoter_all$TSS_signal

ep_dat <- enhancerPromoter_all %>% dplyr::select(GENE_ID,file,comb_sig) %>% spread(file,comb_sig)

### re ordering data so it its in order of replicates
ep_dat <- ep_dat[, c(1:2, (3+1):11, 3,12:ncol(ep_dat))]

ep_dat[is.na(ep_dat)] <- 0.0

### averaging the replicates by gene. 
mean_cols <- sapply(seq(2, ncol(ep_dat), by = 2), function(i) rowMeans(ep_dat[, i:(i+1)])) %>% as.data.frame()

colnames(mean_cols) <- colnames(ep_dat)[seq(3, ncol(ep_dat), by = 2)]
mean_cols$GENE <- ep_dat$GENE_ID

ep_dat <- mean_cols
ep_dat <- ep_dat %>% dplyr::select(GENE, everything())

### Identifying the rank of gene ChIP signal for p300 and IRF4

ep_dat$rank_ep300 <- rank(-ep_dat$p300H3K27ac_34_DMSO_compare) 
ep_dat$rank_IRF4 <- rank(-ep_dat$IRF4H3K27ac_32_DMSO_compare) 

ep_dat$tot_rank <- ep_dat$rank_ep300 + ep_dat$rank_IRF4

### setting DUSP22 to IRF4
ep_dat <- ep_dat %>% filter(GENE != "IRF4")
ep_dat$GENE[ep_dat$GENE == "DUSP22"] <- "IRF4"

ep_dat_p300_IRF4 <- ep_dat

rm(mean_cols)

########

### reading in DM3 reads for spike-in ChIP control. 
MM_p300_internal_QC_report <- read.csv("../ChIP_data/MM_p300_internal_QC_report.csv") 
MM_p300_internal_QC_report$Scale_factor <- MM_p300_internal_QC_report$Total_Reads_dm3/1000000
MM_p300_internal_QC_report$Scale_factor <- 1/MM_p300_internal_QC_report$Scale_factor

MM_p300_internal_QC_report$Scale_factor[21:34] <- 1

MM_p300_internal_QC_report$X <- NULL
MM_p300_internal_QC_report$X.1 <- NULL
MM_p300_internal_QC_report$X.2 <- NULL

##### This was used to create internal read stats found as a supplementary table .
#write_csv(MM_p300_internal_QC_report,"./MM_p300_internal_QC_report_supp_table.csv")

numbers <- gsub("_.*", "", MM_p300_internal_QC_report$Sample) %>% as.numeric()

MM_p300_internal_QC_report$sample_num <- numbers
########

enhancerPromoter_DMSO <- enhancerPromoter_all

enhancerPromoter_DMSO$sample_num <- enhancerPromoter_DMSO$file

enhancerPromoter_DMSO$sample_num <- gsub("_DMSO_compare","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k27ac_","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k18ac_","",enhancerPromoter_DMSO$sample_num)
enhancerPromoter_DMSO$sample_num <- gsub("H3k27me3_","",enhancerPromoter_DMSO$sample_num)

### identifying sample number in order to match against QC report for scaling.
enhancerPromoter_DMSO$sample_num <- as.numeric(enhancerPromoter_DMSO$sample_num)

MM_p300_internal_QC_report <- MM_p300_internal_QC_report %>% dplyr::select(sample_num,Scale_factor)

enhancerPromoter_DMSO <- merge(enhancerPromoter_DMSO,MM_p300_internal_QC_report)

### normalizing signal
enhancerPromoter_DMSO$comb_sig <- enhancerPromoter_DMSO$comb_sig*enhancerPromoter_DMSO$Scale_factor

ep_dat <- enhancerPromoter_DMSO %>% dplyr::select(GENE_ID,file,comb_sig) %>% spread(file,comb_sig)

ep_dat <- ep_dat[, c(1:2, (3+1):11, 3,12:ncol(ep_dat))]

ep_dat[is.na(ep_dat)] <- 0.0

mean_cols <- sapply(seq(2, ncol(ep_dat), by = 2), function(i) rowMeans(ep_dat[, i:(i+1)])) %>% as.data.frame()

colnames(mean_cols) <- colnames(ep_dat)[seq(3, ncol(ep_dat), by = 2)]
mean_cols$GENE <- ep_dat$GENE_ID

ep_dat <- mean_cols
ep_dat_DMSO <- ep_dat %>% dplyr::select(GENE, everything())


#### writing out supplementary table again for DMSO signal used. 
#write.csv(ep_dat_DMSO,"./ep_dat_DMSO.csv")

### log FC to track change of H3K27ac signal at time points

ep_dat_DMSO$H3K27ac_1Hr <- log(ep_dat_DMSO$H3k27ac_14_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_3Hr <- log(ep_dat_DMSO$H3k27ac_16_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_6Hr <- log(ep_dat_DMSO$H3k27ac_18_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)
ep_dat_DMSO$H3K27ac_24Hr <- log(ep_dat_DMSO$H3k27ac_20_DMSO_compare/ep_dat_DMSO$H3k27ac_12_DMSO_compare)

ep_dat_DMSO <- ep_dat_DMSO %>% dplyr::select(GENE,H3k27ac_12_DMSO_compare,H3K27ac_1Hr)

rm(enhancerPromoter_DMSO,mean_cols)


########

####normalizing data gene names in order to combine data together.
### RNA
dat$GENE <- rownames(dat)
dat <- genematch(dat)

###context specificity
result_df  <- genematch(result_df )

### percent dependency - not used in final figure
percent_dependent <- genematch(percent_dependent)

#### p300 and IRF4 signal 
ep_dat <- genematch(ep_dat_p300_IRF4)

### H3K27ac signal 
ep_dat_DMSO <- genematch(ep_dat_DMSO)

### p300H3K27ac_34_DMSO_compare p300 specific signal at H3K27ac DMSO loci
### IRF4H3K27ac_32_DMSO_compare p300 specific signal at H3K27ac DMSO loci
### tot_rank is the sum of the rank of p300 and IRF4 signal at H3K27ac DMSO loci
t1 <- ep_dat %>% dplyr::select(GENE,p300H3K27ac_34_DMSO_compare,IRF4H3K27ac_32_DMSO_compare,tot_rank)
### RNA
t2 <- dat %>% dplyr::select(GENE,LOG2_FOLD_CHANGE)
### DepMap
t3 <- percent_dependent %>% dplyr::select(GENE,percentage)
t4 <- result_df %>% dplyr::select(GENE,p_value_adj)
### H3K27ac signal logFC and total at DMSO 
t5 <- ep_dat_DMSO %>% dplyr::select(GENE,H3K27ac_1Hr,H3k27ac_12_DMSO_compare)

heatmap_dat <- merge(t1,t2)
heatmap_dat <- merge(heatmap_dat,t3)
heatmap_dat <- merge(heatmap_dat,t4)
heatmap_dat <- merge(heatmap_dat,t5)
heatmap_dat <- heatmap_dat %>% filter(!is.nan(H3K27ac_1Hr ))

### ordering heatmap to tot rank
order_indices <- order(heatmap_dat$tot_rank)
heatmap_dat <- heatmap_dat[order_indices, ]

### setting up color scales for heatmap
### p300 AUC ChIP signal
col_fun1 = colorRamp2(c(0,10000), c( "#FFFFFF00", "#00FFFF99"))

### IRF4 AUC ChIP signal
col_fun2 = colorRamp2(c(0,10000), c( "#FFFFFF00", "magenta"))

### H3K27ac AUC ChIP signal
col_fun3 = colorRamp2(c(0,10000), c( "#FFFFFF00", "black"))

### logFC RNA-seq
col_fun4 = colorRamp2(c(-1.5,-0.8, 0.8, 1.5), c("blue", "#FFFFFF00","#FFFFFF00", "red"),transparency = 0.5)
### percent dependent - not used
col_fun5 = colorRamp2(c(0,85,90), c("#FFFFFF00","#FFFFFF00","green"),transparency = 0.5)
### context specificity 
col_fun6 = colorRamp2(c(0,4), c("#FFFFFF00", "purple"),transparency = 0.5) ### cutoff/value shown is~p e-5, due to so many elements this was used in the graphic. 
### H3K27ac 1 Hr logFC
col_fun7 = colorRamp2(c(median(heatmap_dat$H3K27ac_1Hr[1:1000]),median(heatmap_dat$H3K27ac_1Hr[1001:nrow(heatmap_dat)]),0), c("slateblue4","#FFFFFF00","#FFFFFF00"))

### SE MM TRN TFS 
picked_genes <- c("ATF4","PRDM10","IRF4","PRDM1","MEF2C","MYC","POU2F1","TCF3","ZNF217")

heatmap_dat$label[heatmap_dat$GENE %in% picked_genes] <- heatmap_dat$GENE[heatmap_dat$GENE %in% picked_genes]
rownames(heatmap_dat)[heatmap_dat$GENE %in% picked_genes] <- heatmap_dat$GENE[heatmap_dat$GENE %in% picked_genes]

heatmap_dat_temp <- heatmap_dat %>% dplyr::select(p300H3K27ac_34_DMSO_compare,IRF4H3K27ac_32_DMSO_compare,H3k27ac_12_DMSO_compare,LOG2_FOLD_CHANGE,percentage,p_value_adj,H3K27ac_1Hr)
heatmap_dat_temp <- as.matrix(heatmap_dat_temp)
heatmap_dat_temp_rows <- heatmap_dat %>% dplyr::select(p300H3K27ac_34_DMSO_compare)

rownames(heatmap_dat_temp ) <- as.factor(rownames(heatmap_dat_temp ))

###splitting heatmap based on top 1000, middle and bottom 1000

labels = c(rep("1.Top", 1000), rep("2.Middle", nrow(heatmap_dat_temp) - 2000), rep("3.Bottom", 1000))

top_1000 <- heatmap_dat[1:1000,'GENE']
bot_1000 <- heatmap_dat[(nrow(heatmap_dat)-999):nrow(heatmap_dat),'GENE']

gene_lists_boxplots <- data_frame(top_1000,bot_1000 )

write_csv(gene_lists_boxplots,"./data/gene_list_boxplots.csv")

ht1 = Heatmap( heatmap_dat_temp[,1],  col = col_fun1, name = "p300 ChIP",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht2 = Heatmap( heatmap_dat_temp[,2], col = col_fun2, name = "IRF4 ChIP",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht3 = Heatmap( heatmap_dat_temp[,3], col = col_fun3, name = "H3K27ac DMSO",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht4 = Heatmap( heatmap_dat_temp[,4], col = col_fun4, name = "RNAseqFC_3Hr",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht5 = Heatmap( heatmap_dat_temp[,5], col = col_fun5, name = "Dependency_tot",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht6 = Heatmap(-log( heatmap_dat_temp[,6]), col = col_fun6, name = "Dependency_cont",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)
ht7 = Heatmap(heatmap_dat_temp[,7], col = col_fun7, name = "H3K27acFC_1Hr",show_row_names = F,cluster_rows = FALSE,cluster_columns =FALSE,row_split = labels)

p <- ht1 + ht2 + ht3 + ht4 + ht7 + ht6 + ht5

pdf("./fig3_heatmap.pdf",width = 5,height = 12)
p
dev.off()

### writing out supplementary table for heatmap data
write_csv(heatmap_dat,"./data/fig3_heatmap_dat_supp.csv")

### calculating statistics used in heatmap

cor.test(heatmap_dat$p300H3K27ac_34_DMSO_compare,
         heatmap_dat$IRF4H3K27ac_32_DMSO_compare)
# Pearson's product-moment correlation
# 
# data:  heatmap_dat$p300H3K27ac_34_DMSO_compare and heatmap_dat$IRF4H3K27ac_32_DMSO_compare
# t = 292.06, df = 8965, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9492503 0.9531890
# sample estimates:
#       cor 
# 0.9512584 

p<-cor.test(heatmap_dat$tot_rank,
         heatmap_dat$H3k27ac_12_DMSO_compare,
         method = "spearman")
p
###correlation being switched to positive

# data:  heatmap_dat$tot_rank and heatmap_dat$H3k27ac_12_DMSO_compare
# S = 2.3321e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.9406661 

cor.test(heatmap_dat$tot_rank,
         heatmap_dat$LOG2_FOLD_CHANGE,method = "spearman")

# data:  heatmap_dat$tot_rank and heatmap_dat$LOG2_FOLD_CHANGE
# S = 7.8695e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3451292 

cor.test(heatmap_dat$tot_rank,
         heatmap_dat$H3K27ac_1Hr,method = "spearman")

# Spearman's rank correlation rho
# 
# data:  heatmap_dat$tot_rank and heatmap_dat$H3K27ac_1Hr
# S = 6.5053e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4586525 

pathways <- gmtPathways("~/data/GO_pathways/c2.cp.kegg.v7.4.symbols.gmt")

context_deps <- result_df %>% dplyr::filter(GENE %in% heatmap_dat$GENE) %>% filter(p_value_adj < 0.01)
core_ess <- percent_dependent %>% dplyr::filter(GENE %in% heatmap_dat$GENE) %>% filter(percentage > 90)

pathways <- append(pathways,list(context_dep = c(context_deps$GENE)))
pathways <- append(pathways,list(core_ess = c(core_ess$GENE)))

ranks <- heatmap_dat %>% pull(-tot_rank)
names(ranks) <- heatmap_dat %>% pull(GENE)
p_val <- c()
p_val2 <- c()


### originally looking at other pathways, only used context dependency 
pathways <- pathways[c("context_dep")]

## with some zeros, we re ran it, but min p-value is still 1-e50.
for (i in 1:1000){

  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks)

  p_vals <- fgseaRes %>% filter(pathway == "context_dep") %>% pull(pval)
  p_val <- c(p_val,p_vals)

}

max(p_val)
min(p_val)
median(p_val)
mean(p_val)

# > max(p_val)
# [1] 1e-50
# 
# > min(p_val)
# [1] 1e-50
# 
# > median(p_val)
# [1] 1e-50
# 
# > mean(p_val)
# [1] 1e-50


#############

### This code is used to create the Arc plots

crc <- read.csv("../ChIP_data/MM1S_crc_all_chip.csv")
crc <- crc %>% filter(From %in% picked_genes,To %in% picked_genes,file %in% c("11_0EWX_0229Kronos_MM1S-No-Treatment-1_H3K27Ac_hs-dm_i20_EDGE_LIST.txt",
                                                                              "12_0EWY_0229Kronos_MM1S-No-Treatment-2_H3K27Ac_hs-dm_i21_EDGE_LIST.txt"))# %>% dplyr::select(From,file) %>% unique()

crc$comb <- paste(crc$From,crc$To)
tab <- table(crc$comb)
tab <- tab[tab > 1]
crc <- crc %>% filter(comb %in% names(tab)) %>% dplyr::select(From,To) %>% unique()

crc <- crc %>% filter(From != To)

# Concatenate the columns and identify unique rows
unique_rows <- !duplicated(apply(crc, 1, function(row) paste(sort(row), collapse = ',')))

# Filter the data frame
crc <- crc[unique_rows, ]

#### used for a supplemental table
#write_csv(crc,"./MM1.S_DMSO_crc.csv")

### data frame with groups, degree, labels and id, matches label order in heatmap
nodes <- data.frame(name = heatmap_dat$GENE,label = heatmap_dat$label)
nodes$label <- factor(nodes$label, levels = unique(nodes$label))

nodes$label_size <- 0
nodes$label_size[!is.na(nodes$label)] <- 20

nodes <- as_tibble(nodes)

edges <- as_tibble(crc)
edges$color = "gray"

edges$size = 0.2

net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE, node_key = "label")

#### creating the arcs on the heatmap

p2 <- ggraph(net.tidy, layout = "linear") + 
  geom_edge_arc(alpha = 0.5,force_flip = TRUE,fold = TRUE,aes(col = color,edge_width = size)) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_point(aes(size = label_size)) +
  scale_edge_color_identity(c("gray")) +
  #scale_edge_color_identity(c("gray","magenta")) +
  
  scale_colour_manual(values= edges$color) +
  geom_node_text(aes(label = label), angle = 90, hjust = 1, nudge_y = -0.2, size = 4) +   
  coord_cartesian(clip = "off") + 
  theme_graph() +
  theme(legend.position = "none") 

pdf("./fig3_arcs.pdf",width = 12,height = 4)
p2
dev.off()

