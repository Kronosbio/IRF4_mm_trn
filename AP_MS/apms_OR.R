rm(list = ls())

if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

gc()

setwd("~/ash_2024/")

library(tidyverse)
library(scales)
library(patchwork)
library(eulerr)
library(gridExtra)
library(igraph)
library(org.Hs.eg.db)
'%!in%' <- function(x,y)!('%in%'(x,y))

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

p300_sig <- read.csv("~/ash_2024/apms/p300_3R24M_noMito_average_score.csv")
IRF4_sig <- read.csv("~/ash_2024/apms/IRF4_3R24M_noMito_average_score.csv")

p300_all <- read.csv("~/ash_2024/apms/MM1S_3R24M_p300_noMito_FP_byMSstats_DE.csv")
IRF4_all <- read.csv("~/ash_2024/apms/MM1S_3R24M_IRF4_noMito_FP_byMSstats_DE.csv")
IKZF1_all <- read.csv("~/ash_2024/apms/MM1S_B5ran_IKZF_noMito_FP_byMSstats_DE_IKZF1.csv")
IKZF3_all <- read.csv("~/ash_2024/apms/MM1S_B5ran_IKZF_noMito_FP_byMSstats_DE_IKZF3.csv")

### Values set by our proteomics lead Yupeng
IKZF1_sig <- IKZF1_all %>% filter(logadjpval > 1.3, log2FC > 0.5)
IKZF3_sig <- IKZF3_all %>% filter(logadjpval > 1.3, log2FC > 0.5)

dat <- data_frame(gene = c("p300","IRF4","IKZF1","IKZF3"),sig = c(nrow(p300_sig),nrow(IRF4_sig),nrow(IKZF1_sig),nrow(IKZF3_sig)))

p300_genes <- p300_sig %>% pull(Prey) %>% unique()
p300_bg_gene <- p300_all %>% pull(Gene) %>% unique()

IRF4_genes <- IRF4_sig %>% pull(Prey) %>% unique()
IRF4_bg_gene <- IRF4_all %>% pull(Gene) %>% unique()

IKZF1_bg_genes <- IKZF1_all %>% pull(Gene) %>% unique()
IKZF3_bg_genes <- IKZF3_all %>% pull(Gene) %>% unique()

background_set <- c(p300_bg_gene,IRF4_bg_gene,IKZF1_bg_genes,IKZF3_bg_genes) %>% unique()

IKZF1_genes <- IKZF1_sig %>% pull(Gene) %>% unique() 
IKZF3_genes <- IKZF3_sig %>% pull(Gene) %>% unique() 

IKZF1_df <- IKZF1_sig %>% dplyr::select(Gene)
colnames(IKZF1_df)[1] <- "GENE"
IKZF1_df$tf <- "IKZF1"
IKZF3_df <- IKZF3_sig %>% dplyr::select(Gene)
colnames(IKZF3_df)[1] <- "GENE"
IKZF3_df$tf <- "IKZF3"
EP300_df <- p300_sig %>% dplyr::select(Prey)
colnames(EP300_df)[1] <- "GENE"
EP300_df$tf <- "EP300"
IRF4_df <- IRF4_sig %>% dplyr::select(Prey)
colnames(IRF4_df)[1] <- "GENE"
IRF4_df$tf <- "IRF4"
combined_interactors <- rbind(IKZF1_df,IKZF3_df,EP300_df,IRF4_df) %>% unique()
combined_interactors <- genematch(combined_interactors)


wilcox_test_function <- function(df) {
  result <- wilcox.test(Exp ~ MM, data = df)
  return(data.frame(GENE = df$GENE[1], p_value = result$p.value))
}

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

String_ppi <- PPI_links %>% filter(GENE1 %in% c("IRF4","EP300","IKZF1","IKZF3"))
String_ppi$GENE <- String_ppi$GENE2
String_ppi <- genematch(String_ppi)

metadata <-  read.csv("~/irf4_mm_trn_code/public_data/Model.csv")
Chronos <-  read.csv("~/irf4_mm_trn_code/public_data/CRISPR_(DepMap_Public_22Q4+Score,_Chronos).csv")

metadata_MM <- metadata %>% filter(OncotreeSubtype == "Plasma Cell Myeloma")

Chronos <- Chronos %>% gather(GENE,Exp,-X)
Chronos$MM <- 0 
Chronos$MM[Chronos$X %in% metadata_MM$ModelID] <- 1

result_df_MM <- Chronos %>%
  group_by(GENE) %>%
  do(wilcox_test_function(.))

result_df_MM  <- genematch(result_df_MM)

### chose the length of p300 bg list as the number of top "MM" dependencies, all background sets seemed to be around ~4k genes.
### since not all genes are available for protein measurement, the bg set was used to identify total "MM context dependencies"
### context dependency cut off is ~0.045 p_value for context dependency score 

result_df_MM <- result_df_MM %>% arrange(p_value) %>% head(length(p300_bg_gene))

ppi_bg_mm <- PPI_links %>% filter(GENE1 %in% result_df_MM$GENE) %>% dplyr::select(GENE1) %>% unique() %>% nrow()
p300_bg_mm <- p300_bg_gene[p300_bg_gene %in% result_df_MM$GENE]
IRF4_bg_mm <- IRF4_bg_gene[IRF4_bg_gene %in% result_df_MM$GENE]

String_ppi$MM_dep <- "sig"
String_ppi$MM_dep[String_ppi$GENE %in% result_df_MM$GENE] <- "sig_mm"
combined_interactors$MM_dep <- "sig"
combined_interactors$MM_dep[combined_interactors$GENE %in% result_df_MM$GENE] <- "sig_mm"

odds <- c()
ci_lowers <- c()
ci_uppers <- c()
pvals <- c()

for (gene in c("IRF4","IKZF1","IKZF3","EP300")){
  interactors_t_1 <- combined_interactors %>% filter(tf %in% c(gene),MM_dep == "sig_mm") %>% nrow()
  interactors_t_0 <- combined_interactors %>% filter(tf %in% c(gene),MM_dep != "sig_mm") %>% nrow()
  string_ppi_t_1 <- String_ppi %>% filter(GENE1 %in% c(gene),MM_dep == "sig_mm") %>% nrow()
  string_ppi_t_0 <- String_ppi %>% filter(GENE1 %in% c(gene),MM_dep != "sig_mm") %>% nrow()
  
  group1 <- c(interactors_t_1,interactors_t_0)  # Group 1: [Count of 1's, Count of 0's]
  group2 <- c(string_ppi_t_1,string_ppi_t_0)  # Group 2: [Count of 1's, Count of 0's]
  odds_ratio <- (interactors_t_1 /interactors_t_0) / (string_ppi_t_1 / string_ppi_t_0)
  
  # SE
  se_ln_or <- sqrt(1/interactors_t_1 + 1/interactors_t_0 + 1/string_ppi_t_1 + 1/string_ppi_t_0)
  
  # Log-OR Confidence Intervals
  z <- 1.96 # for 95% confidence
  ln_or <- log(odds_ratio)
  ci_lower_ln <- ln_or - z * se_ln_or
  ci_upper_ln <- ln_or + z * se_ln_or
  
  # Back-transform to original scale
  ci_lower <- exp(ci_lower_ln)
  ci_upper <- exp(ci_upper_ln)
  
  contingency_table <- matrix(c(group1, group2), nrow = 2, byrow = TRUE)

  chisq <- chisq.test(contingency_table, correct = FALSE)
  chisq$p.value
  
  cat("Odds Ratio:", odds_ratio, "\n")
  cat("95% CI: (", ci_lower, ", ", ci_upper, ")\n")
  
  odds <- c(odds, odds_ratio)
  ci_lowers <- c(ci_lowers, ci_lower)
  ci_uppers <- c(ci_uppers, ci_upper)
  pvals <- c(pvals, chisq$p.value)

}

odds_df <- data.frame(gene = c("IRF4","IKZF1","IKZF3","EP300"), 
                      odds = odds, 
                      ci_lower = ci_lowers, 
                      ci_upper = ci_uppers,
                      p_value = pvals)

p <- odds_df %>% ggplot(aes(y = gene, x = odds)) + geom_point(size = 10) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, size = 2)) + theme_classic(base_size = 20) +
 theme(legend.position = "none") + xlab("Odds Ratio") + ylab("Gene") +
  ggtitle("APMS MM Context Dependency Hit rate") + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept =1, linetype = "dashed", color = "black", size = 0.8) 

pdf("odds_ratio.pdf",height = 5, width = 10)
p
dev.off()

String_IRF4_MM <- String_ppi %>% filter(GENE1 == "IRF4",MM_dep == "sig_mm") %>% droplevels()
String_EP300_MM <- String_ppi %>% filter(GENE1 == "EP300",MM_dep == "sig_mm") %>% droplevels()
interactors_p300_MM <- combined_interactors %>% filter(tf == "EP300",MM_dep == "sig_mm") %>% droplevels()
interactors_IRF4_MM <- combined_interactors %>% filter(tf == "IRF4",MM_dep == "sig_mm") %>% droplevels()

String_IRF4 <- String_ppi %>% filter(GENE1 == "IRF4",MM_dep != "sig_mm") %>% droplevels()
String_EP300 <- String_ppi %>% filter(GENE1 == "EP300",MM_dep != "sig_mm") %>% droplevels()
interactors_p300 <- combined_interactors %>% filter(tf == "EP300",MM_dep != "sig_mm") %>% droplevels()
interactors_IRF4 <- combined_interactors %>% filter(tf == "IRF4",MM_dep != "sig_mm") %>% droplevels()

GENE1 = c("p300","p300","p300","p300","IRF4","IRF4")
GENE2 = c("MM_unique_p300","non_MM_p300","shared_MM","shared","MM_unique_IRF4","non_MM_IRF4")
shared_mm <- String_EP300_MM$GENE[String_EP300_MM$GENE %in% String_IRF4_MM$GENE] %>% length()
shared_o <- String_EP300$GENE[String_EP300$GENE %in% String_IRF4$GENE] %>% length()
p300_mm <- String_EP300_MM$GENE[String_EP300_MM$GENE %!in% String_IRF4_MM$GENE] %>% length()
p300_o <- String_EP300$GENE[String_EP300$GENE %!in% String_IRF4$GENE] %>% length()
IRF4_mm <- String_IRF4_MM$GENE[String_IRF4_MM$GENE %!in% String_EP300_MM$GENE] %>% length()
IRF4_o <- String_IRF4$GENE[String_IRF4$GENE %!in% String_EP300$GENE] %>% length()
count = c(p300_mm,p300_o,shared_mm,shared_o,IRF4_mm,IRF4_o)
string_shared_int <- data_frame(GENE1,GENE2,count)

GENE1 = c("p300","p300","p300","p300","IRF4","IRF4")
GENE2 = c("MM_unique_p300","non_MM_p300","shared_MM","shared","MM_unique_IRF4","non_MM_IRF4")
shared_mm <- interactors_p300_MM$GENE[interactors_p300_MM$GENE %in% interactors_IRF4_MM$GENE] %>% length()
shared_o <- interactors_p300$GENE[interactors_p300$GENE %in% interactors_IRF4$GENE] %>% length()
p300_mm <- interactors_p300_MM$GENE[interactors_p300_MM$GENE %!in% interactors_IRF4_MM$GENE] %>% length()
p300_o <- interactors_p300$GENE[interactors_p300$GENE %!in% interactors_IRF4$GENE] %>% length()
IRF4_mm <- interactors_IRF4_MM$GENE[interactors_IRF4_MM$GENE %!in% interactors_p300_MM$GENE] %>% length()
IRF4_o <- interactors_IRF4$GENE[interactors_IRF4$GENE %!in% interactors_p300$GENE] %>% length()
count = c(p300_mm,p300_o,shared_mm,shared_o,IRF4_mm,IRF4_o)
apms_shared_int <- data_frame(GENE1,GENE2,count)

write_delim(string_shared_int,"string_shared_int.csv",delim = ",")
write_delim(apms_shared_int,"apms_shared_int.csv",delim = ",")

venn_data_interactors <- list(
  IRF4 = interactors_IRF4_MM$GENE,
  p300 = interactors_p300_MM$GENE
)

venn_data_STRING <- list(
  IRF4 = String_IRF4_MM$GENE,
  p300 = String_EP300_MM$GENE
)

interactors_tot <- interactors_IRF4_MM %>% filter(GENE %in% interactors_p300_MM$GENE)

write_delim(interactors_tot,"./interactors_apms_MM_deps.csv",delim = ",")

euler_STRING <- euler(venn_data_STRING)
euler_interactors <- euler(venn_data_interactors)

common_params <- list(
  fills = c("lavender", "lavender","#54278F"),
  alpha = 0.5,
  quantities = list(fontsize = 18), # Same font size for quantities
  labels = list(fontsize = 18)     # Same font size for labels
)

# Plot diagrams with the same settings
pdf("euler_STRING.pdf",height = 5, width = 7.5)
plot(euler_STRING, main = "Small Diagram",   fills = c("lavender", "lavender","#54278F"),  alpha = 0.5, quantities = list(fontsize = 18), labels = list(fontsize = 18))
dev.off()

a <- 22
b <- 324
c <- 11
d <- ppi_bg_mm - a - b - c

t <- chisq.test(matrix(c(a,b,c,d), nrow = 2, byrow = TRUE), correct = TRUE)
t$p.value


pdf("euler_APMS.pdf",height = 5, width = 7.5)
plot(euler_interactors, main = "Small Diagram",   fills = c("lavender", "lavender","#54278F"),  alpha = 0.5, quantities = list(fontsize = 18), labels = list(fontsize = 18))
dev.off()

a <- 58
b <- 82
c <- 146
d <- length(p300_bg_mm) - a - b - c

t <- chisq.test(matrix(c(a,b,c,d), nrow = 2, byrow = TRUE), correct = TRUE)
t$p.value


