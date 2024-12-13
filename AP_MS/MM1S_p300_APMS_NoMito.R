library(MSstats)
library(tidyverse)
library(dplyr)
library(stringr)
library(UniProt.ws)
library(ggplot2)
library(ggrepel)

#Annotation
Human_TF <- read_csv("C:\\Users\\ypzhe\\Documents\\R_code\\Complex_annotation\\Human_TF_1639.csv")
Mito_compartment <-read_csv("C:\\Users\\ypzhe\\Documents\\R_code\\Complex_annotation\\Mito_annot_NoDup.csv")

############################################################
## Read Fragpipe MSstats export, filter out mito proteins
###########################################################

raw_FP <- read.csv("msstats_MM1S_B124re_S24A15_24MBR_default_p300.csv") |> 
  mutate(Protein_Acc = str_extract(ProteinName, "(?<=\\|).*(?=\\|)")) |> 
  anti_join(Mito_compartment, by = "Protein_Acc")
  
study_name <- "MM1S_3R24M_p300_noMito"

##############################
## Make MSstats required format
##############################
input_FP <- FragPipetoMSstatsFormat( 
                             input = raw_FP) 

##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################
## censoredInt='NA' for FP
processed_quant <- dataProcess(input_FP,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               featureSubset = "highQuality",
                               remove_uninformative_feature_outlier = TRUE,
                               min_feature_count = 2,
                               censoredInt="NA",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

Protein_hq <- processed_quant$ProteinLevelData
Peptide_hq <- processed_quant$FeatureLevelData
##################################################################
#protein ranking plot
################################################################
#filter out non-human protein  
#positive look behind assertion, ?<=
#positive look ahead assertion, ?=
pr_hq <- Protein_hq |>
  filter(str_detect(Protein,"HUMAN")) |>  
  mutate(Protein_Acc = str_extract(Protein, "(?<=\\|).*(?=\\|)"))

ID_mapping <- mapUniProt(
  from="UniProtKB_AC-ID",
  to='Gene_Name',
  query=unique(pr_hq$Protein_Acc)
) |> 
  rename (Protein_Acc = "From", Gene = "To") |> 
  aggregate(Gene ~ Protein_Acc, paste, collapse = ";")

ge_hq <- pr_hq |>
  left_join(ID_mapping, by = "Protein_Acc") |> 
  relocate(Gene, .before = 1)

ID_map <- ge_hq |> 
  dplyr::select(Gene, Protein) |> 
  unique() |> 
  mutate(Protein = as.character(Protein))

ge_hq$Gene <- str_replace_all(ge_hq$Gene,"(?=;).*","")
ge_hq$Gene <- coalesce(ge_hq$Gene, ge_hq$Protein_Acc)

write.csv(ge_hq,file=str_c(study_name, "_gene_long.csv"))

##############################################
## Model-based comparison + adjust p-value
#############################################

comparison <- matrix(c(-1,1,0,0,0,0,
                       0,0,-1,1,0,0,
                       0,0,0,0,-1,1
                       ), nrow=3, byrow = TRUE)

row.names(comparison) <- c("R1_p300-IgG",
                           "R2_p300-IgG",
                           "R3_p300-IgG"
)
groups <- levels(processed_quant$ProteinLevelData$GROUP)
colnames(comparison) <- groups

test_contrast <- groupComparison(contrast.matrix=comparison, data=processed_quant)

##################################################
## ID mapping from Uniprot accession to Gene name
##################################################
test_res <- test_contrast$ComparisonResult |> 
  mutate(Protein_Acc = str_extract(Protein, "(?<=\\|).*(?=\\|)"))

res <- test_res |>
  filter(is.na(issue)) |> 
  left_join(ID_mapping, by = "Protein_Acc") |> 
  relocate(c("Gene","Protein_Acc"), .before = 1) |> 
  dplyr::select(-Protein) |> 
  mutate (logpval = -log10(pvalue))|>
  mutate(logadjpval = -log10(adj.pvalue))
  
res$Gene <- coalesce(res$Gene, res$Protein_Acc) 

res2 <- res |> 
  group_by(Label) |> 
   mutate(logadjpval = ifelse(is.infinite(logadjpval), 
                             ceiling(max(logadjpval[is.finite(logadjpval)], na.rm = TRUE)), 
                             logadjpval)) |> 
  ungroup() |> 
  left_join(Human_TF, by = "Gene")

write_csv(res2, file=str_c(study_name,'_FP_byMSstats_EqMed_hq_DE.csv'))

by_Label <- res2 |> 
  group_nest(Label)

by_Label <- by_Label |> 
  mutate(path = str_glue("MM1S_5B-{Label}.csv"))
walk2(by_Label$data, by_Label$path, write_csv)

#######################################
#old code
R1_p300 <- res2 |> 
  filter(Label == "R1_p300-IgG")
R2_p300 <- res2 |> 
  filter(Label == "R2_p300-IgG")
R3_p300 <- res2 |> 
  filter(Label == "R3_p300-IgG")

write_csv(res2, file=str_c(study_name,'_FP_byMSstats_DE.csv'))
write_csv(R1_p300, file=str_c(study_name,'_FP_byMSstats_DE_R1_p300.csv'))
write_csv(R2_p300, file=str_c(study_name,'_FP_byMSstats_DE_R2_p300.csv'))
write_csv(R3_p300, file=str_c(study_name,'_FP_byMSstats_DE_R3_p300.csv'))

res2 <- read.csv("MM1S_3R24M_p300_noMito_FP_byMSstats_DE.csv")
################################################
#Calculating AP-MS stoichiometry 
################################################
APMS_ratio_ref <- ge_hq |> 
  group_by(RUN) |> 
  mutate(ratio_p300 = 2^(LogIntensities - LogIntensities[Gene=="EP300"])
  ) |> 
  ungroup()

APMS_ratio <- APMS_ratio_ref |> 
  group_by(GROUP,Gene, Protein_Acc) |> 
  summarise(p300_ratio_mean = mean(ratio_p300,  na.rm = TRUE),
            p300_ratio_cv = sd(ratio_p300,na.rm = TRUE)/p300_ratio_mean
  ) |> 
  ungroup()

write.csv(APMS_ratio, str_c(study_name,"_APMS_protein_ratio_by_GROUP.csv"))

########################################################
#Load lysate stoichiometry 
#1. protein fraction to total with each run
#2. mean within GROUP
#3. Stoichiometry to target protein
#######################################################
Sto_lysate_p300 <- read_csv("MM1S_4batches_lysate_profiling_lysate_sto.csv") |> 
  dplyr::select(Gene, Protein_Acc,GROUP,Sto_p300) |> 
  mutate(GROUP = case_match(GROUP,
                            "MM1S_B1_SAN24h" ~ "R1",
                            "MM1S_B2_SAN24h" ~ "R2",
                            "MM1S_B4_SAN24h" ~ "R3",
                            "MM1S_B5_SAN24h" ~ "B5"
  )) |>
  rename(Run = GROUP, Lysate_sto = Sto_p300) |> 
  filter(Run != "B5")

###############################################################################
#AP-MS and lysate stoichiometry
###############################################################################
Sto_APMS_p300 <- APMS_ratio |> 
  separate_wider_delim(GROUP,delim = "_", names = c('Run', 'Ab'))|> 
  filter(Ab == 'p300') |> 
  dplyr::select('Gene','Protein_Acc','Run','p300_ratio_mean') |> 
  rename(APMS_sto = p300_ratio_mean)

Sig_APMS_p300 <- res2 |> 
  separate_wider_delim(Label,delim = "_", names = c('Run','APMS')) |> 
  dplyr::select('Gene','Protein_Acc','Run','log2FC','logadjpval','ImputationPercentage', 'TF_1639')

p300_sto <- Sto_APMS_p300 |> 
  left_join(Sto_lysate_p300, by = c('Gene','Protein_Acc','Run')) |> 
  left_join(Sig_APMS_p300, by = c('Gene','Protein_Acc','Run')) 

write.csv(p300_sto, str_c(study_name, "_p300_stoichiometry.csv"))

p300_sto_sig <- p300_sto |> 
  filter(logadjpval > 2 & log2FC > 1)

p300_sto_sig <- p300_sto_sig |> 
  group_by(Run) |> 
  mutate(Lysate_sto = replace_na(Lysate_sto, 0.01),
         TF = ifelse(!is.na(TF_1639), "TF","Cofactor"))

write.csv(p300_sto_sig, "MM1S_B1234_SAN24h_Ab15m_16MBR_B3B4rerun_p300_sto_sig.csv")


group_counts <- p300_sto_sig |> 
  group_by(Run,TF) |> 
  count()

new_labels_Run <- c("R1" = str_c("R1 [", 
                                 group_counts[group_counts$Run=="R1",]$n[1],
                                 ", ",
                                 group_counts[group_counts$Run=="R1",]$n[2],
                                 "]"),
                    "R2" = str_c("R2 [", 
                                 group_counts[group_counts$Run=="R2",]$n[1],
                                 ", ",
                                 group_counts[group_counts$Run=="R2",]$n[2],
                                 "]"),
                    "R3" = str_c("R3 [", 
                                 group_counts[group_counts$Run=="R3",]$n[1],
                                 ", ",
                                 group_counts[group_counts$Run=="R3",]$n[2],
                                 "]"))

p300_sto_sig |> 
  ggplot (aes (log10(Lysate_sto),log10(APMS_sto)))+ 
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(color = 'grey25', alpha = 0.3)+
  geom_point(data = p300_sto_sig|> 
               filter(Gene %in% c('EP300','CREBBP', 'KHDRBS1', 'SMARCC2','SMARCE1','RUNX2','ARID3B','ZBTB38','ARID2','TFAP4','MEF2D','TCF12','CXXC1','ZNF207')), 
             color = 'blue')+
  geom_point(data = p300_sto_sig|> 
               filter(Gene == "EP300"),
             color = "red", shape = 15, size = 3)+
  geom_text_repel(data = p300_sto_sig|> 
                    filter(Gene %in% c('EP300','CREBBP', 'KHDRBS1', 'SMARCC2','SMARCE1','RUNX2','ARID3B','ZBTB38','ARID2','TFAP4','MEF2D','TCF12','CXXC1','ZNF207')), 
                  aes (label = Gene), size = 4, max.overlaps = 400)+
  facet_grid(TF~Run,labeller = labeller(Run = new_labels_Run))+
  labs (x = "Abundance stoichiometry, % of p300", 
        y = "Interaction stoichiometry, % of p300",
        title = "p300 interactors in MM1S, Adj.Pval<0.01, LFC>1")+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw()+ 
  theme(legend.position = 'none',
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  scale_x_continuous(breaks=c(-2, -1, 0, 1), 
                     labels=c("-2" = "1%", "-1" = "10%", "0" = "1x", "1" = "10x"))+
  scale_y_continuous(breaks=c(-2, -1, 0,1), 
                     labels=c("-2" = "1%", "-1" = "10%", "0" = "1x", "1" = "10x"))    
ggsave("p300_stoichiometry_MM1S_B124r_Batch_vs_TF_labeled_LFC1_NewLabel.png", width = 8, height = 6)




