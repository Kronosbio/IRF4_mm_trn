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

raw_FP <- read.csv("msstats_MM1S_B5ran_IKZF1_IKZF3_12MBR.csv") |>
  mutate(Protein_Acc = str_extract(ProteinName, "(?<=\\|).*(?=\\|)")) |> 
  anti_join(Mito_compartment, by = "Protein_Acc")

study_name <- "MM1S_B5ran_IKZF_noMito"

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

comparison <- matrix(c(-1,1,0,
                       -1,0,1
), nrow=2, byrow = TRUE)

row.names(comparison) <- c("B5r_IKZF1-IgG",
                           "B5r_IKZF3-IgG")
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

#######################################
#old code
B5_IKZF1 <- res2 |> 
  filter(Label == "B5r_IKZF1-IgG")
B5_IKZF3 <- res2 |> 
  filter(Label == "B5r_IKZF3-IgG")

write_csv(res2, file=str_c(study_name,'_FP_byMSstats_DE.csv'))
write_csv(B5_IKZF1, file=str_c(study_name,'_FP_byMSstats_DE_B1_IKZF1.csv'))
write_csv(B5_IKZF3, file=str_c(study_name,'_FP_byMSstats_DE_B2_IKZF3.csv'))

################################################
#Calculating AP-MS stoichiometry 
################################################
APMS_ratio_ref <- ge_hq |> 
  group_by(RUN) |> 
  mutate(ratio_IKZF1 = 2^(LogIntensities - LogIntensities[Gene=="IKZF1"]),
         ratio_IKZF3 = 2^(LogIntensities - LogIntensities[Gene=="IKZF3"])) |> 
  ungroup()

APMS_ratio <- APMS_ratio_ref |> 
  group_by(GROUP,Gene, Protein_Acc) |> 
  summarise(IKZF1_ratio_mean = mean(ratio_IKZF1,  na.rm = TRUE),
            IKZF3_ratio_mean = mean(ratio_IKZF3,  na.rm = TRUE),
            IKZF1_ratio_cv = sd(ratio_IKZF1,na.rm = TRUE)/IKZF1_ratio_mean,
            IKZF3_ratio_cv = sd(ratio_IKZF3,na.rm = TRUE)/IKZF3_ratio_mean) |> 
  ungroup()

write.csv(APMS_ratio, str_c(study_name,"_APMS_protein_ratio_by_GROUP.csv"))

########################################################
#Load lysate stoichiometry 
#1. protein fraction to total with each run
#2. mean within GROUP
#3. Stoichiometry to target protein
#######################################################
Sto_lysate_IKZF <- read_csv("MM1S_4batches_lysate_profiling_lysate_sto.csv") |> 
  dplyr::select(Gene, Protein_Acc,GROUP,Sto_IKZF1, Sto_IKZF3) |> 
  mutate(GROUP = case_match(GROUP,
                            "MM1S_B1_SAN24h" ~ "R1",
                            "MM1S_B2_SAN24h" ~ "R2",
                            "MM1S_B4_SAN24h" ~ "R3",
                            "MM1S_B5_SAN24h" ~ "B5"
  )) |>
  filter(GROUP == "B5") |> 
  dplyr::select(-GROUP) |> 
  pivot_longer(!c('Gene','Protein_Acc'), names_to = "Bait", values_to = "Lysate_sto") |> 
  separate_wider_delim(Bait,delim = "_", names = c(NA, 'Bait'))

###############################################################################
#AP-MS and lysate stoichiometry
###############################################################################
Sto_APMS_IKZF1 <- APMS_ratio |> 
  separate_wider_delim(GROUP,delim = "_", names = c('Run', 'Ab'))|> 
  filter(Ab == 'IKZF1') |> 
  dplyr::select('Gene','Protein_Acc','Ab','IKZF1_ratio_mean') |> 
  rename(APMS_sto = IKZF1_ratio_mean)

Sto_APMS_IKZF3 <- APMS_ratio |> 
  separate_wider_delim(GROUP,delim = "_", names = c('Run', 'Ab'))|> 
  filter(Ab == 'IKZF3') |> 
  dplyr::select('Gene','Protein_Acc','Ab','IKZF3_ratio_mean') |> 
  rename(APMS_sto = IKZF3_ratio_mean)

Sto_APMS_IKZF <-rbind(Sto_APMS_IKZF1, Sto_APMS_IKZF3) |> 
  rename(Bait = Ab)

Sig_APMS_IKZF <- res2 |> 
  separate_wider_delim(Label,delim = "_", names = c(NA,'APMS')) |>
  separate_wider_delim(APMS,delim = "-", names = c('Bait',NA))|>
  dplyr::select('Gene','Protein_Acc','Bait','log2FC','logadjpval','ImputationPercentage', 'TF_1639')

IKZF_sto <- Sto_APMS_IKZF |> 
  left_join(Sto_lysate_IKZF, by = c('Gene','Protein_Acc','Bait')) |> 
  left_join(Sig_APMS_IKZF, by = c('Gene','Protein_Acc','Bait')) 

write.csv(IKZF_sto, str_c(study_name, "_IKZF_stoichiometry.csv"))

IKZF_sto_sig <- IKZF_sto |> 
  filter(logadjpval > 1.3 & log2FC > 0.5)

IKZF_sto_sig <- IKZF_sto_sig |> 
  group_by(Bait) |> 
  mutate(Lysate_sto = replace_na(Lysate_sto, 0.01),
         TF = ifelse(!is.na(TF_1639), "TF","Cofactor"))


write.csv(IKZF_sto_sig, "MM1S_B5ran_SAN24h_Ab15m_12MBR_IKZF_sto_sig.csv")


group_counts <- IKZF_sto_sig |> 
  group_by(Bait,TF) |> 
  count()

new_labels_Run <- c("IKZF1" = str_c("IKZF1 [", 
                                    group_counts[group_counts$Bait=="IKZF1",]$n[1],
                                    ", ",
                                    group_counts[group_counts$Bait=="IKZF1",]$n[2],
                                    "]"),
                    "IKZF3" = str_c("IKZF3 [", 
                                    group_counts[group_counts$Bait=="IKZF3",]$n[1],
                                    ", ",
                                    group_counts[group_counts$Bait=="IKZF3",]$n[2],
                                    "]"))

IKZF_sto_sig |> 
  ggplot (aes (log10(Lysate_sto),log10(APMS_sto)))+ 
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(color = 'grey25', alpha = 0.3)+
  geom_point(data = IKZF_sto_sig|> 
               filter(Gene %in% c('IKZF1','IKZF2','IKZF3', 'IKZF4', 'IKZF5','RUNX1','RUNX2','RUNX3','ZBTB38','ZNF217','MEF2C','NFATC2','WDR5','KMT2C','KAT5','EP400')), 
             color = 'blue')+
  geom_point(data = IKZF_sto_sig|> 
               filter((Bait == "IKZF1"&Gene == "IKZF1") | (Bait == "IKZF3"&Gene == "IKZF3")),
             color = "red", shape = 15, size = 3)+
  geom_text_repel(data = IKZF_sto_sig|> 
                    filter(Gene %in% c('IKZF1','IKZF2','IKZF3', 'IKZF4', 'IKZF5','RUNX1','RUNX2','RUNX3','ZBTB38','ZNF217','MEF2C','NFATC2','WDR5','KMT2C','KAT5','EP400')), 
                  aes (label = Gene), size = 4, max.overlaps = 400)+
  facet_grid(TF~Bait,labeller = labeller(Bait = new_labels_Run))+
  labs (x = "Abundance stoichiometry, % of Bait", 
        y = "Interaction stoichiometry, % of Bait",
        title = "IKZF1 and IKZF3 interactors in MM1S, Adj.Pval<0.05, LFC>0.5")+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw()+ 
  theme(legend.position = 'none',
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  scale_x_continuous(breaks=c(-2, -1, 0,1), 
                     labels=c("-2" = "1%", "-1" = "10%", "0" = "1x", "1" = "10x"))+
  scale_y_continuous(breaks=c(-3,-2, -1, 0), 
                     labels=c("-3" = "0.1%", "-2" = "1%", "-1" = "10%", "0" = "1x"))    
ggsave("IKZF_stoichiometry_MM1S_B5ran_Bait_vs_TF_labeled_LFCp5.png", width = 6, height = 6)



