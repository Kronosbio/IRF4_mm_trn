library(MSstats)
library(tidyverse)
library(dplyr)
library(stringr)
library(UniProt.ws)

#Annotation
Human_TF <- read_csv("C:\\Users\\ypzhe\\Documents\\R_code\\Complex_annotation\\Human_TF_1639.csv")
Mito_compartment <-read_csv("C:\\Users\\ypzhe\\Documents\\R_code\\Complex_annotation\\Mito_annot_NoDup.csv")

##############################
## Read Fragpipe MSstats export
##############################

raw_FP <- read_csv("msstats_MM1S_B1245_SAN24h_lysate_profiling_6MBR.csv")
study_name <- "MM1S_4batches_lysate_profiling"

##############################
## Make MSstats required format
##############################
input_FP <- FragPipetoMSstatsFormat( 
                             input = raw_FP) 

#################################################################################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
#################################################################################

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

###########################################################################################
#Gene long format
#filter out non-human protein

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

###################################
#filter out CTNNB1 and PLCH1 (single feature high abundance, likely artifact)

ge_hq <- ge_hq |> 
  filter(!Gene %in% c('CTNNB1','PLCH1'))

write.csv(ge_hq,file=str_c(study_name, "_gene_long.csv"))

########################################################
#Calculate lysate stoichiometry from protein level data
#######################################################
#Two ways to calculate lysate stoichiometry
# 1.calculate ratio of each protein against target protein (IRF4, p300, IKZF1,3, etc)for each run, then take mean from 3 bio x 2 tech
# 2.Calculate fraction of each protein against total, take mean from replicates, calculate stoichiometry from mean 
# Decide to use approach 2.  Rational is to 1) reduce variation from different extraction, 2) help with unstable measurement of low abundant target such as IKZF1
# in the AP-MS target protein is usually high abundant, so calculate stoichiometry directly against bait for each run

Protein_long_abudance <- ge_hq |> 
  mutate(abudance = 2^LogIntensities) |> 
  group_by(RUN) |> 
  mutate(fraction = abudance/sum(abudance)) |> 
  ungroup()

Mean_protein_frac <- Protein_long_abudance |> 
  group_by(Gene,Protein_Acc,GROUP) |> 
  summarise(mean_frac = mean(fraction, na.rm = TRUE),
            cv_frac = sd(fraction, na.rm =TRUE)/mean_frac) |> 
  ungroup()

Lysate_sto <- Mean_protein_frac |> 
  group_by(GROUP) |> 
  mutate(Sto_IRF4 = mean_frac/mean_frac[Gene=="IRF4"],
         Sto_p300 = mean_frac/mean_frac[Gene=="EP300"],
         Sto_IKZF1 = mean_frac/mean_frac[Gene=="IKZF1"],
         Sto_IKZF3 = mean_frac/mean_frac[Gene=="IKZF3"]
  ) |> 
  ungroup()

write_csv(Lysate_sto, str_c(study_name, "_lysate_sto.csv"))










