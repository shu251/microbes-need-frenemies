# Compile large datasets for use in R. 
library(data.table)
library(tximport)
library(readr)
library(tidyverse)

cat("Importing files\n")
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/dfs_gr_mcr_feb2025.RData", verbose = TRUE)
# lenscaled_TPM, gr_mcr_lenscaled_TPM_df, mean_gr_mcr_TPM_df

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

tx2gene_in <- taxfxn %>% 
  separate_longer_delim(KEGG_ko, delim = ",") %>% # Remove comma separated KEGG IDs
  select(SequenceID, seed_ortholog, COG_category, GOs, KEGG = KEGG_ko, PFAMs, classification_level, full_classification, classification, max_pid) %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))

cat("Modified taxfxn:\n")
head(tx2gene_in)
# head(rownames(mean_gr_mcr_TPM_df))

cat("Files imported.\n")

cat("View mean_gr_mcr_TPM_df\n\n")
head(mean_gr_mcr_TPM_df)

cat("Pivot to long format:\n")
# Make longer data frame
long_frenemies <- mean_gr_mcr_TPM_df %>%
  rowwise() %>%
  # Add column with number of samples with a 0.
  mutate(NUM_ZERO = sum(c_across(starts_with("mean.")) == 0)) %>%
  rownames_to_column(var = "SEQ_ID") %>%
  pivot_longer(cols = starts_with("mean"), values_to = "scaledTPM") %>%
  # Remove zeroes
  filter(scaledTPM > 0) %>%
  separate(name, c("meanfield", "LIBRARY_NUM", "fieldyear", "LOCATION", "SAMPLETYPE", "SAMPLEID"), "_",
           remove = FALSE) %>%
  # Isolate the in situ samples only
  filter(SAMPLETYPE == "insitu") %>% 
  select(-meanfield, -name) %>%
  mutate(VENT_FIELD = case_when(grepl("Piccard", fieldyear) ~ "Piccard",
                                grepl("VonDamm", fieldyear) ~ "Von Damm",
                                grepl("Axial", fieldyear) ~ "Axial",
                                grepl("Gorda", fieldyear) ~ "Gorda Ridge")) %>%
  mutate(VENT_BIN = case_when(
    (LOCATION == "Background" | LOCATION == "Plume" | LOCATION == "BSW") ~ "Non-vent",
    grepl("IntlDistrict", LOCATION) ~ "Non-vent",
    grepl("ASHES", LOCATION) ~ "Non-vent",
    TRUE ~ "Vent"
  ))



# Add in annotations, revise taxa
as_is <- c("Amoebozoa", "Apusozoa", "Excavata", "Hacrobia", "Archaeplastida")

# head(mean_gr_mcr_TPM_df)
# head(taxfxn)

long_frenemies_annot <- long_frenemies %>%
  left_join(tx2gene_in, by = "SEQ_ID") %>%
  separate(full_classification, c("Domain", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus_spp"), sep = "; ", remove = FALSE) %>%
  mutate(SUPERGROUP_18S = case_when(
    Phylum == "Ciliophora" ~ "Alveolata-Ciliophora",
    Phylum == "Dinophyta" ~ "Alveolata-Dinoflagellata",
    # Phylum == "Perkinsea" ~ "Protalveolata",
    # Phylum == "Colponemidia" ~ "Protalveolata",
    # Phylum == "Chromerida" ~ "Protalveolata",
    Supergroup == "Alveolata" ~ "Other Alveolata",
    Supergroup %in% as_is ~ Supergroup,
    Supergroup == "Haptista" ~ "Hacrobia",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    (Supergroup == "Rhizaria" & Phylum != "Radiolaria" & Phylum != "Cercozoa") ~ "Rhizaria",
    Order == "Bigyra" ~ "Stramenopiles-Opalozoa;Sagenista",
    Class == "Ochromonadales" ~ "Stramenopiles-Ochrophyta",
    Supergroup == "Stramenopiles" ~ "Stramenopiles",
    Supergroup == "Opisthokonta" ~ "Opisthokonta",
    (is.na(Supergroup) | Supergroup == "Eukaryota incertae sedis") ~ "Unknown Eukaryota",
    TRUE ~ "Other-metaT only"))

cat("Completed data anneal, output:\n")

head(long_frenemies_annot)

# Save output R object files
save(long_frenemies, long_frenemies_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/frenemies-longdfs-gr-mcr.RData")
cat("SAVED\n")

cat("switching to making matrices for use:")
# head(mean_gr_mcr_TPM_df)
# head(long_frenemies)
matrix_frenemies <- long_frenemies %>% 
  unite(SAMPLE_NAME, VENT_FIELD, VENT_BIN, LOCATION, sep = "_") %>% 
  select(SEQ_ID, SAMPLE_NAME, scaledTPM, COG_category, GOs, KEGG_ko, PFAMs, full_classification) %>% 
  distinct() %>% 
  group_by(SEQ_ID, SAMPLE_NAME) %>% 
    summarise(SUM_TPM = sum(scaledTPM)) %>% 
  pivot_wider(names_from = SAMPLE_NAME, values_from = SUM_TPM, values_fill = 0) %>% 
  as.matrix

matrix_frenemies_ventonly <- long_frenemies %>% 
  unite(SAMPLE_NAME, VENT_FIELD, VENT_BIN, LOCATION, sep = "_") %>% 
  filter(VENT_BIN = "Vent") %>% 
  select(SEQ_ID, SAMPLE_NAME, scaledTPM, COG_category, GOs, KEGG_ko, PFAMs, full_classification) %>% 
  distinct() %>% 
  group_by(SEQ_ID, SAMPLE_NAME) %>% 
  summarise(SUM_TPM = sum(scaledTPM)) %>% 
  pivot_wider(names_from = SAMPLE_NAME, values_from = SUM_TPM, values_fill = 0) %>% 
  as.matrix

matrix_frenemies_ciliates <- long_frenemies_annot %>%
  filter(SUPERGROUP_18S == "Alveolata-Ciliophora") %>% 
  unite(SAMPLE_NAME, VENT_FIELD, VENT_BIN, LOCATION, sep = "_") %>% 
  select(SEQ_ID, SAMPLE_NAME, scaledTPM, COG_category, GOs, KEGG_ko, PFAMs, full_classification, SUPERGROUP_18S) %>% 
  distinct() %>% 
  group_by(SEQ_ID, SAMPLE_NAME) %>% 
    summarise(SUM_TPM = sum(scaledTPM)) %>% 
  pivot_wider(names_from = SAMPLE_NAME, values_from = SUM_TPM, values_fill = 0) %>% 
  as.matrix

matrix_frenemies_ciliates_ventonly <- long_frenemies_annot %>%
  filter(SUPERGROUP_18S == "Alveolata-Ciliophora") %>% 
  filter(VENT_BIN = "Vent") %>% 
  unite(SAMPLE_NAME, VENT_FIELD, VENT_BIN, LOCATION, sep = "_") %>% 
  select(SEQ_ID, SAMPLE_NAME, scaledTPM, COG_category, GOs, KEGG_ko, PFAMs, full_classification, SUPERGROUP_18S) %>% 
  distinct() %>% 
  group_by(SEQ_ID, SAMPLE_NAME) %>% 
  summarise(SUM_TPM = sum(scaledTPM)) %>% 
  pivot_wider(names_from = SAMPLE_NAME, values_from = SUM_TPM, values_fill = 0) %>% 
  as.matrix

# Ciliate taxa key
# head(long_frenemies_annot)
key_taxa <- long_frenemies_annot %>% 
  select(SEQ_ID, full_classification, SUPERGROUP_18S) %>% 
  distinct() %>% 
  separate(full_classification, c("Domain", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus_spp"), sep = "; ", remove = FALSE)

cat("Saving data files:\n")
save(matrix_frenemies, matrix_frenemies_ventonly, key_taxa, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/matrix_frenemies.RData")

save(matrix_frenemies_ciliates, matrix_frenemies_ciliates_ventonly, key_taxa, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/matrix_ciliates_frenemies.RData")
  
