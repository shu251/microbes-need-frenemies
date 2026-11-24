# QC_pfam_panel_check.R
# Can the trophic-mode panel use Pfam where no KEGG KO exists?
#
# DIAGNOSTIC ONLY: changes nothing in the pipeline and writes four QC tables.
# Run from the project root after QC of the gene lists:
#   source("scripts/QC_pfam_panel_check.R")
#
# Questions (one table each):
#   1. For each candidate Pfam domain, how many eukaryote transcripts carry it,
#      and do they also carry a KO? Transcripts with no KO are removed by the
#      shared prep (R/prep_paired.R), so a Pfam-only gene can enter the panel
#      only if its transcripts are kept some other way.
#   2. Which KOs are co-assigned to each Pfam? A domain shared by many proteins
#      (e.g. SapB in prosaposin and aspartic proteases) is not a specific marker.
#   3. Which lineages carry the candidate KOs (K04641 rhodopsin, K07093 for PhoX)?
#   4. Expression totals, so a gene with almost no signal can be dropped.
#
# NOTES FOR THE PANEL (2026-11-08)
#   * Panel rule: a gene is matched by its KEGG KO; when no KO exists it is
#     matched by its Pfam domain instead.
#   * PhoX (PF05787; listed as "DUF839" in the eggNOG PFAMs column): KEGG
#     associates this domain with K07093, which ko00001 defines only as
#     "uncharacterized protein" (09194 Poorly characterized). Because that KO
#     does not identify PhoX, PhoX is matched by Pfam.
#   * Amoebapore / saposin-like lytic peptides: matched by Pfam (eggNOG names
#     SapB_1 = PF05184, SapB_2 = PF03489). KEGG has no amoebapore KO. The saposin
#     KO, K12382 (PSAP), is prosaposin (lipid metabolism), and the same domain sits
#     inside aspartic proteases (phytepsin K08245, cathepsin D K01379), so a Pfam
#     match alone is not specific. Amoebapores, granulysin and NK-lysin are small
#     single-domain saposins, so the specific set is transcripts whose ONLY Pfam
#     domains are SapB_1/SapB_2 (`saposin_only` below).
#     K06124 and K10809 were proposed for saposin B and granulysin, but in
#     ko00001 they are PHOSPHO1 (phosphoethanolamine phosphatase) and TG
#     (thyroglobulin); neither is used.
#   * The source list gives PF00827 / "Saposin_B" for amoebapore; eggNOG writes
#     no "Saposin_B" name. Confirm the accession on InterPro before citing it.
#   * Rhodopsin: matched by K04641 (bop; bacteriorhodopsin), the only
#     light-driven proton-pump KO in ko00001. K04642 (halorhodopsin, chloride
#     pump) and K23306 (brz, haloarchaeal regulator) are not used.
#   * eggNOG-mapper writes Pfam NAMES, not PF accessions, in PFAMs, so matching
#     uses the names in `pfam_candidates` below.

source("scripts/_setup.R")

fields_in_scope <- c("Gorda Ridge", "Von Damm", "Piccard")

# Pfam name (as eggNOG writes it) -> panel gene it would stand for
pfam_candidates <- tribble(
  ~pfam_name,       ~pfam_acc,  ~panel_gene,
  "DUF839",         "PF05787",  "phoX (alkaline phosphatase PhoX)",
  "PhoX",           "PF05787",  "phoX (alkaline phosphatase PhoX)",
  "SapB_1",         "PF05184",  "amoebapore (saposin-like)",
  "SapB_2",         "PF03489",  "amoebapore (saposin-like)",
  "Bac_rhodopsin",  "PF01036",  "microbial rhodopsin"
)

ko_candidates <- tribble(
  ~ko_id,   ~panel_gene,
  "K04641", "microbial rhodopsin (bop)",
  "K07093", "phoX? (KEGG: uncharacterized protein)",
  "K08515", "VAMP7",
  "K04409", "PAK1"
)

# --- 1. load and apply the shared prep filters EXCEPT the KO requirement -------
long <- load_long_annot()

# same lineage rule as prep_paired(), computed once per distinct lineage
lin_ok <- long |>
  distinct(full_classification) |>
  filter(!is.na(full_classification)) |>
  mutate(.l2 = lineage_depth(full_classification, 2L)) |>
  filter(!is.na(.l2), .l2 != "Eukaryota",
         !grepl("Opisthokonta", full_classification) |
           grepl("Choanoflagell", full_classification)) |>
  pull(full_classification)

euk <- long |>
  filter(field %in% fields_in_scope, full_classification %in% lin_ok) |>
  mutate(has_ko = !is.na(ko_id) & !ko_id %in% c("", "-"),
         PFAMs  = replace_na(PFAMs, "")) |>
  add_curated_taxonomy()

rm(long); invisible(gc())

# --- 2. Pfam candidates: transcripts, KO status, expression -------------------
pfam_hits <- pfam_candidates |>
  mutate(rows = map(pfam_name, \(p) euk |>
                      filter(str_detect(PFAMs, paste0("(^|,)", p, "(,|$)"))) |>
                      select(SEQ_ID, ko_id, has_ko, PFAMs, Supergroup, Phylum,
                             cur_group, scaledTPM))) |>
  unnest(rows) |>
  # TRUE when every Pfam on the transcript is one of this gene's candidate domains
  mutate(domain_only = map2_lgl(PFAMs, panel_gene, \(p, g) {
    all(str_split_1(p, ",") %in% pfam_candidates$pfam_name[pfam_candidates$panel_gene == g])
  }))

pfam_summary <- pfam_hits |>
  summarise(n_transcripts       = n_distinct(SEQ_ID),
            n_transcripts_no_ko = n_distinct(SEQ_ID[!has_ko]),
            tpm_total           = sum(scaledTPM),
            tpm_no_ko           = sum(scaledTPM[!has_ko]),
            n_transcripts_domain_only = n_distinct(SEQ_ID[domain_only]),
            tpm_domain_only     = sum(scaledTPM[domain_only]),
            .by = c(panel_gene, pfam_name, pfam_acc)) |>
  mutate(frac_tpm_lost_to_prep = tpm_no_ko / tpm_total)

pfam_summary

# --- 3. which KOs share each Pfam (specificity) --------------------------------
pfam_co_ko <- pfam_hits |>
  filter(has_ko) |>
  separate_longer_delim(ko_id, ",") |>
  distinct(panel_gene, pfam_name, SEQ_ID, ko_id) |>
  count(panel_gene, pfam_name, ko_id, name = "n_transcripts", sort = TRUE)

pfam_co_ko

# --- 4. KO candidates: lineages and Pfam names -------------------------------
ko_hits <- ko_candidates |>
  mutate(rows = map(ko_id, \(k) euk |>
                      filter(str_detect(ko_id, k)) |>
                      select(SEQ_ID, PFAMs, Supergroup, Phylum, cur_group,
                             scaledTPM))) |>
  select(-ko_id) |>
  unnest(rows)

ko_lineages <- ko_hits |>
  summarise(n_transcripts = n_distinct(SEQ_ID), tpm_total = sum(scaledTPM),
            pfams = paste(sort(unique(PFAMs)), collapse = " | "),
            .by = c(panel_gene, cur_group)) |>
  arrange(panel_gene, desc(tpm_total))

ko_lineages

# --- write --------------------------------------------------------------------
write_csv(pfam_summary, "output-tables/qc-pfam-panel-summary.csv")
write_csv(pfam_co_ko,   "output-tables/qc-pfam-panel-co-ko.csv")
write_csv(ko_lineages,  "output-tables/qc-pfam-panel-ko-lineages.csv")
write_csv(pfam_hits |> count(panel_gene, pfam_name, cur_group, name = "n_rows"),
          "output-tables/qc-pfam-panel-lineages.csv")

message("Wrote output-tables/qc-pfam-panel-*.csv")
