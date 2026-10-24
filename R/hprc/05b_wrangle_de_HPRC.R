## =============================================================================
## 05b — Reduce and classify the ALDEx2 results on HPRC, then ship a small
##       bundle to the laptop
## Reads the outputs of 05a_aldex_lean_HPRC.R. Never touches the tximport bundle.
## =============================================================================
##
## WHAT THIS IS FOR
## ----------------
## 05a turns the Monte Carlo whale into a few small tables. This script does the
## wrangling that is cheapest where the data already lives, so that what crosses
## to the laptop is a curated set rather than a full transcriptome:
##
##   PART 1 (always runs, needs no annotation)
##     - per-coefficient significance and effect-size tiers
##     - a CROSS-AXIS classification: is a transcript a vent responder, a field
##       responder, both, or neither? This is the "across my various axes"
##       question, and it is a pivot, not a new model -- 05a already estimated
##       every coefficient in one fit.
##     - a candidate set, and the CLR medians restricted to it
##
##   PART 2 (optional, runs only if a curated-gene CSV is present on HPRC)
##     - join annotation, aggregate to KO, and re-run with a housekeeping-
##       informed scale model
##
## PART 2 is opt-in because your curated genes live on the laptop. If you never
## upload them, PART 1 still produces everything you need to do the annotation
## join locally -- the DE tables are keyed on transcript_id and are small.
##
##
## A NOTE ON WHICH THRESHOLD TO TRUST
## ----------------------------------
## With ~23 libraries and a transcriptome-sized feature set, BH-adjusted p-values
## are the wrong primary filter -- the earlier unfiltered run returning n_sig = 3
## of 207k is what that looks like. ALDEx2's own recommendation for small n is to
## rank and threshold on the standardised effect size (|effect| > 1), using the
## adjusted p-value as a secondary filter rather than the gate
## (Fernandes et al. 2014, Microbiome 2:15, doi:10.1186/2049-2618-2-15).
##
## Both tiers are computed here and kept side by side. Report the effect-size
## tier as primary and say so in the methods; do not quietly switch between them.
## =============================================================================

library(tidyverse)


## ---------------------------------------------------------------------------
## Paths and knobs
## ---------------------------------------------------------------------------
base_dir <- "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output"
in_dir   <- file.path(base_dir, "05a-aldex-lean")
out_dir  <- file.path(base_dir, "05b-for-local")

# PART 2 input. Upload your curated gene table here to switch PART 2 on.
# Must contain a transcript-level join key and, ideally, ko_id + a taxon column.
annotation_path <- file.path(base_dir, "curated-genes.csv")
annotation_key  <- "transcript_id"   # column in the annotation that matches ours

effect_threshold <- 1.0   # |effect| > 1 -- the ALDEx2 reporting convention
fdr_alpha        <- 0.05  # overwritten from the manifest if present

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------------------------------------------
## Load 05a outputs
## ---------------------------------------------------------------------------
manifest <- readRDS(file.path(in_dir, "run_manifest.rds"))
# Not using %||% here: it is only in base R from 4.4, and HPRC modules lag.
if (!is.null(manifest$fdr_alpha)) fdr_alpha <- manifest$fdr_alpha

message(glue::glue(
  "05a run: {manifest$n_transcripts} transcripts x {manifest$n_libraries} ",
  "libraries, mc.samples = {manifest$mc_samples}, gamma = {manifest$gamma_main}, ",
  "ALDEx2 {manifest$aldex2_version}."))

de_all_terms <- readr::read_csv(file.path(in_dir, "de_all_terms.csv.gz"),
                                show_col_types = FALSE)
clr_median   <- readRDS(file.path(in_dir, "clr_median.rds"))
design       <- readr::read_csv(file.path(in_dir, "design.csv"),
                                show_col_types = FALSE)

effect_path <- file.path(in_dir, "effect_sizes.csv.gz")
has_effect  <- file.exists(effect_path)

effect_sizes <- if (has_effect) {
  readr::read_csv(effect_path, show_col_types = FALSE)
} else {
  message("No effect_sizes.csv.gz -- effect-size tiers will be skipped.")
  NULL
}


## =============================================================================
## PART 1 — annotation-free reduction and cross-axis classification
## =============================================================================

## ---------------------------------------------------------------------------
## Tier each transcript within each coefficient
## ---------------------------------------------------------------------------
# Attach effect sizes if 05a produced them, otherwise carry NA columns so the
# rest of the script has a stable shape either way. Written as a plain if/else
# rather than an inline lambda in the pipe -- easier to read and to step through.
if (has_effect) {
  de_with_effect <- de_all_terms |>
    left_join(effect_sizes |> select(transcript_id, term, effect, overlap),
              by = c("transcript_id", "term"))
} else {
  de_with_effect <- de_all_terms |>
    mutate(effect = NA_real_, overlap = NA_real_)
}

de_tiered <- de_with_effect |>
  mutate(
    direction = case_when(
      estimate > 0 ~ "up",
      estimate < 0 ~ "down",
      TRUE         ~ "flat"),
    # Two independent calls, deliberately not collapsed into one flag.
    tier_fdr    = !is.na(p_adj)  & p_adj < fdr_alpha,
    tier_effect = !is.na(effect) & abs(effect) > effect_threshold,
    tier = case_when(
      tier_fdr & tier_effect ~ "both",
      tier_effect            ~ "effect_only",
      tier_fdr               ~ "fdr_only",
      TRUE                   ~ "none"))

message("Transcripts per coefficient per tier:")
de_tiered |>
  count(term, tier) |>
  pivot_wider(names_from = tier, values_from = n, values_fill = 0) |>
  print(n = Inf)


## ---------------------------------------------------------------------------
## Cross-axis classification — one row per transcript
## ---------------------------------------------------------------------------
# Which axis does a transcript respond to? The vent coefficient and the field
# coefficients came from the SAME fit, so these are marginal contrasts and can
# be read side by side. A transcript significant on a field coefficient but not
# on vent_bin is a field responder; the reverse is a vent responder; both is
# a transcript whose expression tracks the vent boundary AND differs among
# fields on top of that.
#
# Interpretation caveat carried down from Q1: among-field structure was NS at
# the community level, so a long list of field responders here would be
# individually-detectable genes inside a community-level null, not a
# contradiction of Q1. Say that explicitly if you report it.
vent_terms  <- grep("^vent_bin", manifest$coefficients, value = TRUE)
field_terms <- grep("^field",    manifest$coefficients, value = TRUE)

message(glue::glue("Vent coefficients: {paste(vent_terms, collapse = ', ')}"))
message(glue::glue("Field coefficients: {paste(field_terms, collapse = ', ')}"))

responds <- function(tbl, which_terms) {
  tbl |>
    filter(term %in% which_terms, tier %in% c("both", "effect_only", "fdr_only")) |>
    distinct(transcript_id) |>
    pull(transcript_id)
}

vent_responders  <- responds(de_tiered, vent_terms)
field_responders <- responds(de_tiered, field_terms)

cross_axis <- de_tiered |>
  distinct(transcript_id) |>
  mutate(
    is_vent_responder  = transcript_id %in% vent_responders,
    is_field_responder = transcript_id %in% field_responders,
    axis_class = case_when(
      is_vent_responder &  is_field_responder ~ "both_axes",
      is_vent_responder & !is_field_responder ~ "vent_only",
      !is_vent_responder & is_field_responder ~ "field_only",
      TRUE                                    ~ "no_axis"))

message("Cross-axis classification:")
cross_axis |> count(axis_class) |> arrange(desc(n)) |> print(n = Inf)

# Vent direction, attached where it exists -- the column you will most want
# when the trophic-mode enrichment runs locally.
vent_direction <- de_tiered |>
  filter(term %in% vent_terms) |>
  group_by(transcript_id) |>
  slice_min(p_adj, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(transcript_id,
         vent_estimate = estimate,
         vent_effect   = effect,
         vent_p_adj    = p_adj,
         vent_direction = direction)

cross_axis <- cross_axis |> left_join(vent_direction, by = "transcript_id")


## ---------------------------------------------------------------------------
## Candidate set — what actually needs to travel
## ---------------------------------------------------------------------------
candidate_ids <- cross_axis |>
  filter(axis_class != "no_axis") |>
  pull(transcript_id)

message(glue::glue(
  "Candidate set: {length(candidate_ids)} of {nrow(cross_axis)} transcripts ",
  "({round(100 * length(candidate_ids) / nrow(cross_axis), 1)}%)."))

clr_candidates <- clr_median[rownames(clr_median) %in% candidate_ids, , drop = FALSE]

message(glue::glue(
  "CLR medians for candidates: {nrow(clr_candidates)} x {ncol(clr_candidates)} ",
  "({round(object.size(clr_candidates) / 1e6, 1)} MB)."))


## ---------------------------------------------------------------------------
## Library-level sanity check before anything leaves the cluster
## ---------------------------------------------------------------------------
# CLR values are centred per library by construction, so every column mean
# should sit at ~0. A column that does not is a sign the matrix was subset or
# reordered incorrectly somewhere upstream.
library_check <- tibble(
    sample_rep = colnames(clr_median),
    clr_mean   = colMeans(clr_median),
    clr_sd     = apply(clr_median, 2, sd)) |>
  left_join(design |> select(sample_rep, field, vent_bin), by = "sample_rep")

message("Per-library CLR summary (clr_mean should be ~0):")
print(library_check, n = Inf)

if (max(abs(library_check$clr_mean)) > 1e-6) {
  warning("Some library CLR means are not ~0. Check the subsetting in 05a.")
}


## ---------------------------------------------------------------------------
## Write the PART 1 bundle
## ---------------------------------------------------------------------------
readr::write_csv(de_tiered,   file.path(out_dir, "de_tiered_all_terms.csv.gz"))
readr::write_csv(cross_axis,  file.path(out_dir, "cross_axis_classification.csv.gz"))
readr::write_csv(library_check, file.path(out_dir, "library_clr_check.csv"))
saveRDS(clr_candidates, file.path(out_dir, "clr_median_candidates.rds"),
        compress = "xz")
saveRDS(manifest, file.path(out_dir, "run_manifest.rds"))


## =============================================================================
## PART 2 — optional: annotation join, KO aggregation, housekeeping scale model
## =============================================================================
## Runs only if you have scp'd your curated gene table to `annotation_path`.
## Everything above is complete without it.

if (!file.exists(annotation_path)) {

  message(glue::glue(
    "\nPART 2 skipped: no annotation at {annotation_path}.\n",
    "Everything in {out_dir} is keyed on transcript_id, so the join can be ",
    "done on the laptop instead. Upload the curated gene CSV here if you want ",
    "the KO aggregation and the housekeeping scale model to run on HPRC."))

} else {

  message(glue::glue("PART 2: joining annotation from {annotation_path}"))

  annotation <- readr::read_csv(annotation_path, show_col_types = FALSE) |>
    rename(transcript_id = all_of(annotation_key))

  # Coverage gate, in the spirit of load_annotation()'s 0.60 floor: a silent
  # low-coverage join is worse than a stop.
  coverage <- mean(unique(de_tiered$transcript_id) %in% annotation$transcript_id)
  message(glue::glue("Annotation covers {round(100 * coverage, 1)}% of tested transcripts."))
  if (coverage < 0.60) {
    warning("Join coverage below 0.60 -- check that annotation_key names the ",
            "right column and that both sides use the same transcript ids.")
  }

  de_annotated <- de_tiered |> left_join(annotation, by = "transcript_id")

  readr::write_csv(de_annotated,
                   file.path(out_dir, "de_tiered_annotated.csv.gz"))

  ## -- KO-level aggregation -------------------------------------------------
  # MULTI-KO RULE (locked earlier in this project): counts cannot be split
  # fractionally, so this does NOT divide anything by the number of KOs, and it
  # does NOT duplicate a transcript's counts across its KOs. The aggregation
  # below is over TEST RESULTS, not over counts -- each transcript contributes
  # its own effect size to each KO set it belongs to, which is exactly the
  # gene-set framing. Do not read the n_transcripts column as independent.
  if ("ko_id" %in% names(de_annotated)) {

    ko_summary <- de_annotated |>
      filter(!is.na(ko_id), ko_id != "") |>
      separate_longer_delim(ko_id, delim = ",") |>
      mutate(ko_id = str_trim(ko_id)) |>
      group_by(term, ko_id) |>
      summarise(
        n_transcripts = n(),
        n_sig_fdr     = sum(tier_fdr, na.rm = TRUE),
        n_sig_effect  = sum(tier_effect, na.rm = TRUE),
        median_effect = median(effect, na.rm = TRUE),
        median_est    = median(estimate, na.rm = TRUE),
        .groups = "drop") |>
      arrange(term, desc(abs(median_effect)))

    readr::write_csv(ko_summary, file.path(out_dir, "ko_summary.csv.gz"))
    message(glue::glue("KO summary: {nrow(ko_summary)} term x KO rows."))

  } else {
    message("No ko_id column in the annotation -- KO aggregation skipped.")
  }

  ## -- Housekeeping-informed scale model ------------------------------------
  # The blind gamma sweep in 05a asks "how much would my call move if the
  # system scale were uncertain by gamma?". A better question, when you have
  # annotation, is "what does the scale actually look like?" -- estimated from
  # genes that should NOT be differentially expressed. Gloor, Nixon & Silverman
  # demonstrate exactly this for metatranscriptomes: build the scale model from
  # known housekeeping genes rather than assuming a normalisation
  # (2025, NAR Genom Bioinform 7(3):lqaf108, doi:10.1093/nargab/lqaf108).
  #
  # This re-runs aldex.clr once and so needs the same memory budget as 05a.
  # It is written as a separate block so you can submit it as its own job.
  housekeeping_ko <- c(
    # KEGG ribosomal proteins and core translation machinery. Replace this seed
    # list with a proper KEGG module join before submission.
    "K02863", "K02886", "K02906", "K02931", "K02933", "K02946",
    "K02950", "K02952", "K02967", "K02982", "K02988", "K02994")

  if ("ko_id" %in% names(de_annotated)) {

    housekeeping_ids <- de_annotated |>
      filter(!is.na(ko_id)) |>
      separate_longer_delim(ko_id, delim = ",") |>
      mutate(ko_id = str_trim(ko_id)) |>
      filter(ko_id %in% housekeeping_ko) |>
      distinct(transcript_id) |>
      pull(transcript_id)

    message(glue::glue(
      "Housekeeping transcripts identified: {length(housekeeping_ids)}."))

    if (length(housekeeping_ids) < 20) {
      message("Too few housekeeping transcripts to anchor a scale model; ",
              "skipping. Widen housekeeping_ko first.")
    } else {
      message("Housekeeping set is usable. Run the scale-model pass as its own ",
              "job: load counts.RData, build the scale matrix with ",
              "ALDEx2::aldex.makeScaleMatrix() anchored on these ids, and pass ",
              "it to aldex.clr(gamma = <that matrix>).")
      saveRDS(housekeeping_ids,
              file.path(out_dir, "housekeeping_transcript_ids.rds"))
    }
  }
}


## ---------------------------------------------------------------------------
## What to copy down
## ---------------------------------------------------------------------------
message("\nBundle ready. Copy this whole folder to the laptop:")
tibble(file = list.files(out_dir, full.names = TRUE)) |>
  mutate(mb = round(file.size(file) / 1e6, 2), file = basename(file)) |>
  arrange(desc(mb)) |>
  print(n = Inf)

message(glue::glue(
  "Total: {round(sum(file.size(list.files(out_dir, full.names = TRUE))) / 1e6, 1)} MB"))

message(glue::glue(
  "\n  rsync -avP <netid>@<hprc>:{out_dir}/ ",
  "~/Desktop/Projects/microbes-need-frenemies/data/processed/05-aldex/"))

sessionInfo()
