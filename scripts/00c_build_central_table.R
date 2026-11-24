#!/usr/bin/env Rscript
# 00c_build_central_table.R
# Compiles every gene list into ONE central table that downstream scripts can read.
#
# RUN (from the project root, or source() from anywhere in the project)
#   Rscript scripts/00c_build_central_table.R
#
# IN   data/gene-lists/trophic_mode_panel.csv                         (built by 00_build_gene_sets.R)
#      data/gene-lists/protist_phagotrophy_experimental.csv           (hand-curated, one citation per row)
#      data/gene-lists/candidates/symbiosis_host_candidates_20261111.csv
#      data/gene-lists/candidates/parasitism_host_candidates_Bu2022_Dean2014_20261111.csv
#      data/gene-lists/pfam_name_map.csv
# OUT  data/gene-lists/gene_candidates_central.csv     (the central table; one row per gene x KO or Pfam feature)
#      output-tables/00c_central_table_excluded.csv    (rows left out, with the reason)
#
# This writes a NEW file next to trophic_mode_panel.csv. Nothing is overwritten and
# 06 and 07 still read the old panel until they are pointed at this table.
# TODO (reminder): compare the two tables, then slim down to one panel file.
#
# COLUMNS
#   feature_id   KO id, or "pfam:<name>" when a row has no KO (same convention as the panel)
#   match_type   ko | pfam | none
#   ko_id, pfam_acc, pfam_name
#   gene         gene or family name as written in the source list
#   category     Heterotrophy/phagotrophy | Osmotrophy | Symbiosis | Parasitism
#   subcategory  functional grouping inside the mode
#   mechanism    why the gene is included (function and the proposed mechanism)
#   citation, doi
#   notes        caveats from the source list
#   source       which list the row came from
#   status       resolved_ko | resolved_pfam | needs_review | unresolved
#                needs_review: a Pfam-matched row whose domain name is not yet in pfam_name_map
#   shared_modes a feature kept in more than one mode (e.g. "Heterotrophy/phagotrophy + Symbiosis"); shared_flag = TRUE then
#                unresolved:   no KO or Pfam yet; kept visible so it can be filtered out
#   link         direct (included) | context | dropped  (context and dropped are both excluded, and logged)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# --- locate the project root (same contract as scripts/_setup.R) ---------------
find_project_root <- function(start = getwd(),
                              markers = c("config/params.yml", "_quarto.yml")) {
  dir_now <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (any(file.exists(file.path(dir_now, markers)))) return(dir_now)
    parent <- dirname(dir_now)
    if (identical(parent, dir_now)) break
    dir_now <- parent
  }
  stop("Could not find the project root (looked for config/params.yml upward from '",
       start, "').", call. = FALSE)
}
project_root <- find_project_root()
if (!identical(normalizePath(getwd(), winslash = "/"), normalizePath(project_root, winslash = "/"))) {
  setwd(project_root)
  message("[central-table] working directory set to project root: ", project_root)
}

# --- paths -------------------------------------------------------------------
path_panel     <- "data/gene-lists/trophic_mode_panel.csv"
path_exper     <- "data/gene-lists/protist_phagotrophy_experimental.csv"
path_sym       <- "data/gene-lists/candidates/symbiosis_host_candidates_20261111.csv"
path_par       <- "data/gene-lists/candidates/parasitism_host_candidates_Bu2022_Dean2014_20261111.csv"
path_pfam_map  <- "data/gene-lists/pfam_name_map.csv"
path_central   <- "data/gene-lists/gene_candidates_central.csv"
path_excluded  <- "output-tables/00c_central_table_excluded.csv"

stopifnot(file.exists(path_panel), file.exists(path_exper), file.exists(path_sym),
          file.exists(path_par), file.exists(path_pfam_map))

MODES <- c("Heterotrophy/phagotrophy", "Osmotrophy", "Symbiosis", "Parasitism")

# --- mechanism notes written by hand -------------------------------------------
# A gene family matches a row when its name starts with the pattern. Add a line here
# to give a family a fuller "why it is included" statement than the source list's
# function text. Hypotheses are labelled as hypotheses.
mechanism_notes <- tribble(
  ~pattern, ~mechanism_note,
  "Rab GTPases",
  "Hypothesis (Jerlstrom-Hultqvist et al. 2024): the symbiont compartment resembles a phagosome whose maturation is held at an early stage and elaborated into a symbiosome instead of proceeding to lysosomal digestion. Expanded Rab families are candidates for tuning that process. Evidence is gene copy number, not expression.",
  "TBC RabGAPs",
  "Hypothesis (Jerlstrom-Hultqvist et al. 2024): TBC proteins switch Rabs off; TBC-E and TBC-K duplications predate the Anaeramoeba species split, and TBC-K regulates Arf6, which is linked to phagosome maturation.",
  "HOPS/CORVET",
  "Hypothesis (Jerlstrom-Hultqvist et al. 2024): HOPS/CORVET subunits drive Rab conversion and fusion steps; expansion of the HOPS-specific subunits may control how long a symbiont-containing vesicle stays in each maturation state."
)

# --- helpers -------------------------------------------------------------------
split_ids <- function(x) {
  # "K00533;K00534" -> c("K00533", "K00534"); NA or "" -> NA
  x <- str_trim(x)
  if (is.na(x) || x == "" || x == "NA") return(NA_character_)
  str_split_1(x, "\\s*[;,]\\s*")
}

# One candidate/experimental row -> one row per KO (or one row when it has none).
expand_kos <- function(dat, ko_col = "ko_id") {
  dat |>
    mutate(.ko = map(.data[[ko_col]], split_ids)) |>
    unnest(.ko) |>
    mutate(ko_id = .ko) |>
    select(-.ko)
}

pfam_map <- read_csv(path_pfam_map, show_col_types = FALSE)

# --- 1. existing panel (already resolved to KOs / Pfam features) ----------------
panel <- read_csv(path_panel, show_col_types = FALSE) |>
  transmute(
    feature_id, match_type, ko_id,
    pfam_acc, pfam_name,
    gene,
    category,
    subcategory,
    mechanism   = NA_character_,
    citation    = NA_character_,
    doi         = NA_character_,
    notes       = resolution_note,
    source,
    status      = if_else(match_type == "ko", "resolved_ko", "resolved_pfam"),
    link        = NA_character_
  )

# --- 2. candidate lists (symbiosis + parasitism) ---------------------------------
read_candidates <- function(path) {
  dat <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  # the two files name the link column differently
  if ("symbiosis_link" %in% names(dat)) dat <- rename(dat, link = symbiosis_link)
  # curated mechanism text (2026-11-11); kept aside so the computed fallback coalesces under it
  if ("mechanism" %in% names(dat)) dat <- rename(dat, mechanism_csv = mechanism)
  if (!"mechanism_csv" %in% names(dat)) dat$mechanism_csv <- NA_character_
  dat
}

candidates_all <- bind_rows(
  read_candidates(path_sym) |> mutate(source = "symbiosis_host_candidates_20261111"),
  read_candidates(path_par) |> mutate(source = "parasitism_host_candidates_20261111")
)

stopifnot(all(candidates_all$category_proposed %in% MODES),
          all(candidates_all$link %in% c("direct", "context", "dropped")))

# only host-parasite / host-symbiont interface genes go in; context rows are logged
excluded_context <- candidates_all |>
  filter(link != "direct") |>
  transmute(source, gene, category = category_proposed, link,
            reason = if_else(link == "dropped",
                             "dropped by panel decision (see the row's notes)",
                             "link = context (general metabolism or stress gene)"))

candidates <- candidates_all |>
  filter(link == "direct") |>
  mutate(
    # prefix match against the hand-written notes above
    mechanism_note = map_chr(gene, \(g) {
      hit <- mechanism_notes$mechanism_note[str_starts(g, fixed(mechanism_notes$pattern))]
      if (length(hit)) hit[[1]] else NA_character_
    }),
    mechanism = coalesce(
      na_if(str_squish(mechanism_csv), ""),
      if_else(!is.na(mechanism_note), paste0(`function`, ". ", mechanism_note),
              paste0(`function`, ". Evidence: ", evidence))
    ),
    n_ko = map_int(ko_id, \(x) sum(!is.na(split_ids(x))))
  ) |>
  expand_kos() |>
  mutate(
    match_type = case_when(!is.na(ko_id) ~ "ko", !is.na(pfam_acc) & pfam_acc != "" ~ "pfam", TRUE ~ "none"),
    pfam_acc   = na_if(pfam_acc, ""),
    pfam_name  = pfam_map$eggnog_name[match(pfam_acc, pfam_map$pfam_acc)],
    feature_id = case_when(
      match_type == "ko"   ~ ko_id,
      match_type == "pfam" ~ paste0("pfam:", coalesce(pfam_name, pfam_acc)),
      TRUE                 ~ NA_character_
    ),
    status = case_when(
      match_type == "none" ~ "unresolved",
      # multi-KO rows were all reviewed in the 2026-11-11 voice pass -> a kept multi-KO row is resolved
      match_type == "pfam" & is.na(pfam_name) ~ "needs_review",
      match_type == "ko"   ~ "resolved_ko",
      TRUE                 ~ "resolved_pfam"
    ),
    notes = case_when(
      match_type == "pfam" & is.na(pfam_name) ~
        paste0(na_if(notes, ""), " | Pfam accession not in pfam_name_map.csv: add the eggNOG-mapper name before 06 can match it"),
      TRUE ~ na_if(notes, "")
    ),
    category    = category_proposed,
    subcategory = subcategory_proposed,
    citation    = source_paper
  ) |>
  select(feature_id, match_type, ko_id, pfam_acc, pfam_name, gene, category, subcategory,
         mechanism, citation, doi, notes, source, status, link)

# --- 3. experimental list (loss-of-function evidence in protists) -------------------
experimental <- read_csv(path_exper, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
  expand_kos() |>
  mutate(
    match_type = case_when(!is.na(ko_id) ~ "ko", !is.na(pfam_acc) & pfam_acc != "" ~ "pfam", TRUE ~ "none"),
    pfam_acc   = na_if(pfam_acc, ""),
    pfam_name  = pfam_map$eggnog_name[match(pfam_acc, pfam_map$pfam_acc)],
    feature_id = case_when(
      match_type == "ko"   ~ ko_id,
      match_type == "pfam" ~ paste0("pfam:", coalesce(pfam_name, pfam_acc)),
      TRUE                 ~ NA_character_
    )
  ) |>
  transmute(
    feature_id, match_type, ko_id, pfam_acc, pfam_name,
    gene        = gene_label,
    category    = "Heterotrophy/phagotrophy",
    subcategory = "Experimental evidence (protist loss of function)",
    mechanism   = coalesce(na_if(str_squish(mechanism), ""), paste0(organism, ": ", evidence)),
    citation, doi,
    notes       = na_if(ko_match, ""),
    source      = "protist_phagotrophy_experimental",
    status      = case_when(match_type == "none" ~ "unresolved",
                            match_type == "ko"   ~ "resolved_ko",
                            TRUE                 ~ "resolved_pfam"),
    link        = NA_character_
  )

# --- 4. bind, check, write ----------------------------------------------------------
central <- bind_rows(panel, candidates, experimental) |>
  mutate(across(where(is.character), \(x) na_if(str_squish(x), ""))) |>
  arrange(category, subcategory, gene, feature_id)

stopifnot(
  all(central$category %in% MODES),
  all(central$status %in% c("resolved_ko", "resolved_pfam", "needs_review", "unresolved")),
  # a feature id must exist exactly when the row is not unresolved
  all(is.na(central$feature_id) == (central$status == "unresolved")),
  all(str_detect(central$ko_id[!is.na(central$ko_id)], "^K\\d{5}$"))
)

# --- shared-mode key (#27): a feature kept in >1 mode gets a combination label ---
# The row stays in each mode it belongs to; 06 plots it primarily under
# Heterotrophy/phagotrophy and marks the other mode(s) with a symbol.
shared_key <- central |>
  filter(!is.na(feature_id)) |>
  summarise(.by = feature_id,
            shared_modes = { m <- sort(unique(category))
                             if (length(m) > 1) paste(m, collapse = " + ") else NA_character_ })
central <- central |>
  left_join(shared_key, by = "feature_id") |>
  mutate(shared_flag = !is.na(shared_modes))

dir.create(dirname(path_excluded), showWarnings = FALSE, recursive = TRUE)
write_csv(central, path_central, na = "")
write_csv(excluded_context, path_excluded, na = "")

# --- report ---------------------------------------------------------------------------
message(sprintf("[central-table] wrote %s (%d rows, %d distinct features)",
                path_central, nrow(central), n_distinct(na.omit(central$feature_id))))
central |>
  count(category, status, name = "n_rows") |>
  pivot_wider(names_from = status, values_from = n_rows, values_fill = 0) |>
  print(n = Inf)

shared <- central |>
  filter(!is.na(feature_id)) |>
  summarise(n_modes = n_distinct(category), modes = paste(sort(unique(category)), collapse = " + "), .by = feature_id) |>
  filter(n_modes > 1)
if (nrow(shared)) {
  message("[central-table] ", nrow(shared), " feature(s) sit in more than one mode (downstream scripts decide how to count them):")
  print(shared, n = 20)
}
message(sprintf("[central-table] %d context-only candidate rows excluded -> %s",
                nrow(excluded_context), path_excluded))
