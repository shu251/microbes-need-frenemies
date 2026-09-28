#!/usr/bin/env Rscript
# 00_build_het_gene_reference.R
# RUN
#   Rscript scripts/00_build_het_gene_reference.R
#   # or, from an R console anywhere in the project:
#   source("scripts/00_build_het_gene_reference.R")
#
# INPUTS  (curated by hand; this script never edits them)
#   data/gene-lists/marker_panel_expansion_candidates.csv
#   data/gene-lists/phagotrophy_genes_protists_1.csv
#
# OUTPUTS: data/gene-lists/het_gene_reference.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# --- locate the project root, independent of the current working directory ----
# Same contract as scripts/_setup.R: root = nearest ancestor holding
# config/params.yml or _quarto.yml. Every path below is root-relative, so this
# script behaves identically from Rscript, from RStudio, and under quarto render.
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
if (!identical(normalizePath(getwd(), winslash = "/"),
               normalizePath(project_root, winslash = "/"))) {
  setwd(project_root)
  message("[het-reference] working directory set to project root: ", project_root)
}

# --- paths -------------------------------------------------------------------
path_expansion   <- "data/gene-lists/marker_panel_expansion_candidates.csv"
path_phagotrophy <- "data/gene-lists/phagotrophy_genes_protists_1.csv"
path_reference   <- "data/gene-lists/het_gene_reference.csv"
path_unresolved  <- "output-tables/00_het_gene_reference_unresolved.csv"

stopifnot(file.exists(path_expansion), file.exists(path_phagotrophy))

# --- guard: malformed CSV quoting ---------------------------------------------
# An unquoted comma inside a field (e.g. a gene named "PI(4,5)P2 3-kinase" or
# "Glucan 1,4-alpha-glucosidase") shifts every later field right by one, so
# KEGG_KO silently receives a gene-name fragment and that KO disappears from the
# panel. Two such rows were found on 2026-09-26 and fixed by quoting the field.
# readr records these as parsing problems; stop rather than build a quietly wrong
# reference table. Fix by wrapping the offending field in double quotes.
read_panel_checked <- function(path) {
  dat <- read_csv(path, show_col_types = FALSE)
  bad <- readr::problems(dat)
  if (nrow(bad)) {
    print(bad)
    stop("'", path, "' has ", nrow(bad), " malformed row(s) -- almost always an ",
         "unquoted comma inside a field, which shifts KEGG_KO and drops that KO. ",
         "Quote the offending field in the CSV, then re-run.", call. = FALSE)
  }
  dat
}

# --- resolve one raw KEGG field into zero or more KO ids ---------------------
# `KEGG_KO` in the phagotrophy list is not always a clean id. Alongside plain
# ids and delimited lists it holds "K04624-K08424 (family range)",
# "K21348 (representative)", and one row that is not a KO at all
# ("Not in KEGG KO (non-orthologous)"). Each pattern is handled explicitly and
# every row keeps a note, so nothing disappears without a record.
#
# On family ranges: the two boundary ids are kept as members and the range is
# NOT expanded. KEGG ortholog numbers are assigned roughly sequentially rather
# than by family, so K04624-K08424 spans ~3,800 unrelated orthologs; expanding
# it would inject thousands of irrelevant KOs. The unresolved log records the
# range so real family membership can come from a KEGG lookup instead.
resolve_kegg_field <- function(raw, source_row) {
  raw <- str_trim(raw)

  if (is.na(raw) || raw == "") {
    return(tibble(source_row = source_row, ko_id = NA_character_,
                  resolution_note = "empty KEGG field"))
  }

  if (str_detect(raw, regex("not in kegg ko", ignore_case = TRUE))) {
    return(tibble(source_row = source_row, ko_id = NA_character_,
                  resolution_note = raw))
  }

  range_match <- str_match(raw, "^(K\\d{5})\\s*-\\s*(K\\d{5})\\s*\\(family range\\)$")
  if (!is.na(range_match[1, 1])) {
    return(tibble(source_row = source_row,
                  ko_id = c(range_match[1, 2], range_match[1, 3]),
                  resolution_note = paste0(
                    "family_range_boundary (range '", raw,
                    "' not expanded; boundary ids kept, full family needs a KEGG lookup)")))
  }

  rep_match <- str_match(raw, "^(K\\d{5})\\s*\\(representative\\)$")
  if (!is.na(rep_match[1, 1])) {
    return(tibble(source_row = source_row, ko_id = rep_match[1, 2],
                  resolution_note = "representative id for the family"))
  }

  tokens  <- str_trim(str_split(raw, "[,;/]+")[[1]])
  tokens  <- tokens[tokens != ""]
  valid   <- tokens[str_detect(tokens, "^K\\d{5}$")]
  invalid <- setdiff(tokens, valid)

  out <- tibble(source_row = integer(), ko_id = character(), resolution_note = character())
  if (length(valid)) {
    out <- bind_rows(out, tibble(source_row = source_row, ko_id = valid,
                                 resolution_note = NA_character_))
  }
  if (length(invalid)) {
    out <- bind_rows(out, tibble(
      source_row = source_row, ko_id = NA_character_,
      resolution_note = paste("unresolved fragment:", paste(invalid, collapse = " / "))))
  }
  if (!nrow(out)) {
    out <- tibble(source_row = source_row, ko_id = NA_character_,
                  resolution_note = paste("unresolved:", raw))
  }
  out
}

# --- subcategory canonicalization ---------------------------------------------
# The two curated lists name the same concept at different granularity and in
# different case, so the same biology would become two separate fgsea gene sets in
# 07 ("Digestive Hydrolases" with 22 KOs and "digestive hydrolase" with 7). This
# map collapses approved pairs onto one label. `subcategory_raw` keeps the source
# text, so a merge is always traceable.
#
# Approved 2026-09-26. Add a line here to merge another pair -- that is the only
# place this decision lives.
SUBCATEGORY_MERGES <- c(
  "digestive hydrolase" = "Digestive Hydrolases",
  "vesicle trafficking" = "Vesicle Trafficking & Phagosome Maturation"
)

# Candidate pairs NOT merged -- same concept at different granularity, but the call
# is Sarah's. Uncomment a line to fold it in; nothing else needs changing.
#   "acidification"               = "Vacuole / Phagolysosome Acidification",   #  6 ->  6
#   "actin machinery"             = "Adhesion & Cytoskeletal Remodeling",      #  5 -> 16
#   "prey recognition"            = "Prey Recognition & Chemosensation",       #  3 -> 12
#
# Deliberately left separate (distinct biology, not a naming variant):
#   "oxidative burst" vs "Antimicrobial / Prey Killing"
#   "phosphoinositide conversion" vs "Engulfment Signal Transduction"

# --- read both lists onto one schema -----------------------------------------
# The two lists name the same things differently: `KEGG ID` vs `KEGG_KO`,
# CATEGORY_00 / CATEGORY_01 vs a single `Category`. Only the phagotrophy list
# carries KEGG_Pathway, so expansion-candidate rows get NA there.
panel_long <- bind_rows(
  read_panel_checked(path_expansion) |>
    transmute(kegg_raw     = `KEGG ID`,
              category     = CATEGORY_00,
              subcategory  = CATEGORY_01,
              gene         = `short gene name`,
              kegg_pathway = NA_character_,
              source       = "expansion_candidates"),
  read_panel_checked(path_phagotrophy) |>
    transmute(kegg_raw     = KEGG_KO,
              category     = "Phagotrophy",
              subcategory  = Category,
              gene         = Gene_Symbol,
              kegg_pathway = KEGG_Pathway,
              source       = "phagotrophy_protists")) |>
  mutate(source_row = row_number()) |>
  # keep the source text, then canonicalize
  mutate(subcategory_raw = subcategory,
         subcategory = coalesce(SUBCATEGORY_MERGES[subcategory], subcategory))

merged_report <- panel_long |>
  filter(subcategory != subcategory_raw) |>
  count(subcategory_raw, subcategory, name = "n_source_rows")
if (nrow(merged_report)) {
  message("[het-reference] subcategory merges applied:")
  print(merged_report)
}

resolved <- map2(panel_long$kegg_raw, panel_long$source_row, resolve_kegg_field) |>
  list_rbind()

panel_resolved <- panel_long |>
  select(-kegg_raw) |>
  left_join(resolved, by = "source_row", relationship = "one-to-many")

# --- unresolved rows: write the log before filtering them out ----------------
unresolved <- panel_resolved |>
  filter(is.na(ko_id)) |>
  select(source_row, source, category, subcategory, subcategory_raw, gene, resolution_note) |>
  arrange(source, source_row)

write_csv(unresolved, path_unresolved)

# --- KEGG pathways: free text, sometimes several per row ---------------------
# e.g. "ko04145 Phagosome; ko04010 MAPK signaling". Split on ";", then pull the
# ko##### id out. The id is NOT always leading ("Nitrogen metabolism ko00910"),
# so it is matched anywhere in the field and the remainder becomes the name.
# Rows with no parsable pathway keep one row with NA pathway columns, so their KO
# stays in the reference instead of being dropped along with the pathway.
pathway_parsed <- panel_resolved |>
  filter(!is.na(ko_id)) |>
  separate_longer_delim(kegg_pathway, delim = ";") |>
  mutate(kegg_pathway = str_trim(kegg_pathway),
         pathway_id   = str_extract(kegg_pathway, "ko\\d{5}"),
         pathway_name = if_else(is.na(pathway_id), NA_character_,
                                str_squish(str_remove(kegg_pathway, "ko\\d{5}")))) |>
  select(-kegg_pathway)

# One pathway id can carry two different free-text names across rows (ko04142 is
# written both "Lysosome" and "Lysosome (saposin-like family)"). Left as-is those
# become two separate gene sets in 07 for one pathway, so a canonical label is
# assigned per id: the most frequently used name, shortest name breaking ties.
# `pathway_label` is what downstream gene sets should key on; `pathway_name`
# keeps the source text for traceability.
pathway_canonical <- pathway_parsed |>
  filter(!is.na(pathway_id)) |>
  count(pathway_id, pathway_name, name = "n_rows") |>
  arrange(pathway_id, desc(n_rows), nchar(pathway_name), pathway_name) |>
  slice_head(n = 1, by = pathway_id) |>
  transmute(pathway_id, pathway_label = paste(pathway_id, pathway_name))

het_gene_reference <- pathway_parsed |>
  left_join(pathway_canonical, by = "pathway_id") |>
  distinct(ko_id, source, category, subcategory, subcategory_raw, gene,
           pathway_id, pathway_name, pathway_label, resolution_note, source_row) |>
  arrange(ko_id, source, source_row, pathway_id)

write_csv(het_gene_reference, path_reference)

# --- report ------------------------------------------------------------------
n_ko        <- n_distinct(het_gene_reference$ko_id)
n_boundary  <- het_gene_reference |>
  filter(str_detect(coalesce(resolution_note, ""), "family_range_boundary")) |>
  distinct(ko_id) |>
  nrow()
n_pathway   <- n_distinct(na.omit(het_gene_reference$pathway_id))

message(sprintf(
  "[het-reference] %d source rows -> %d unique KOs (%d from family-range boundaries), %d KEGG pathways",
  nrow(panel_long), n_ko, n_boundary, n_pathway))
message(sprintf("[het-reference] %d source rows resolved to no KO -> %s",
                nrow(unresolved), path_unresolved))
message(sprintf("[het-reference] wrote %s (%d rows)",
                path_reference, nrow(het_gene_reference)))

het_gene_reference |>
  distinct(ko_id, category) |>
  count(category, name = "n_ko") |>
  arrange(desc(n_ko)) |>
  print(n = Inf)

# --- assertions --------------------------------------------------------------
stopifnot(
  nrow(het_gene_reference) > 0,
  all(str_detect(het_gene_reference$ko_id, "^K\\d{5}$")),   # every id is a clean KO
  !anyNA(het_gene_reference$ko_id),
  !anyNA(het_gene_reference$category),
  # no source row is silently lost: it either yields a KO or appears in the log
  setequal(panel_long$source_row,
           union(het_gene_reference$source_row, unresolved$source_row))
)

invisible(het_gene_reference)
