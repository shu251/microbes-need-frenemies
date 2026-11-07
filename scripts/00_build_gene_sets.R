#!/usr/bin/env Rscript
# 00_build_gene_sets.R
# Builds every KO gene list the analysis reads, in one pass. Replaces
# 00_build_het_gene_reference.R and 00b_build_trophic_gene_sets.R (merged
# 2026-10-25; originals in scripts/_to_delete/pre_gene_set_merge_20261025/).
#
# RUN (from the project root, or source() from anywhere in the project)
#   Rscript scripts/00_build_gene_sets.R
#
# PART 1  Heterotrophy / phagotrophy reference panel
#   IN   data/gene-lists/marker_panel_expansion_candidates.csv   (hand-curated)
#        data/gene-lists/phagotrophy_genes_protists_1.csv        (hand-curated)
#   OUT  data/gene-lists/het_gene_reference.csv                  (read by 06, 07)
#        output-tables/00_het_gene_reference_unresolved.csv
#
# PART 2  Trophic gene sets (read by 06 and 07)
#     kegg_level_c          KEGG pathway maps (level C) inside four level-B categories:
#                           Transport and catabolism, Signal transduction,
#                           Cell growth and death, Cell motility. Only pathways kept
#                           in data/gene-lists/kegg_levelC_curation.csv.
#     burns2018             Protein families that predict phagocytosis, photosynthesis
#                           and prototrophy across 35 eukaryotic genomes
#                           (Burns et al. 2018 Nat Ecol Evol 2:697-704), mapped from
#                           their UniProt representatives to KOs.
#     protist_experimental  KOs for genes whose loss impairs phagotrophy in a protist
#                           (knockout, knockdown or silencing), one citation per row.
#   IN   data/external/ko00001.keg
#        data/gene-lists/kegg_levelC_curation.csv                (hand-curated)
#        data/external/burns2018/modelGenes{PH,PHOTO,Prototrophy}.txt
#        data/external/burns2018/burns2018_uniprot_ko.csv        (cache; built on first run)
#        data/gene-lists/protist_phagotrophy_experimental.csv    (hand-curated)
#   OUT  data/gene-lists/trophic_gene_sets.csv
#        output-tables/00b_het_reference_ko_names.csv            (curated label vs KEGG name)
#
# NETWORK
#   Only the first run needs the internet: UniProt accession -> KEGG gene -> KO via the
#   KEGG REST API. Cached in data/external/burns2018/burns2018_uniprot_ko.csv.

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
if (!identical(normalizePath(getwd(), winslash = "/"),
               normalizePath(project_root, winslash = "/"))) {
  setwd(project_root)
  message("[gene-sets] working directory set to project root: ", project_root)
}

source("R/read_brite.R")


# ==============================================================================
# PART 1. Heterotrophy / phagotrophy reference panel
# ==============================================================================

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


# ==============================================================================
# PART 2. Trophic gene sets
# ==============================================================================

# --- paths and constants (part 2) ----------------------------------------------
path_keg          <- "data/external/ko00001.keg"
dir_burns         <- "data/external/burns2018"
path_burns_cache  <- file.path(dir_burns, "burns2018_uniprot_ko.csv")
path_experimental <- "data/gene-lists/protist_phagotrophy_experimental.csv"
path_sets_out     <- "data/gene-lists/trophic_gene_sets.csv"
path_names_out    <- "output-tables/00b_het_reference_ko_names.csv"

# Level-B categories asked about (2026-10-21)
LEVEL_B_KEEP <- c("09141" = "Transport and catabolism",
                  "09132" = "Signal transduction",
                  "09143" = "Cell growth and death",
                  "09142" = "Cell motility")

BURNS_FILES <- c(phagocytosis   = "modelGenesPH.txt",
                 photosynthesis = "modelGenesPHOTO.txt",
                 prototrophy    = "modelGenesPrototrophy.txt")

# --- 1. KEGG level-C pathway sets ----------------------------------------------
# Pathways come from the protist curation sheet (read_kegg_curation(), 2026-10-25):
# prokaryote-, plant-, animal- and virus-specific maps are not gene sets.
brite_c <- read_brite_level_c(path_keg, keep_b = names(LEVEL_B_KEEP))

kegg_c_sets <- brite_c |>
  transmute(set_source = "kegg_level_c",
            set_name   = paste(pathway_id, pathway_name),
            category   = b_description,
            ko_id,
            evidence   = "KEGG pathway membership",
            citation   = "KEGG BRITE ko00001")

message(sprintf("[trophic-sets] KEGG level C (curated): %d pathway sets, %d KOs",
                n_distinct(kegg_c_sets$set_name), n_distinct(kegg_c_sets$ko_id)))

# --- 2. Burns et al. 2018 sets: UniProt representative -> KO ---------------------
burns_acc <- imap_dfr(BURNS_FILES, \(file, model) {
  tibble(model = model,
         uniprot = read_lines(file.path(dir_burns, file)) |> str_squish()) |>
    filter(uniprot != "")
})

# KEGG REST returns tab-separated pairs; one request per chunk of ids joined by "+".
kegg_rest_pairs <- function(operation, ids, chunk_size = 50) {
  empty <- tibble(from = character(), to = character())
  if (!length(ids)) return(empty)
  pairs <- split(ids, ceiling(seq_along(ids) / chunk_size)) |>
    map(\(chunk) {
      url <- paste0("https://rest.kegg.jp/", operation, "/", paste(chunk, collapse = "+"))
      Sys.sleep(0.4)                                   # stay well under KEGG's rate limit
      txt <- tryCatch(read_lines(url), error = \(e) character(), warning = \(w) character())
      txt <- txt[str_detect(txt, "\t")]
      if (!length(txt)) return(empty)
      tibble(from = str_extract(txt, "^[^\t]+"), to = str_extract(txt, "(?<=\t).+$"))
    }) |>
    list_rbind()
  if (!nrow(pairs)) empty else pairs
}

if (!file.exists(path_burns_cache)) {
  message("[trophic-sets] no cache: mapping ", n_distinct(burns_acc$uniprot),
          " UniProt accessions to KOs through the KEGG REST API (first run only)")
  up_to_gene <- kegg_rest_pairs("conv/genes", paste0("up:", unique(burns_acc$uniprot))) |>
    transmute(uniprot = str_remove(from, "^up:"), kegg_gene = to)
  gene_to_ko <- kegg_rest_pairs("link/ko", unique(up_to_gene$kegg_gene)) |>
    transmute(kegg_gene = from, ko_id = str_remove(to, "^ko:"))
  burns_map <- up_to_gene |>
    left_join(gene_to_ko, by = "kegg_gene", relationship = "many-to-many") |>
    distinct(uniprot, kegg_gene, ko_id)
  if (nrow(burns_map)) {
    write_csv(burns_map, path_burns_cache)
  } else {
    warning("KEGG REST returned nothing (offline?). Burns 2018 sets are left OUT of this ",
            "build and nothing was cached; re-run with a network connection.", call. = FALSE)
  }
}
burns_map <- if (file.exists(path_burns_cache)) {
  read_csv(path_burns_cache, col_types = cols(.default = "c"))
} else {
  tibble(uniprot = character(), kegg_gene = character(), ko_id = character())
}

burns_sets <- burns_acc |>
  left_join(burns_map |> filter(!is.na(ko_id)) |> distinct(uniprot, ko_id),
            by = "uniprot", relationship = "many-to-many")

burns_coverage <- burns_sets |>
  summarise(n_families = n_distinct(uniprot),
            n_with_ko  = n_distinct(uniprot[!is.na(ko_id)]),
            n_kos      = n_distinct(ko_id, na.rm = TRUE), .by = model)
message("[trophic-sets] Burns 2018 families mapped to a KO:")
print(burns_coverage)

burns_sets <- burns_sets |>
  filter(!is.na(ko_id)) |>
  distinct(model, ko_id) |>
  transmute(set_source = "burns2018",
            set_name   = paste("Burns2018", model),
            category   = model,
            ko_id,
            evidence   = "protein family predictive of trophic mode across 35 eukaryotic genomes",
            citation   = "Burns et al. 2018 Nat Ecol Evol 2:697-704 doi:10.1038/s41559-018-0477-7")

# --- 3. Protist experimental phagotrophy set -----------------------------------
experimental <- read_csv(path_experimental, show_col_types = FALSE)

experimental_sets <- experimental |>
  filter(!is.na(ko_id)) |>
  separate_longer_delim(ko_id, ";") |>
  mutate(ko_id = str_squish(ko_id)) |>
  transmute(set_source = "protist_experimental",
            set_name   = "Protist phagotrophy (experimental)",
            category   = gene_label,
            ko_id,
            evidence   = paste0(organism, ": ", evidence),
            citation   = paste0(citation, " doi:", doi))

no_ko <- experimental |> filter(is.na(ko_id))
if (nrow(no_ko)) {
  message("[trophic-sets] experimental genes with no KEGG KO (kept in the CSV, not in any set): ",
          paste(no_ko$gene_label, collapse = "; "))
}

# --- combine, check against ko00001, write -------------------------------------
ko_names <- tibble(line = read_lines(path_keg)) |>
  filter(str_detect(line, "^D")) |>
  transmute(ko_id   = str_extract(str_squish(str_sub(line, 2)), "^K\\d{5}"),
            ko_name = str_squish(str_remove(str_squish(str_sub(line, 2)), "^K\\d{5}"))) |>
  filter(!is.na(ko_id)) |>
  distinct(ko_id, .keep_all = TRUE)

trophic_sets <- bind_rows(kegg_c_sets, burns_sets, experimental_sets) |>
  distinct(set_source, set_name, ko_id, .keep_all = TRUE) |>
  left_join(ko_names, by = "ko_id")

missing_names <- trophic_sets |> filter(is.na(ko_name)) |> distinct(set_source, ko_id)
if (nrow(missing_names)) {
  message("[trophic-sets] ", nrow(missing_names), " KO(s) absent from this ko00001.keg ",
          "(retired or newer than the file):")
  print(missing_names)
}

stopifnot(
  !anyNA(trophic_sets$ko_id),
  all(str_detect(trophic_sets$ko_id, "^K\\d{5}$")),
  !anyDuplicated(trophic_sets[c("set_source", "set_name", "ko_id")])
)

write_csv(trophic_sets, path_sets_out)
message("[trophic-sets] wrote ", path_sets_out, ": ")
print(trophic_sets |> summarise(n_sets = n_distinct(set_name), n_kos = n_distinct(ko_id),
                                .by = set_source))

# --- 4. het_gene_reference: curated label next to KEGG's name ------------------
# A KO typed one digit off, or copied from a neighbouring row, silently puts a
# different protein in the panel. This table lines up the label each row was curated
# under with KEGG's own name for that KO so the panel can be checked by eye.
# `label_in_kegg_name` is a crude screen (does any 3+ character token of the label
# appear in the KEGG name?) -- FALSE means "read this row", TRUE does not prove a match.
het_names <- het_gene_reference |>
  distinct(ko_id, gene, category, subcategory, source) |>
  left_join(ko_names, by = "ko_id") |>
  mutate(label_in_kegg_name = map2_lgl(gene, ko_name, \(g, k) {
    if (is.na(k)) return(FALSE)
    toks <- str_split_1(g, "[^A-Za-z0-9]+")
    toks <- toks[nchar(toks) >= 3]
    any(str_detect(str_to_lower(k), fixed(str_to_lower(toks))))
  })) |>
  arrange(label_in_kegg_name, ko_id)

write_csv(het_names, path_names_out)
message("[trophic-sets] wrote ", path_names_out, ": ",
        sum(!het_names$label_in_kegg_name), " of ", nrow(het_names),
        " label/KO pairs flagged for a read")

invisible(list(het_gene_reference = het_gene_reference, trophic_sets = trophic_sets))
