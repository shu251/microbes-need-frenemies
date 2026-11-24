# validate_pfam_ids.R
# Check every Pfam accession in one or more gene-list CSVs against InterPro.
# For each accession: does it exist, what are its official name and short name,
# and does the Pfam_Domain_Name written in the CSV match either one?
#
# Run from the project root:
#   source("scripts/validate_pfam_ids.R")
# Needs internet access to www.ebi.ac.uk. Writes output-tables/qc-pfam-id-validation.csv.

library(tidyverse)
library(httr2)

csv_paths <- c(
  "data/gene-lists/phagotrophy_genes_protists_1.csv",
  "data/gene-lists/het_panel_new_rows_v3_20261107.csv",
  "data/gene-lists/symbiosis_parasitism_genes_protists_20261108.csv"   # >SKH: 2026-11-08
)
name_map_path <- "data/gene-lists/pfam_name_map.csv"      # eggNOG older names (e.g. DUF839)
out_path      <- "output-tables/qc-pfam-id-validation.csv"

# ---- 1. collect accessions and the names written beside them -----------------
read_panel <- function(path) {
  read_csv(path, col_types = cols(.default = "c"), show_col_types = FALSE) |>
    mutate(source_file = basename(path))
}

missing_paths <- csv_paths[!file.exists(csv_paths)]
if (length(missing_paths) > 0) {
  stop("Gene-list file(s) not found: ", paste(missing_paths, collapse = ", "))
}

claimed <- csv_paths |>
  map_dfr(read_panel) |>
  select(source_file, Gene_Symbol, Pfam_ID, Pfam_Domain_Name) |>
  filter(!is.na(Pfam_ID)) |>
  mutate(pfam_acc = str_extract_all(Pfam_ID, "PF\\d{5}")) |>
  unnest(pfam_acc)

# ---- 2. query InterPro (cached, polite pacing) ----------------------------------
query_pfam <- function(acc) {
  Sys.sleep(1)
  resp <- request("https://www.ebi.ac.uk/interpro/api/entry/pfam") |>
    req_url_path_append(acc) |>
    req_retry(max_tries = 4, backoff = ~ 10) |>
    req_error(is_error = \(r) FALSE) |>
    req_perform()

  if (resp_status(resp) != 200) {
    return(tibble(pfam_acc = acc, found = FALSE,
                  interpro_name = NA_character_, interpro_short = NA_character_,
                  http_status = resp_status(resp)))
  }
  meta <- resp_body_json(resp)$metadata
  tibble(pfam_acc = acc, found = TRUE,
         interpro_name  = meta$name$name  %||% NA_character_,
         interpro_short = meta$name$short %||% NA_character_,
         http_status = 200L)
}

interpro <- claimed |>
  distinct(pfam_acc) |>
  pull(pfam_acc) |>
  map_dfr(query_pfam)

# ---- 3. compare written names with InterPro ------------------------------------
older_names <- if (file.exists(name_map_path)) {
  read_csv(name_map_path, show_col_types = FALSE) |> select(pfam_acc, eggnog_name)
} else {
  tibble(pfam_acc = character(), eggnog_name = character())
}

norm <- \(x) str_to_lower(str_replace_all(x, "[^A-Za-z0-9]", ""))

result <- claimed |>
  left_join(interpro,   by = "pfam_acc") |>
  left_join(older_names, by = "pfam_acc") |>
  mutate(
    name_match = pmap_lgl(
      list(Pfam_Domain_Name, interpro_short, interpro_name, eggnog_name),
      \(written, short, long, old) {
        w <- norm(str_split(written, "[;,]|\\s{2,}")[[1]])
        w <- w[nzchar(w)]
        any(w %in% norm(na.omit(c(short, long, old))))
      }
    ),
    status = case_when(
      !found      ~ "ACCESSION NOT FOUND",
      name_match  ~ "ok",
      TRUE        ~ "exists, name differs: review"
    )
  ) |>
  arrange(status == "ok", pfam_acc)

write_csv(result, out_path)
count(result, status) |> print()
result |> filter(status != "ok") |> print(n = Inf)
