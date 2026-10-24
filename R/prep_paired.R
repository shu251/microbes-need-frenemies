#' Shared transcript prep for scripts 02-07
#'
#' One rule for every analysis script, so the methods describe a single data set
#' (decision 2026-09-25; curated taxonomy added 2026-10-12):
#'   * fields limited to `fields` (Axial excluded),
#'   * transcripts need a lineage resolved past "Eukaryota" and at least one KO,
#'   * NO per-row TPM cutoff. Presence/absence work (05-07) decides detection on
#'     unit totals with [unit_expression()] (>= 1 scaledTPM AND >= 10 per million
#'     of the unit's retained expression); community matrices use
#'     prevalence_filter(..., 4),
#'   * multi-KO transcripts split with TPM divided equally (mass-conserving),
#'   * Eukaryota-only and Opisthokonta lineages dropped on the FULL lineage, so the
#'     same transcripts enter every taxonomic grain,
#'   * taxonomy follows the 02 manual curation ([add_curated_taxonomy()]): every
#'     row gets `cur_supergroup`, `cur_group`, `cur_phylum`, `cur_L4` and
#'     `resolved_L4`, and `taxon` is the curated lineage at `depth`.
#'
#' @param long Output of load_long_annot().
#' @param fields Vent fields in scope.
#' @param depth Taxonomic grain for the `taxon` column ([curated_lineage()]);
#'   NULL skips it (02's ladder assigns its own depth per step).
#' @param drop_lineage_pattern Regex matched against full_classification; matching
#'   transcripts are removed.
#' @return Long table, one row per transcript x KO x library.
#' @export
prep_paired <- function(long,
                        fields = c("Gorda Ridge", "Von Damm", "Piccard"),
                        depth = 4L,
                        drop_lineage_pattern = "Opisthokonta") {
  # lineages with at least one real rank below Eukaryota (placeholders removed)
  lin_ok <- long |>
    dplyr::distinct(full_classification) |>
    dplyr::filter(!is.na(full_classification)) |>
    dplyr::mutate(.l2 = lineage_depth(full_classification, 2L)) |>
    dplyr::filter(!is.na(.l2), .l2 != "Eukaryota",
                  !grepl(drop_lineage_pattern, full_classification)) |>
    dplyr::pull(full_classification)

  out <- long |>
    dplyr::filter(field %in% fields,
                  full_classification %in% lin_ok,
                  !is.na(ko_id), !ko_id %in% c("", "-")) |>
    dplyr::mutate(.tx = dplyr::row_number()) |>
    tidyr::separate_longer_delim(ko_id, delim = ",") |>
    dplyr::mutate(ko_id = stringr::str_trim(ko_id)) |>
    dplyr::filter(!ko_id %in% c("", "-")) |>
    dplyr::mutate(scaledTPM = scaledTPM / dplyr::n(), .by = .tx) |>
    dplyr::select(-.tx) |>
    add_curated_taxonomy()

  if (!is.null(depth)) {
    key <- dplyr::distinct(out, full_classification, SUPERGROUP_0, SUPERGROUP_PHYLUM) |>
      dplyr::mutate(taxon = curated_lineage(full_classification, SUPERGROUP_0,
                                            SUPERGROUP_PHYLUM, depth))
    out <- dplyr::left_join(out, key,
                            by = c("full_classification", "SUPERGROUP_0", "SUPERGROUP_PHYLUM"))
  }

  stopifnot(all(unique(out$field) %in% fields), !anyNA(out$cur_L4))
  out
}
