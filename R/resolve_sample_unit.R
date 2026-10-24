#' Resolve the physical sample unit and its biological-replicate grouping
#'
#' Establishes one row per physical metatranscriptome sample and the two ids
#' every downstream step keys on:
#'
#' * `sample_id` — unique per **physical sample** (the individual biological
#'   replicate). Duplicate sites are disambiguated with a `_r1/_r2` suffix.
#' * `sample_unit` — the **biological unit** a set of replicates belongs to
#'   (field|site, no suffix). Replicates of the same vent share one `sample_unit`.
#'
#' Decision (Hu, 2026-07-30, refined 2026-07-31): the repeated vents (Ravelin2,
#' Lotsoshrimp) and Gorda Plume are **true biological replicates**. They are
#' never silently collapsed. Instead, both ids are carried so each analysis can
#' choose, on its own terms, whether to use the replicates **individually** (the
#' statistic needs the replication / degrees of freedom — PERMANOVA, DESeq2,
#' betadisper, within-gradient RDA) or **averaged to `sample_unit`** (the
#' statistic wants one observation per biological unit — concentration/Gini,
#' cross-dataset ordination concordance). The averaging itself lives in
#' [build_matrix()] (`aggregate = "mean"`) and [collapse_replicates()]; this
#' function only supplies the grouping key and the replicate count.
#'
#' @param meta A data frame with, at minimum, `field`, `habitat`, and a
#'   `raw_label` column identifying the physical sample. Any additional columns
#'   (geochem, depth) are carried through unchanged.
#' @param field_col,habitat_col,site_col Column names for field, habitat
#'   (vent/non-vent token), and the site/vent name. Defaults match the
#'   curated `config/sample_metadata.csv` schema.
#'
#' @return `meta` with added columns: `sample_id` (unique per physical sample;
#'   duplicate sites get a `_r1/_r2` suffix), `sample_unit` (field|site grouping
#'   shared by replicates), `n_reps` (replicate count in that unit), `site`, and
#'   a `vent_bin` factor from [bin_vent()]. Row order is preserved.
#' @examples
#'   m <- data.frame(field = c("Von Damm", "Von Damm"),
#'                   habitat = c("Ravelin2", "Ravelin2"),
#'                   raw_label = c("VD_Ravelin2_A", "VD_Ravelin2_B"))
#'   resolve_sample_unit(m)[, c("sample_id", "sample_unit", "n_reps")]
#'   #> sample_id            sample_unit        n_reps
#'   #> VonDamm|Ravelin2_r1  VonDamm|Ravelin2   2
#'   #> VonDamm|Ravelin2_r2  VonDamm|Ravelin2   2
#' @seealso [bin_vent()], [build_matrix()]
#' @export
# for debugging
# meta <- readr::read_csv("config/sample_metadata.csv", show_col_types = FALSE)
# meta
resolve_sample_unit <- function(meta,
                                field_col = "field",
                                habitat_col = "habitat",
                                site_col = "site") {
  stopifnot(all(c(field_col, habitat_col) %in% names(meta)))

  field <- as.character(meta[[field_col]])
  habitat <- as.character(meta[[habitat_col]])
  site <- as.character(meta[[site_col]])

  field_slug <- gsub("\\s+", "", field)

  # biological-replicate grouping key (shared by replicates of the same vent)
  unit <- paste(field_slug, site, sep = "|")

  # replicate index within each unit so duplicates get distinct ids
  rep_index <- stats::ave(seq_along(unit), unit, FUN = seq_along)
  n_in_group <- stats::ave(seq_along(unit), unit, FUN = length)
  suffix <- ifelse(n_in_group > 1, paste0("_r", rep_index), "")

  meta$site <- site
  meta$sample_unit <- unit
  meta$n_reps <- n_in_group
  meta$sample_id <- paste0(unit, suffix)
  meta$vent_bin <- bin_vent(site)

  if (anyDuplicated(meta$sample_id)) {
    stop("sample_id not unique after resolution — check raw_label / site columns")
  }
  
  meta
}
