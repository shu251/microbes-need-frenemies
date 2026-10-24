#' Load the handoff table and attach replicate-aware sample ids
#'
#' The single data-load every Q file calls. Loads 02's saved RData
#' (`long_frenemies_annot_tax`) and derives the ids downstream code keys on from
#' the long table's OWN physical-library key (`LIBRARY_NUM` + `SAMPLEID`), NOT a
#' `sample_meta` join on `raw_label`.
#'
#' Why derive from the long table itself: three vents are TRUE biological
#' replicates that SHARE a `raw_label` — `Von Damm_Vent_Ravelin2`,
#' `Piccard_Vent_Lotsoshrimp`, and `Gorda Ridge_Non-vent_Plume` (two libraries
#' each; confirmed in config/sample_metadata.csv). A `left_join(sample_meta,
#' by = "raw_label")` would silently DUPLICATE every transcript in those units.
#' `LIBRARY_NUM` + `SAMPLEID` (carried down from 01's sample-name parse) separate
#' the libraries cleanly, so the replicate structure is recoverable from the data
#' itself — no fragile row-order match to the metadata csv.
#'
#' Multi-KO transcripts (comma-joined `ko_id`) are deliberately NOT split here:
#' that is a per-question choice (Q1 explodes with mass-conserving allocation;
#' Q2 needs integer counts), so it lives in the question file, not this shared
#' loader.
#'
#' @param path RData path; default the 02 handoff object.
#' @param long Optional, an already-loaded `long_frenemies_annot_tax`. When
#'   given, `path` is not read (saves a second ~420 MB load in 03, which also
#'   needs the other objects in the same RData).
#' @return `long_frenemies_annot_tax` with added columns:
#'   * `sample_id`   — unique per physical library (biological replicate).
#'   * `sample_unit` — the biological unit (= `raw_label`), shared by replicates.
#'   * `n_reps`      — libraries in that unit.
#'   * `field`, `vent_bin` — design factors, copied from VENT_FIELD / VENT_BIN.
#' @seealso [build_matrix()], [resolve_sample_unit()]
#' @export
load_long_annot <- function(path = "data/processed/long_frenemies_annot_tax.RData",
                            long = NULL) {
  if (is.null(long)) {
    if (!file.exists(path)) {
      stop("Handoff RData not found at '", path, "'.\n",
           "Run scripts/02_taxonomy_grain_ladder.qmd with WRITE_HANDOFF <- TRUE ",
           "first.", call. = FALSE)
    }
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    if (!exists("long_frenemies_annot_tax", envir = e)) {
      stop("Object 'long_frenemies_annot_tax' not found in ", path, call. = FALSE)
    }
    long <- e$long_frenemies_annot_tax
  }

  need <- c("raw_label", "VENT_FIELD", "VENT_BIN",
            "LIBRARY_NUM", "SAMPLEID", "ko_id", "scaledTPM")
  miss <- setdiff(need, names(long))
  if (length(miss)) {
    stop("long_frenemies_annot_tax is missing column(s): ",
         paste(miss, collapse = ", "), ".\n",
         "LIBRARY_NUM / SAMPLEID come from 01's sample-name parse; ",
         "ko_id / raw_label from 02. Check those steps ran.", call. = FALSE)
  }

  long |>
    dplyr::mutate(
      sample_unit = raw_label,
      sample_id   = paste(raw_label, LIBRARY_NUM, SAMPLEID,
                          sep = "|"),
      field       = VENT_FIELD,
      vent_bin    = VENT_BIN
    ) |>
    dplyr::group_by(sample_unit) |>
    dplyr::mutate(n_reps = dplyr::n_distinct(.data$sample_id)) |>
    dplyr::ungroup()
}
