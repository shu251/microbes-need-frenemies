#' Build a sample x feature matrix and its Aitchison distance
#'
#' Turns a long transcript table into an abundance matrix (summing `scaledTPM`
#' within feature), then applies the robust-CLR transform and Euclidean distance
#' = robust Aitchison distance. Using one builder for both the taxonomy and
#' function matrices guarantees Q1's two matrices are constructed identically —
#' the contrast is only meaningful if they are.
#'
#' **Replicate handling.** Rows are built at the level given by `id_col`
#' (default `sample_id`, i.e. one row per physical biological replicate). Set
#' `aggregate = "mean"` to collapse replicates to their biological unit
#' (`unit_col`, default `sample_unit`) *after* per-sample abundances are formed:
#' each unit's feature vector is the mean over its replicates, **counting
#' structural zeros** (a feature absent from one replicate contributes 0 to that
#' unit's mean, not a dropped row). The robust-CLR transform and distance are
#' then computed on whichever level was requested. Use `"none"` when the
#' statistic needs the replication (PERMANOVA, betadisper); use `"mean"` when it
#' wants one observation per unit (e.g. cross-dataset ordination concordance).
#'
#' Compositional handling follows Gloor et al. 2017 (Front. Microbiol. 8:2224):
#' metatranscriptome abundances are compositional, so a CLR-family transform is
#' not optional. `rclr` tolerates the zeros that plague sparse KO tables.
#'
#' @param long A long data frame: one row per feature x sample, with columns
#'   `id_col`, a feature column, `scaledTPM`, and (if `aggregate = "mean"`)
#'   `unit_col`.
#' @param feature_col Name of the feature column (e.g. "SUPERGROUP_PHYLUM" or "KO").
#' @param min_tpm Drop rows below this `scaledTPM` before aggregating.
#' @param transform `decostand` method; default "rclr".
#' @param id_col Physical-sample id column; default "sample_id".
#' @param aggregate One of "none" (rows = replicates) or "mean" (rows = averaged
#'   biological units). Default "none".
#' @param unit_col Biological-unit column used when `aggregate = "mean"`;
#'   default "sample_unit".
#'
#' @return A list with `mat` (RAW abundances, rows x features, rows ordered by
#'   id), `z` (the same matrix after `transform`, which `dist` is built from),
#'   `dist` (a `dist` object, robust Aitchison), and `level` ("sample" or
#'   "unit"). Align metadata to `rownames(mat)` before modelling. Use `z`, not
#'   `mat`, for anything computing per-feature statistics.
#' @examples
#'   \dontrun{
#'     tax     <- build_matrix(long_df, "SUPERGROUP_PHYLUM")                 # replicates
#'     tax_avg <- build_matrix(long_df, "SUPERGROUP_PHYLUM", aggregate = "mean") # units
#'   }
#' @seealso [run_permanova()], [resolve_sample_unit()]
#' @export
build_matrix <- function(long, feature_col, min_tpm = 1, transform = "rclr",
                         id_col = "sample_id",
                         aggregate = c("none", "mean"),
                         unit_col = "sample_unit") {
  aggregate <- match.arg(aggregate)
  requireNamespace("dplyr"); requireNamespace("tidyr"); requireNamespace("tibble")
  requireNamespace("vegan")
  stopifnot(id_col %in% names(long))
  if (aggregate == "mean" && !unit_col %in% names(long)) {
    stop("aggregate = 'mean' needs the '", unit_col,
         "' column in `long` (add it via resolve_sample_unit() / the 02 join).")
  }

  # per-physical-sample feature abundances (sum scaledTPM within sample x feature)
  per_sample <- long |>
    dplyr::filter(.data$scaledTPM >= min_tpm,
                  !is.na(.data[[feature_col]]), .data[[feature_col]] != "NA") |>
    dplyr::group_by(.id = .data[[id_col]], feature = .data[[feature_col]]) |>
    dplyr::summarise(tpm = sum(.data$scaledTPM), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "feature", values_from = "tpm",
                       values_fill = 0)

  if (aggregate == "mean") {
    # map each physical sample to its biological unit, then take the mean feature
    # vector over that unit's replicates (structural zeros already filled above,
    # so the mean correctly divides by the number of replicates in the unit).
    id2unit <- long |>
      dplyr::distinct(.id = .data[[id_col]], .unit = .data[[unit_col]])
    per_sample <- per_sample |>
      dplyr::left_join(id2unit, by = ".id") |>
      dplyr::select(-".id") |>
      dplyr::group_by(.data$.unit) |>
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean),
                       .groups = "drop") |>
      dplyr::rename(.id = ".unit")
  }

  mat <- per_sample |>
    dplyr::arrange(.data$.id) |>
    tibble::column_to_rownames(".id") |>
    as.matrix()

  z <- if (identical(transform, "rclr")) rclr_transform(mat)
       else vegan::decostand(mat, method = transform)
  d <- vegan::vegdist(z, method = "euclidean")
  # `mat` is RAW abundance; `z` is the transformed matrix `dist` was built from.
  # Per-feature statistics want `z` (see 07), never `mat`.
  list(mat = mat, z = z, dist = d,
       level = if (aggregate == "mean") "unit" else "sample")
}


#' Robust CLR with zeros kept at zero (no imputation)
#'
#' For each sample (row): log of every non-zero value minus the mean log of
#' that row's non-zero values; zeros stay 0. Euclidean distance on the result
#' is the robust Aitchison distance (Gloor et al. 2017, Front. Microbiol.
#' 8:2224), and it keeps every dimension of the data.
#'
#' **Why this is written out instead of calling `vegan::decostand("rclr")`.**
#' vegan <= 2.6-4 computed exactly this. Newer vegan versions changed the
#' default to `impute = TRUE`, which fills the zeros by OptSpace matrix
#' completion at rank `ropt = 3` and returns the rank-3 reconstruction of the
#' WHOLE matrix (the RPCA approach of Martino et al. 2019, mSystems
#' 4:e00016-19). That silently turned every distance into a 3-axis
#' approximation (found 2026-09-23: Procrustes variance captured = 1.0 at
#' k = 3). Writing the transform here makes the result independent of the
#' installed vegan version. Same convention as the per-unit rclr in
#' `05_trophic-enrichment.qmd`.
#'
#' @param mat Numeric matrix, samples in rows, non-negative abundances.
#' @return Matrix of the same shape; rows sum to 0 over their non-zero entries.
#' @export
rclr_transform <- function(mat) {
  mat <- as.matrix(mat)
  if (any(mat < 0, na.rm = TRUE)) stop("rclr needs non-negative data.", call. = FALSE)
  out <- t(apply(mat, 1, function(row) {
    nz  <- row > 0
    res <- numeric(length(row))
    res[nz] <- log(row[nz]) - mean(log(row[nz]))
    res
  }))
  dimnames(out) <- dimnames(mat)
  out
}
