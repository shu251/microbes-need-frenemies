#' Keep only features present in at least `min_samples` samples
#'
#' A prevalence floor for compositional matrices. At fine taxonomic (or
#' functional) grain most features are singletons — present in a single library —
#' and singleton-dominated robust-Aitchison distances inflate among-sample
#' separation, and with it PERMANOVA R2 (Gloor et al. 2017, Front. Microbiol.
#' 8:2224). Filtering to features seen in at least `min_samples` libraries tests
#' fine grains on shared signal rather than noise. Apply BEFORE [build_matrix()],
#' on the same feature column and id you will build on.
#'
#' @param long Long transcript table.
#' @param feature_col Feature column name (e.g. a truncated lineage or "ko_id").
#' @param min_samples Minimum distinct samples a feature must appear in (default 2).
#' @param id_col Sample-id column (default "sample_id").
#' @param min_tpm Count a feature as present in a sample only above this TPM
#'   (default 0, i.e. any non-zero expression).
#' @return `long` filtered to features meeting the prevalence floor.
#' @seealso [build_matrix()], [lineage_depth()]
#' @export
prevalence_filter <- function(long, feature_col, min_samples = 2,
                              id_col = "sample_id", min_tpm = 0) {
  keep <- long |>
    dplyr::filter(.data$scaledTPM > min_tpm, !is.na(.data[[feature_col]])) |>
    dplyr::distinct(.feat = .data[[feature_col]], .id = .data[[id_col]]) |>
    dplyr::count(.feat, name = "n_samp") |>
    dplyr::filter(.data$n_samp >= min_samples) |>
    dplyr::pull(.feat)
  long |> dplyr::filter(.data[[feature_col]] %in% keep)
}
