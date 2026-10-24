#' Cliff's delta: a sign-safe nonparametric effect size for two-group contrasts
#'
#' Ranks features (KOs) by how different they are between two groups of sample
#' units, for use as the `stats` argument to `fgsea::fgsea()`.
#'
#' delta = (# pairs where focal > other - # pairs where focal < other) /
#' (n_focal * n_other), bounded [-1, 1]. Positive means the focal group ranks
#' higher; the magnitude is the net share of cross-group unit pairs that agree
#' on the direction (Cliff 1993, Psychol Bull 114:494-509; Vargha & Delaney
#' 2000, J Educ Behav Stat 25:101-132).
#'
#' Chosen over a Wilcoxon-derived rank-biserial r, to which it is algebraically
#' equivalent, because `wilcox.test()`'s internal U convention makes the sign of
#' the result depend on factor level order - a silent way to invert every
#' enrichment direction.
#'
#' Why an effect size and not a per-feature test: the KO-level permutation
#' differential was dropped from the paper (2026-09-25) because at 13 vent and 6
#' non-vent units no permutation p survives BH. fgsea never tests a KO on its
#' own - it asks whether a gene set sits toward one end of a ranked list - so it
#' needs a ranking that is directionally honest, not one filtered to
#' individually significant features (Subramanian et al. 2005, PNAS
#' 102:15545-15550).
#'
#' @param x_focal,x_other Numeric vectors, e.g. one KO's rclr values in each
#'   group. Non-finite values are dropped before comparison.
#' @return Single numeric delta, or NA when either side is empty after dropping
#'   non-finite values.
#' @examples
#'   cliffs_delta(c(3, 4, 5), c(1, 2))    #> 1
#'   cliffs_delta(c(1, 2), c(3, 4, 5))    #> -1
#' @seealso [contrast_ranks()], [rclr_transform()]
#' @export
cliffs_delta <- function(x_focal, x_other) {
  x_focal <- x_focal[is.finite(x_focal)]
  x_other <- x_other[is.finite(x_other)]
  if (length(x_focal) == 0 || length(x_other) == 0) return(NA_real_)
  mean(outer(x_focal, x_other, function(a, b) sign(a - b)))
}


#' Per-feature Cliff's delta, focal group versus everything else
#'
#' Applies [cliffs_delta()] column-wise to a units x features matrix, splitting
#' rows into `focal` versus all other rows. Returns one signed value per
#' feature, ready to hand to `fgsea(stats = ...)` after dropping non-finite
#' entries.
#'
#' Works for any one-versus-rest split of the units: vent vs. non-vent, or one
#' vent field vs. the remaining fields. Levels other than `focal` are pooled
#' into the comparison side, so a three-field design needs one call per field
#' rather than a reshaped matrix.
#'
#' @param z Numeric matrix, sample units in (named) rows, features in (named)
#'   columns. Pass transformed values - `rclr_transform(build_matrix(...)$mat)`,
#'   not `$mat`, which is raw abundance.
#' @param group Named vector or factor of group labels, named by the same unit
#'   ids as `rownames(z)`.
#' @param focal The group level treated as positive.
#' @return Named numeric vector over `colnames(z)`; NA where a feature has no
#'   finite value on one side of the split.
#' @seealso [cliffs_delta()], [build_matrix()], [rclr_transform()]
#' @export
contrast_ranks <- function(z, group, focal) {
  stopifnot(is.matrix(z), !is.null(rownames(z)), !is.null(colnames(z)))
  missing_units <- setdiff(rownames(z), names(group))
  if (length(missing_units) > 0) {
    stop("`group` has no label for ", length(missing_units), " unit(s) in `z`: ",
         paste(utils::head(missing_units, 3), collapse = ", "), call. = FALSE)
  }
  is_focal <- group[rownames(z)] == focal
  if (!any(is_focal) || all(is_focal)) {
    stop("`focal` must split `z`'s rows into two non-empty sides (focal = \"",
         focal, "\").", call. = FALSE)
  }
  purrr::set_names(
    purrr::map_dbl(seq_len(ncol(z)),
                   function(j) cliffs_delta(z[is_focal, j], z[!is_focal, j])),
    colnames(z)
  )
}
