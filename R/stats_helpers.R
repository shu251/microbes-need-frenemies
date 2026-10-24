# Statistical functions R scripts
## Thin, documented wrappers so every question uses the same conventions
## (marginal PERMANOVA terms, seeded permutations, Aitchison distances) and the
## rigor fixes from the 2026-07-30 review are baked in rather than re-typed.

#' @name stats_helpers
NULL

## Marginal PERMANOVA on one distance matrix

#' @param d A `dist` object (robust Aitchison).
#' @param meta Sample metadata, rows aligned to the labels of `d`.
#' @param formula_rhs Right-hand side, default `~ vent_bin + field`.
#' @param by `adonis2` term test; default "margin" (order-independent).
#' @param permutations,seed Permutation count and RNG seed.
#' @return A tidy data frame of terms with `Df`, `R2`, `F`, `p`.
#' @export
run_permanova <- function(d, meta, 
                          formula_rhs = ~ vent_bin + field,
                          by = "margin", 
                          permutations = 9999,
                          seed = 20260730) {
  requireNamespace("vegan"); requireNamespace("tibble") # require
  set.seed(seed)
  f <- stats::as.formula(paste("d", paste(deparse(formula_rhs), collapse = "")))
  # adonis 2 runs the permutational mutlivariate analysis of variance from a distace matrix
  a <- vegan::adonis2(f, data = meta, by = by, permutations = permutations)
  as_tibble(as.data.frame(a), rownames = "term") |>
    rename(p = dplyr::any_of("Pr(>F)"))
}

#' Dispersion check (betadisper) — report alongside every PERMANOVA
#'
#' PERMANOVA confounds location with dispersion, and that confound is worst
#' under the unbalanced design here (Anderson & Walsh 2013, Ecol. Monogr.
#' 83:557). A significant term with heterogeneous dispersion is ambiguous —
#' so this is a main-text result, not a footnote.
#'
#' @param d A `dist` object.
#' @param group Grouping factor (e.g. `meta$field`).
#' @return One-row tibble: `dispersion_p` and an interpretation `note`.
#' @export
betadisper_check <- function(d, group, permutations = 9999, seed = 20260730) {
  requireNamespace("vegan"); requireNamespace("tibble")
  set.seed(seed)
  bd <- vegan::betadisper(d, group)
  p <- vegan::permutest(bd, permutations = permutations)$tab$`Pr(>F)`[1]
  tibble::tibble(
    dispersion_p = p,
    note = if (isTRUE(p < 0.05))
      "groups differ in dispersion — interpret PERMANOVA location effect with caution"
    else
      "dispersion homogeneous — PERMANOVA location effect is clean"
  )
}

#' Procrustes / co-structure test between two ordinations (Q1, Q4)
#'
#' The honest test of whether two matrices share structure. Preferred over the
#' Mantel test, which has low power and inflated type-I error for comparing
#' distance matrices (Legendre, Fortin & Borcard 2015, MEE 6:1239;
#' Peres-Neto & Jackson 2001, Oecologia 129:169).
#'
#' @param d1,d2 Two `dist` objects on the SAME samples (same label order).
#' @param k Number of PCoA axes to superimpose.
#' @return The `protest` object (report `$t0` correlation and `$signif`).
#' @export
procrustes_test <- function(d1, d2, k = 2, permutations = 9999, seed = 20260730) {
  requireNamespace("vegan")
  set.seed(seed)
  o1 <- stats::cmdscale(d1, k = k)
  o2 <- stats::cmdscale(d2, k = k)
  vegan::protest(o1, o2, permutations = permutations)
}

#' Prune collinear predictors by variance inflation factor (Q5 RDA)
#'
#' Drops predictors above `thresh` VIF one at a time. Use before forward
#' selection so the vent-only RDA (~13 samples, ~7 geochem predictors) is not
#' overfit or dominated by collinear terms (Blanchet, Legendre & Borcard 2008,
#' Ecology 89:2623).
#'
#' @param env Data frame of numeric environmental predictors.
#' @param thresh VIF cutoff (default 10; 5 is stricter).
#' @return Character vector of retained predictor names.
#' @export
vif_select <- function(env, thresh = 10) {
  env <- env[stats::complete.cases(env), , drop = FALSE]
  keep <- names(env)
  repeat {
    if (length(keep) < 2) break
    vifs <- vapply(keep, function(v) {
      f <- stats::as.formula(paste(v, "~", paste(setdiff(keep, v), collapse = "+")))
      r2 <- summary(stats::lm(f, data = env))$r.squared
      1 / (1 - r2)
    }, numeric(1))
    if (max(vifs) < thresh) break
    keep <- setdiff(keep, names(which.max(vifs)))
  }
  keep
}
