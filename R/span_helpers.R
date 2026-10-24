# Helpers for 05-07: detection, habitat occupancy, and where / among which taxa
# a KO is expressed. Decisions 2026-10-12:
#   * detection  = unit total >= 1 scaledTPM AND >= 10 per million of the unit's
#                  retained expression (params$detect_min_tpm, detect_min_per_million).
#                  scaledTPM carries library size (Soneson, Love & Robinson 2015,
#                  F1000Research 4:1521), and retained protist expression differs
#                  ~170-fold among units, so a TPM floor alone is depth-dependent.
#   * occupancy  = fraction of units in a habitat where a feature is detected;
#                  associated = occupancy >= params$occupancy_min_frac (0.30:
#                  >= 4 of 13 vent units, >= 2 of 6 non-vent units).
#   * vent-only  = vent-associated AND detected in 0 non-vent units (strict).
#   * carriers   = curated grain-4 taxa whose taxon x KO pair is detected in
#                  >= PAIR_MIN_UNITS vent units.
#   * span       = conservative lowest common ancestor of the carriers.


#' Per-unit expression of one or more feature columns, with detection
#'
#' Sums scaledTPM within unit x feature, divides by the unit's replicate count
#' (zero-aware mean over libraries), and expresses it as a share of the unit's
#' total retained expression. The unit total is taken over ALL rows of `df`, so
#' pass the full prepped table, not a feature-filtered subset.
#'
#' @param df Prepped long table (prep_paired()) with `scaledTPM`, `n_reps`,
#'   `vent_bin`, `field` and the unit column.
#' @param by Character vector of feature columns (e.g. "ko_id", or
#'   c("cur_L4", "ko_id") for taxon x KO pairs).
#' @param min_tpm,min_per_million Detection floors.
#' @param unit_col Biological-unit column.
#' @return Tibble: unit, vent_bin, field, the `by` columns, tpm, unit_total,
#'   share, per_million, detected. Only non-zero unit x feature rows; use
#'   [complete_units()] before averaging over units.
#' @export
unit_expression <- function(df, by, min_tpm = 1, min_per_million = 10,
                            unit_col = "sample_unit") {
  stopifnot(all(by %in% names(df)), unit_col %in% names(df))
  df <- dplyr::mutate(df, unit = .data[[unit_col]])
  totals <- df |>
    dplyr::summarise(unit_total = sum(scaledTPM) / dplyr::first(n_reps),
                     .by = c(unit, vent_bin, field))
  df |>
    dplyr::summarise(tpm = sum(scaledTPM) / dplyr::first(n_reps),
                     .by = dplyr::all_of(c("unit", by))) |>
    dplyr::left_join(totals, by = "unit") |>
    dplyr::mutate(share       = tpm / unit_total,
                  per_million = share * 1e6,
                  detected    = tpm >= min_tpm & per_million >= min_per_million) |>
    dplyr::relocate(unit, vent_bin, field)
}


#' Zero-fill a unit x feature table over every unit
#'
#' Means and medians across units must count units where a feature is absent.
#' @param ue Output of [unit_expression()].
#' @param units Tibble of all units with `unit`, `vent_bin`, `field`.
#' @param by Feature columns.
#' @export
complete_units <- function(ue, units, by) {
  ue |>
    dplyr::select(dplyr::all_of(c("unit", by, "tpm", "share", "detected"))) |>
    tidyr::complete(unit = units$unit, tidyr::nesting(!!!rlang::syms(by)),
                    fill = list(tpm = 0, share = 0, detected = FALSE)) |>
    dplyr::left_join(units, by = "unit")
}


#' Occupancy and share of each feature in each habitat, plus its habitat class
#'
#' @param ue Output of [unit_expression()].
#' @param units Tibble of all units (`unit`, `vent_bin`, `field`).
#' @param by Feature column(s).
#' @param min_frac Occupancy needed to call a feature associated with a habitat.
#' @return One row per feature: n_det / occ / mean_share / median_share /
#'   max_share for vent and non-vent, and `habitat_class`:
#'   vent_only (vent-associated, detected in no non-vent unit),
#'   vent_skewed (vent-associated, detected in some non-vent units but below
#'   `min_frac`), shared (associated with both), nonvent_skewed, nonvent_only,
#'   and below_filter (associated with neither habitat).
#' @export
habitat_occupancy <- function(ue, units, by, min_frac = 0.30) {
  n_units <- dplyr::count(units, vent_bin, name = "n_units")
  complete_units(ue, units, by) |>
    dplyr::summarise(n_det        = sum(detected),
                     mean_share   = mean(share),
                     median_share = stats::median(share),
                     max_share    = max(share),
                     .by = dplyr::all_of(c(by, "vent_bin"))) |>
    dplyr::left_join(n_units, by = "vent_bin") |>
    dplyr::mutate(occ = n_det / n_units,
                  hab = dplyr::if_else(vent_bin == "Vent", "vent", "nonvent")) |>
    dplyr::select(-vent_bin, -n_units) |>
    tidyr::pivot_wider(names_from = hab,
                       values_from = c(n_det, occ, mean_share, median_share, max_share)) |>
    dplyr::mutate(habitat_class = classify_habitat(occ_vent, occ_nonvent, n_det_vent,
                                                   n_det_nonvent, min_frac))
}


#' Habitat class from occupancy (see [habitat_occupancy()])
#' @export
classify_habitat <- function(occ_vent, occ_nonvent, n_det_vent, n_det_nonvent,
                             min_frac = 0.30) {
  v <- occ_vent >= min_frac
  n <- occ_nonvent >= min_frac
  factor(dplyr::case_when(
    v & n                       ~ "shared",
    v & n_det_nonvent == 0      ~ "vent_only",
    v                           ~ "vent_skewed",
    n & n_det_vent == 0         ~ "nonvent_only",
    n                           ~ "nonvent_skewed",
    TRUE                        ~ "below_filter"),
    levels = c("vent_only", "vent_skewed", "shared", "nonvent_skewed",
               "nonvent_only", "below_filter"))
}


#' Field span of each feature within vents
#'
#' @param ue Output of [unit_expression()] restricted to vent units.
#' @param units Vent units (`unit`, `field`).
#' @param by Feature column.
#' @return One row per feature: n_vent_units, field_span (number of fields with
#'   >= 1 detected unit), and occ_<field> (fraction of that field's units).
#' @export
field_span <- function(ue, units, by = "ko_id") {
  per_field <- complete_units(ue, units, by) |>
    dplyr::summarise(n_det = sum(detected), n = dplyr::n(),
                     .by = dplyr::all_of(c(by, "field")))
  per_field |>
    dplyr::summarise(n_vent_units = sum(n_det), field_span = sum(n_det > 0),
                     .by = dplyr::all_of(by)) |>
    dplyr::left_join(
      per_field |>
        dplyr::mutate(field = paste0("occ_", janitor_name(field)), occ = n_det / n) |>
        dplyr::select(-n_det, -n) |>
        tidyr::pivot_wider(names_from = field, values_from = occ),
      by = by)
}

# "Gorda Ridge" -> "gorda_ridge"
janitor_name <- function(x) stringr::str_replace_all(tolower(x), "[^a-z0-9]+", "_")


#' Taxonomic span of each KO across its carrier taxa
#'
#' @param carriers Tibble of carrier pairs: one row per KO x carrier taxon with
#'   `ko_id`, `cur_L4` (carrier), `cur_phylum`, `cur_supergroup`, `resolved_L4`,
#'   and `pair_tpm` (the pair's summed vent expression).
#' @param kos Character vector of every KO to report (KOs with no carrier get
#'   n_L4 = 0 and span_rank "no_carrier").
#' @return One row per KO: n_L2, n_L3, n_L4, top carrier and its share, the
#'   share of the KO's carrier TPM in its largest group at each rank
#'   (dom_L2-dom_L4), and `span_rank`, the deepest rank at which ALL carriers
#'   share one taxon (lowest common ancestor). Conservative rule: a carrier not
#'   resolved at a rank (an "(unres.)" class, or an "-Unannotated" phylum) caps
#'   the LCA at the next coarser rank.
#' @export
taxon_span <- function(carriers, kos) {
  span_levels <- c("class", "phylum", "supergroup", "cross_supergroup", "no_carrier")
  out <- carriers |>
    dplyr::mutate(phylum_resolved = !stringr::str_detect(cur_phylum, "-Unannotated$")) |>
    dplyr::summarise(
      n_L2 = dplyr::n_distinct(cur_supergroup),
      n_L3 = dplyr::n_distinct(cur_phylum),
      n_L4 = dplyr::n_distinct(cur_L4),
      top_carrier       = cur_L4[which.max(pair_tpm)],
      top_carrier_share = max(pair_tpm) / sum(pair_tpm),
      dom_L2 = max(tapply(pair_tpm, cur_supergroup, sum)) / sum(pair_tpm),
      dom_L3 = max(tapply(pair_tpm, cur_phylum, sum)) / sum(pair_tpm),
      dom_L4 = max(tapply(pair_tpm, cur_L4, sum)) / sum(pair_tpm),
      span_rank = dplyr::case_when(
        n_L4 == 1 & all(resolved_L4)                       ~ "class",
        n_L3 == 1 & all(phylum_resolved)                    ~ "phylum",
        n_L2 == 1                                           ~ "supergroup",
        TRUE                                                ~ "cross_supergroup"),
      carriers = paste(sort(unique(taxon_label(cur_L4))), collapse = " | "),
      .by = ko_id)
  tibble::tibble(ko_id = kos) |>
    dplyr::left_join(out, by = "ko_id") |>
    dplyr::mutate(dplyr::across(c(n_L2, n_L3, n_L4), \(x) tidyr::replace_na(x, 0L)),
                  span_rank = factor(tidyr::replace_na(span_rank, "no_carrier"),
                                     levels = span_levels))
}


#' Average taxonomic distinctness (Delta+) of each KO's carrier set
#'
#' Clarke & Warwick 1998 (J Appl Ecol 35:523-531), via vegan::taxondive() on
#' presence/absence. Pairwise carrier distance = number of rank steps to their
#' common ancestor over supergroup / phylum / class (equal steps, scaled 0-100
#' by vegan::taxa2dist). The expected value (EDplus) is the mean over all
#' carrier taxa in the pool and does not depend on how many carriers a KO has;
#' only its SD does, so z = (Dplus - EDplus) / sd.Dplus compares KOs with
#' different carrier counts on one scale. z < -2: carriers more taxonomically
#' clustered than random draws from the pool; z > 2: more spread.
#'
#' Conservative handling of unresolved carriers (decision 2026-10-12): an
#' "(unres.)" class or an "-Unannotated" phylum is made unique to that carrier,
#' so it never counts as sharing a class/phylum with another carrier.
#'
#' @param carriers As for [taxon_span()].
#' @return One row per KO with >= 2 carriers: n_carriers, dplus, edplus,
#'   sd_dplus, lambda_plus, dplus_z.
#' @export
distinctness <- function(carriers) {
  pool <- carriers |>
    dplyr::distinct(cur_L4, cur_phylum, cur_supergroup) |>
    dplyr::mutate(phylum = dplyr::if_else(stringr::str_detect(cur_phylum, "-Unannotated$"),
                                          paste(cur_phylum, cur_L4), cur_phylum))
  stopifnot(!anyDuplicated(pool$cur_L4))
  tax <- data.frame(class = pool$cur_L4, phylum = pool$phylum,
                    supergroup = pool$cur_supergroup, row.names = pool$cur_L4)
  tdist <- vegan::taxa2dist(tax, varstep = FALSE, check = FALSE)

  comm <- carriers |>
    dplyr::distinct(ko_id, cur_L4) |>
    dplyr::mutate(present = 1L) |>
    tidyr::pivot_wider(names_from = cur_L4, values_from = present, values_fill = 0L) |>
    tibble::column_to_rownames("ko_id") |>
    as.matrix()
  comm <- comm[, pool$cur_L4, drop = FALSE]
  comm <- comm[rowSums(comm) >= 2, , drop = FALSE]

  td <- vegan::taxondive(comm, tdist)
  tibble::tibble(ko_id       = rownames(comm),
                 n_carriers  = unname(td$Species),
                 dplus       = unname(td$Dplus),
                 edplus      = unname(td$EDplus),
                 sd_dplus    = unname(td$sd.Dplus),
                 lambda_plus = unname(td$Lambda),
                 dplus_z     = (dplus - edplus) / sd_dplus)
}
