#' Curated, nested taxonomy built from the 02 manual curation
#'
#' 02 (`manual curation of taxa`) adds `SUPERGROUP_0` and `SUPERGROUP_PHYLUM` to
#' `long_frenemies_annot_tax`. Those columns fix naming inconsistencies in the
#' EUKulele reference (Haptista -> Hacrobia, Rhodelphea -> Archaeplastida,
#' incertae sedis -> Picozoa / Telonemia), regroup Stramenopiles into Gyrista and
#' Bigyra (Obiol et al. 2024, FEMS Microbiol Ecol 100:fiae130), and split
#' Alveolata into Dinophyta, Ciliophora, Myzozoa, Syndiniales and Unannotated.
#' Syndiniales is kept as a phylum-level group, parallel to Dinophyta
#' (decision 2026-10-12), so it no longer sits inside the class "Dinophyceae" as
#' it does in the raw EUKulele ranks.
#'
#' This function turns those columns into ONE nested hierarchy, so every taxon at
#' a finer rank belongs to exactly one taxon at each coarser rank:
#'
#'   depth 1  "Eukaryota"
#'   depth 2  curated supergroup   (SUPERGROUP_0 with Alveolata re-merged)
#'   depth 3  curated phylum       (SUPERGROUP_PHYLUM, e.g. "Alveolata-Dinophyta")
#'   depth 4  class-level taxon    (rank-4 name from lineage_depth(), inside its
#'                                  curated phylum; "(unres.)" when the lineage
#'                                  has fewer than 4 real ranks)
#'   depth 5+ deeper real ranks from lineage_depth(), appended
#'
#' Each level is returned as a "; "-joined path, so the last token is the label
#' and the string is unique within its parent. Depth 4 is the analysis grain
#' chosen in 03.
#'
#' @param full_classification,supergroup_0,supergroup_phylum Character vectors of
#'   equal length (columns of the long table).
#' @param depth Integer >= 1.
#' @return Character vector of curated lineage paths truncated to `depth`.
#' @seealso [lineage_depth()], [add_curated_taxonomy()]
#' @export
curated_lineage <- function(full_classification, supergroup_0, supergroup_phylum, depth) {
  stopifnot(length(depth) == 1L, depth >= 1L,
            length(full_classification) == length(supergroup_0),
            length(full_classification) == length(supergroup_phylum))

  key <- tibble::tibble(fc = full_classification, sg0 = supergroup_0, sgp = supergroup_phylum)
  uniq <- dplyr::distinct(key)

  # real (non-placeholder) ranks of the raw lineage, as a list of token vectors
  real_tokens <- stringr::str_split(lineage_depth(uniq$fc, 99L), stringr::fixed("; "))
  n_real      <- lengths(real_tokens)

  supergroup <- stringr::str_split_i(uniq$sg0, "-", 1)
  rank4      <- vapply(real_tokens, \(p) p[min(4L, length(p))], character(1))
  rank4      <- dplyr::if_else(n_real >= 4L, rank4, paste(rank4, "(unres.)"))

  path <- purrr::pmap_chr(
    list(supergroup, uniq$sgp, rank4, real_tokens),
    function(sg, ph, r4, toks) {
      parts <- c("Eukaryota", sg, ph, r4, if (length(toks) > 4L) toks[-(1:4)])
      paste(parts[seq_len(min(depth, length(parts)))], collapse = "; ")
    })

  uniq$path <- path
  dplyr::left_join(key, uniq, by = c("fc", "sg0", "sgp"))$path
}


#' Add the curated rank columns used by 05-07
#'
#' Joins on the distinct lineages only (about 1,400), so it is cheap on the full
#' long table. Adds:
#'   * `cur_supergroup` (depth 2), `cur_group` (= SUPERGROUP_0, the 02 display
#'     grouping used for "summed to supergroup" totals and the 02 palette),
#'     `cur_phylum` (depth 3), `cur_L4` (depth 4), each as a full path;
#'   * `resolved_L4`: TRUE when the raw lineage has >= 4 real ranks.
#'
#' Stops if one lineage maps to more than one curated group (it should not:
#' the 02 curation is a deterministic function of the rank columns).
#'
#' @param df Long table carrying `full_classification`, `SUPERGROUP_0`,
#'   `SUPERGROUP_PHYLUM`.
#' @return `df` with the added columns.
#' @export
add_curated_taxonomy <- function(df) {
  need <- c("full_classification", "SUPERGROUP_0", "SUPERGROUP_PHYLUM")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("add_curated_taxonomy() needs column(s) ", paste(miss, collapse = ", "),
         " from 02's manual curation.", call. = FALSE)
  }
  key <- dplyr::distinct(df, full_classification, SUPERGROUP_0, SUPERGROUP_PHYLUM)
  if (anyDuplicated(key$full_classification)) {
    stop("A full_classification maps to more than one curated group; check 02.",
         call. = FALSE)
  }
  key <- key |>
    dplyr::mutate(
      cur_supergroup = curated_lineage(full_classification, SUPERGROUP_0, SUPERGROUP_PHYLUM, 2),
      cur_group      = SUPERGROUP_0,
      cur_phylum     = curated_lineage(full_classification, SUPERGROUP_0, SUPERGROUP_PHYLUM, 3),
      cur_L4         = curated_lineage(full_classification, SUPERGROUP_0, SUPERGROUP_PHYLUM, 4),
      resolved_L4    = !stringr::str_detect(cur_L4, "\\(unres\\.\\)$")) |>
    dplyr::select(-SUPERGROUP_0, -SUPERGROUP_PHYLUM)
  dplyr::left_join(df, key, by = "full_classification")
}


#' Short display label for a curated path (its last token)
#' @export
taxon_label <- function(path) stringr::str_remove(path, "^.*; ")

#' Display labels that stay unique: a class name shared by two curated phyla
#' (e.g. "Dinophyceae" under Dinophyta and under Syndiniales) gets its phylum
#' appended in brackets.
#' @export
taxon_label_unique <- function(path) {
  u  <- unique(path)
  lu <- taxon_label(u)
  d  <- lu %in% lu[duplicated(lu)]
  ph <- stringr::str_remove(stringr::str_split_i(u[d], "; ", 3), "^[^-]*-")
  lu[d] <- paste0(lu[d], " [", ph, "]")
  lu[match(path, u)]
}
