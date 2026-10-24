#' Truncate a full lineage string to a fixed taxonomic depth
#'
#' Builds a NESTED taxonomic grain from EUKulele-style `full_classification`
#' strings ("Eukaryota; Excavata; Discoba; Euglenozoa; ..."). Truncating to
#' depth k keeps the first k REAL (non-placeholder) ranks and collapses
#' everything deeper. Because a lineage with fewer than k real ranks returns
#' its full available string (never NA-padded), every transcript keeps a
#' value at every depth, so the transcript set is CONSTANT across the grain
#' series. That constancy is what makes an R2-vs-grain accumulation curve
#' interpretable rather than confounded with a shifting sample of transcripts
#' (Ramond, Galand & Logares 2025, FEMS Microbiol. Rev. 49:fuae031).
#'
#' **Placeholder-padding bug (found 2026-09-04, grain-ladder figure QC;
#' reported by S. Hu as a suspected protist-taxonomy naming artifact and
#' confirmed against the real annotation table).** Checked directly against
#' `data/processed/metat-heterotrophy-annotated.csv` (n = 99,057 transcripts,
#' with `Domain`/`Supergroup`/`Phylum`/`Class`/`Order`/`Family`/`Genus_spp`
#' available as separate columns alongside `full_classification`): each
#' `full_classification` token IS correctly aligned to its named rank column
#' (no skipped or shifted ranks) — but EUKulele/the reference database fills
#' many Order and Family cells with a SYNTHETIC placeholder rather than
#' leaving them `NA`, most commonly `"SG_<Supergroup>_Indeterminate"` or
#' `"<ParentRank>_XX"`, and occasionally by literally repeating the parent
#' rank's own value. These are real database conventions for "unresolved
#' below this point," not typos, and in this table they are common, not
#' rare: of non-`NA` values, **57.8% of Order and 84.5% of Family are
#' synthetic placeholders**, vs. 10.4% of Genus_spp (most real transcripts
#' skip straight from a resolved Class/Order to a resolved genus/species).
#'
#' Naive token counting (the original version of this function) treated a
#' placeholder as a real, distinguishing rank. That fabricates spurious,
#' often near-unique taxon bins at intermediate grains (Order/Family), which
#' then fail the prevalence filter — while the SAME transcripts, once a real
#' deeper rank (Genus_spp) is reached, resolve to shared, often-prevalent
#' bins. That is the mechanism behind the non-monotonic (temporarily RISING)
#' `resolved_frac_tpm` observed between grain 6 and grain 7 in
#' `output-tables/resolution-by-taxdepth.csv` — nesting requires resolved TPM
#' to be non-increasing in depth, and it wasn't, because grain 6 fragmented
#' placeholder-tagged transcripts into rare synthetic bins that grain 7
#' correctly re-pooled onto real genus/species identities.
#'
#' The fix: any token matching a known placeholder pattern
#' (`_XX`/`_Indeterminate`/`_indet` suffix, case-insensitive) or exactly
#' repeating the token immediately before it is dropped, WHEREVER in the
#' string it occurs — not only at the end — before truncating to `depth`.
#' `depth` therefore now counts real informative ranks, not raw string
#' position. This is a pattern match against observed conventions in this
#' project's annotation table, not a ground-truth flag from the annotation
#' pipeline itself (no such flag exists in the source table) — if a future
#' annotation refresh introduces a different placeholder convention, extend
#' `placeholder_re` accordingly and re-run `qc_run_all()`
#' (metat-annotation-qc skill) before trusting a new grain ladder.
#'
#' @param x Character vector of lineage strings.
#' @param depth Integer >= 1: number of REAL (non-placeholder) ranks to keep.
#' @param sep Rank separator; default ";" (surrounding whitespace is trimmed).
#' @return Character vector of depth-truncated lineages (`NA` only when no
#'   real rank parses at all, e.g. an empty or fully-placeholder string).
#' @examples
#'   lineage_depth("Eukaryota; Sar; Alveolata; Dinophyta", 2)
#'   #> "Eukaryota; Sar"
#'
#'   # a mid-string placeholder is dropped, not treated as a real rank, so a
#'   # real deeper rank (Genus_spp) is reachable at its true grain:
#'   lineage_depth(
#'     "Eukaryota; Alveolata; Dinophyta; Dinophyceae; Gonyaulacales; Dinophyceae_XX; Azadinium spinosum 3D9",
#'     depth = 6
#'   )
#'   #> "Eukaryota; Alveolata; Dinophyta; Dinophyceae; Gonyaulacales; Azadinium spinosum 3D9"
#'   #  (only 6 real ranks exist here, so depth = 7 returns the same string)
#' @seealso [build_matrix()], [prevalence_filter()]
#' @export
lineage_depth <- function(x, depth, sep = ";") {
  stopifnot(length(depth) == 1L, depth >= 1L)

  placeholder_re <- stringr::regex("(_XX|_Indeterminate|_indet\\.?)$", ignore_case = TRUE)

  parts <- stringr::str_split(x, stringr::fixed(sep))

  vapply(parts, function(p) {
    p <- stringr::str_trim(p)
    p <- p[nzchar(p) & !is.na(p)]
    if (length(p) == 0L) return(NA_character_)

    is_pattern    <- stringr::str_detect(p, placeholder_re)
    is_dup_of_prev <- c(FALSE, p[-1] == p[-length(p)])
    real <- p[!(is_pattern | is_dup_of_prev)]

    if (length(real) == 0L) return(NA_character_)
    paste(real[seq_len(min(depth, length(real)))], collapse = "; ")
  }, character(1))
}
