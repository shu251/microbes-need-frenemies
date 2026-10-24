#' KEGG BRITE ko00001: KO -> level-B category map
#'
#' Parses the KEGG BRITE KO hierarchy (`ko00001.keg`) into one row per KO x
#' level-B category. Line prefixes: A = top class, B = level-B category,
#' C = pathway, D = KO. A and B labels are carried down onto each KO row.
#' A KO in two pathways of the same level-B category counts once there;
#' `n_cat` is the number of distinct level-B categories per KO, used to split
#' a KO's weight equally across them (same mass-conserving rule as 04).
#'
#' Shared by 04 (functional investment, KEGG tier) and 06 (breadth-class
#' composition), so both read the hierarchy identically.
#'
#' @param path Path to ko00001.keg (download once:
#'   curl -o data/external/ko00001.keg "https://rest.kegg.jp/get/br:ko00001").
#' @param exclude_a Top-class ids dropped before mapping. Default drops
#'   09150 Organismal Systems, 09160 Human Diseases, 09180 BRITE hierarchies,
#'   09190 Not included in pathway or BRITE (as in 04).
#' @return Tibble: ko_id, ko_name, kegg_a, a_description, kegg, description, n_cat.
#' @export
read_brite_level_b <- function(path = "data/external/ko00001.keg",
                               exclude_a = c("09150", "09160", "09180", "09190")) {
  if (!file.exists(path)) {
    stop("`", path, "` not found. Download it once with\n",
         "  curl -o ", path, " \"https://rest.kegg.jp/get/br:ko00001\"", call. = FALSE)
  }
  tibble::tibble(line = readr::read_lines(path)) |>
    dplyr::filter(stringr::str_detect(line, "^[ABD]")) |>
    dplyr::mutate(level  = stringr::str_sub(line, 1, 1),
                  text   = stringr::str_squish(stringr::str_sub(line, 2)),
                  a_text = dplyr::if_else(level == "A", text, NA_character_),
                  b_text = dplyr::if_else(level == "B" & text != "", text, NA_character_)) |>
    tidyr::fill(a_text, b_text) |>
    dplyr::filter(level == "D") |>
    dplyr::transmute(
      ko_id         = stringr::str_extract(text, "^K\\d{5}"),
      ko_name       = stringr::str_squish(stringr::str_remove(text, "^K\\d{5}")),
      kegg_a        = stringr::str_extract(a_text, "^\\d{5}"),
      a_description = stringr::str_squish(stringr::str_remove(a_text, "^\\d{5}")),
      kegg          = stringr::str_extract(b_text, "^\\d{5}"),
      description   = stringr::str_squish(stringr::str_remove(b_text, "^\\d{5}"))) |>
    dplyr::filter(!is.na(ko_id), !is.na(kegg), !kegg_a %in% exclude_a) |>
    dplyr::distinct(ko_id, kegg, .keep_all = TRUE) |>
    dplyr::mutate(n_cat = dplyr::n(), .by = ko_id)
}

#' KEGG BRITE ko00001: KO -> level-C pathway map
#'
#' Same file and line grammar as `read_brite_level_b()`, one level down: one row
#' per KO x pathway (C line), with its level-B category carried along. Only C
#' lines that are KEGG pathway maps (`[PATH:ko#####]`) are kept, so BRITE-only
#' groupings never become gene sets.
#'
#' @param path Path to ko00001.keg.
#' @param keep_b Level-B ids to keep (default: all).
#' @return Tibble: ko_id, ko_name, kegg_b, b_description, pathway_id, pathway_name.
#' @export
read_brite_level_c <- function(path = "data/external/ko00001.keg", keep_b = NULL) {
  if (!file.exists(path)) {
    stop("`", path, "` not found. Download it once with\n",
         "  curl -o ", path, " \"https://rest.kegg.jp/get/br:ko00001\"", call. = FALSE)
  }
  tibble::tibble(line = readr::read_lines(path)) |>
    dplyr::filter(stringr::str_detect(line, "^[BCD]")) |>
    dplyr::mutate(level  = stringr::str_sub(line, 1, 1),
                  text   = stringr::str_squish(stringr::str_sub(line, 2)),
                  b_text = dplyr::if_else(level == "B" & text != "", text, NA_character_),
                  c_text = dplyr::if_else(level == "C", text, NA_character_),
                  # a new B line starts a new block: blank its C so a stale
                  # pathway is never carried into the next category
                  c_text = dplyr::if_else(level == "B", "", c_text)) |>
    tidyr::fill(b_text, c_text) |>
    dplyr::filter(level == "D", c_text != "") |>
    dplyr::transmute(
      ko_id         = stringr::str_extract(text, "^K\\d{5}"),
      ko_name       = stringr::str_squish(stringr::str_remove(text, "^K\\d{5}")),
      kegg_b        = stringr::str_extract(b_text, "^\\d{5}"),
      b_description = stringr::str_squish(stringr::str_remove(b_text, "^\\d{5}")),
      pathway_id    = stringr::str_extract(c_text, "(?<=\\[PATH:)ko\\d{5}(?=\\])"),
      pathway_name  = stringr::str_squish(stringr::str_remove_all(c_text, "^\\d{5}|\\[PATH:ko\\d{5}\\]"))) |>
    dplyr::filter(!is.na(ko_id), !is.na(pathway_id),
                  is.null(keep_b) | kegg_b %in% keep_b) |>
    dplyr::distinct(ko_id, pathway_id, .keep_all = TRUE)
}
