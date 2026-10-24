#' Bin samples into Vent vs Non-vent
#'
#' The single, authoritative vent/non-vent definition for the whole project.
#' Previously this logic was re-implemented inline in `metaT-preprocess-all.R`
#' and `vent-metaT-DataAnalysis.qmd`, which let the binning drift between
#' analyses. Source this one function everywhere instead.
#'
#' @param habitat Character vector of habitat labels (e.g. the third token of
#'   `SAMPLE_NAME`, such as "Background", "Plume", "MtEdwards").
#' @param nonvent_tokens Character vector of habitat tokens that count as
#'   Non-vent. Defaults to the project config value.
#'
#' @return A factor with levels `c("Non-vent", "Vent")`, aligned to `habitat`.
#' @examples
#'   bin_vent(c("Background", "Plume", "MtEdwards"))
#'   #> Non-vent Non-vent Vent
#' @export
bin_vent <- function(habitat,
                     nonvent_tokens = c("Background", "Plume", "BSW",
                                        "IntlDistrict", "ASHES")) {
  stopifnot(is.character(habitat))
  is_nonvent <- habitat %in% nonvent_tokens
  factor(ifelse(is_nonvent, "Non-vent", "Vent"),
         levels = c("Non-vent", "Vent"))
}
