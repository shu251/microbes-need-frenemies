#' Save and load the shared paired table built at the end of 02
#'
#' 05, 05b, 06 and 07 all built the same table:
#'   long   <- load_long_annot()                       # ~420 MB read
#'   paired <- prep_paired(long, fields = FIELDS_IN_SCOPE, depth = CHOSEN_DEPTH)
#'   rm(long)
#' 02 now builds it once and saves it to `PAIRED_PATH`; the later scripts load
#' that file with `load_paired()` instead of reading the full handoff again.
#'
#' The saved file records which fields and depth it was built with, so
#' `load_paired()` stops if a script asks for a different scope. It also stops
#' if the file is older than the handoff RData or the R code that builds it
#' (`prep_paired.R`, `taxonomy_curated.R`, `lineage_depth.R`). Re-render 02 to
#' rebuild it.

PAIRED_PATH  <- "data/processed/paired_L4.RData"
HANDOFF_PATH <- "data/processed/long_frenemies_annot_tax.RData"
PAIRED_INPUTS <- c(HANDOFF_PATH, "R/prep_paired.R", "R/taxonomy_curated.R",
                   "R/lineage_depth.R")

#' TRUE when the saved paired table is missing or older than any of its inputs.
paired_is_stale <- function(path = PAIRED_PATH, inputs = PAIRED_INPUTS) {
  if (!file.exists(path)) return(TRUE)
  inputs <- inputs[file.exists(inputs)]
  any(file.mtime(inputs) > file.mtime(path))
}

#' Save the paired table with a record of how it was built.
save_paired <- function(paired, fields, depth, path = PAIRED_PATH) {
  paired_info <- list(fields = fields,
                      depth  = as.integer(depth),
                      n_rows = nrow(paired),
                      built  = Sys.time())
  save(paired, paired_info, file = path)
  message("[save_paired] wrote ", path, " (", format(nrow(paired), big.mark = ","),
          " rows; depth ", depth, "; ", paste(fields, collapse = ", "), ")")
  invisible(path)
}

#' Load the paired table, checking it matches the scope the script expects.
#'
#' @param fields,depth The scope the calling script uses (FIELDS_IN_SCOPE,
#'   CHOSEN_DEPTH). Must match what 02 saved.
#' @return The paired table (same as `prep_paired(load_long_annot(), fields, depth)`).
load_paired <- function(fields, depth, path = PAIRED_PATH) {
  if (!file.exists(path)) {
    stop("Paired table not found at '", path, "'.\n",
         "Render scripts/02_taxonomy_grain_ladder.qmd first; its last chunk writes it.",
         call. = FALSE)
  }
  if (paired_is_stale(path)) {
    stop("'", path, "' is older than the handoff RData or the prep code.\n",
         "Re-render scripts/02_taxonomy_grain_ladder.qmd to rebuild it.",
         call. = FALSE)
  }
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  info <- e$paired_info
  if (!setequal(info$fields, fields) || info$depth != as.integer(depth)) {
    stop("Saved paired table was built for depth ", info$depth, " and fields ",
         paste(info$fields, collapse = ", "), ";\n",
         "this script asked for depth ", depth, " and fields ",
         paste(fields, collapse = ", "), ".\n",
         "Change the scope in 02's save-paired chunk, or call ",
         "prep_paired(load_long_annot(), ...) here.", call. = FALSE)
  }
  message("[load_paired] ", path, " (built ", format(info$built, "%Y-%m-%d %H:%M"),
          "; ", format(info$n_rows, big.mark = ","), " rows)")
  e$paired
}
