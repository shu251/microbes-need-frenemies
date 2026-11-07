# _setup.R — sourced at the top of every .qmd. 
## Loads params, helpers, and sample metadata. 
## Does NOT import or re-run heavy steps (that lives in 01_metat / 02_metat). Underscore prefix = shared helper, not a numbered step.

# Paths in this project are PROJECT-ROOT-RELATIVE. Every qmd starts by sourcing this file (see the locator line in each). Use terminal `quarto render` to render.

# locate the project root, independent of the current working directory
# Root = nearest ancestor directory that holds config/params.yml (also marked by
# _quarto.yml / the .Rproj). Searches upward from the current wd.
.find_proj_root <- function(start = getwd(),
                            markers = c("config/params.yml", "_quarto.yml")) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break            # reached the filesystem root
    d <- parent
  }
  stop("Could not find the project root (looked for config/params.yml upward ",
       "from '", start, "'). Open frenemies-euk-focused.Rproj, or setwd() to the ",
       "project folder, then source this file again.", call. = FALSE)
}

PROJ_ROOT <- .find_proj_root()
if (!identical(normalizePath(getwd(), winslash = "/"),
               normalizePath(PROJ_ROOT, winslash = "/"))) {
  setwd(PROJ_ROOT)
  message("[_setup] working directory set to project root: ", PROJ_ROOT)
}

# --- packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  # core: used in essentially every qmd
  library(tidyverse)   # dplyr / tidyr / ggplot2 / readr / purrr / tibble / stringr / forcats
  library(vegan)       # Aitchison distance, PERMANOVA, betadisper, dbRDA, Procrustes

  # broadly reused across the question qmds
  library(glue)        # readable labels / captions / messages
  library(patchwork)   # compose multi-panel figures (03, 04)
  # broom / janitor / scales dropped 2026-09-27: no bare calls anywhere. 04 uses
  # scales:: namespaced, which needs the package installed but not attached.
})
# NOTE: step-specific Bioconductor packages stay in the step that needs them,
# not here: tximport (01_metat), DESeq2 + apeglm (Q2), arrow (only if the
# annotation contract is delivered as parquet). Keeping them local avoids
# loading heavy deps for the compositional-only questions.

# --- params + helpers (paths are root-relative; wd is now the root) -----------
params <- yaml::read_yaml("config/params.yml")
set.seed(params$seed)

# reusable, tested helper functions
purrr::walk(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

# replicate-aware sample metadata: resolve_sample_unit() adds sample_id (unique
# per physical replicate), sample_unit (biological-unit grouping shared by
# replicates), n_reps, site, and vent_bin. Analyses pick individual replicates
# or averaged units per statistic -- see build_matrix(aggregate=).
sample_meta <- readr::read_csv("config/sample_metadata.csv", show_col_types = FALSE) |>
  resolve_sample_unit()
# sample_meta
