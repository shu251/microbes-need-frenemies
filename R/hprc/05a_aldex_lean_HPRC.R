## =============================================================================
## 05a — Memory-lean ALDEx2 across all design axes (HPRC)
## Successor to counts-aldex-gamma.R. Same structure, same style, three fixes.
## ALDEx2 >= 1.42.0
## =============================================================================
##
## WHY THE OLD SCRIPT WROTE A 17 GB FILE
## -------------------------------------
## The final line was:
##     save(res, clr, glm_res, gamma_sweep, file = "filtered-glm-results.RData")
##
## `clr` is an aldex.clr S4 object. It holds, for EVERY transcript, EVERY
## library, a full set of Monte Carlo Dirichlet instances:
##
##     bytes ~= n_transcripts * n_libraries * mc.samples * 8
##
## At 207,231 transcripts x 23 libraries x 128 instances that is ~48 GB live in
## RAM and ~17 GB on disk after compression. Everything you actually need from
## it is a few tens of MB. `clr` is scratch. It should never be saved, and it
## should be dropped from memory the moment the quantities derived from it exist.
##
## This script extracts three things while `clr` is alive, then deletes it:
##   1. the tidy GLM table for EVERY model coefficient (all your axes at once),
##   2. ALDEx2 standardised effect sizes for the binary terms,
##   3. a CLR point-estimate matrix (per-transcript, per-library median),
## and writes them as small .rds / .csv.gz files.
##
## Expected output footprint: well under 150 MB total, all of it local-friendly.
##
##
## THREE BUGS / GAPS CARRIED OVER FROM counts-aldex-gamma.R
## --------------------------------------------------------
## (1) The gamma sweep looked for a column matching "vent_bin.*BH$". ALDEx2 names
##     that column `vent_binVent:pval.padj`, so the grep returned NA, `gr[[NA]]`
##     gave NULL, and `sum(NULL < alpha)` is 0. Every gamma row reported n_sig = 0
##     regardless of the data. Fixed here by reusing the same pinned-column
##     helper as the main pass.
## (2) `counts_mat` was built (as.matrix + round) but `aldex.clr()` was then
##     called on `counts`. They happened to be identical, so it worked by luck.
##     One object here, asserted integer.
## (3) Only the vent_bin coefficient was extracted. The `~ field + vent_bin` fit
##     already contains every field contrast; they were being computed and then
##     thrown away. All coefficients are kept here.
##
##
## SIZING YOUR SLURM REQUEST
## -------------------------
## Peak RAM is driven by the Monte Carlo array, not by the counts matrix:
##
##     GB ~= (n_transcripts * n_libraries * mc.samples * 8) / 1e9 * 2.5
##                                                            ^ construction overhead
##
## Worked example, 207k transcripts x 23 libraries:
##     mc.samples = 128  ->  ~122 GB   (needs a bigmem node)
##     mc.samples =  64  ->   ~61 GB
##     mc.samples =  32  ->   ~31 GB
##
## The three levers, in the order you should reach for them:
##   a. Tighten the prevalence filter. Halving n_transcripts halves peak RAM and
##      halves the multiple-testing burden at the same time. Best value by far.
##   b. Lower mc.samples. 128 is the ALDEx2 recommendation for n < 10 per group;
##      at your group sizes 64 is defensible, and the sweep below lets you show
##      the result does not move.
##   c. Ask for a bigmem node.
##
## Do NOT try to chunk the transcripts across jobs. `denom = "all"` computes the
## CLR geometric-mean denominator per library across ALL retained features, so
## splitting the rows silently changes the geometry of every chunk and the
## results are not comparable.
##
## The script refuses to run if predicted peak RAM exceeds `mem_budget_gb`.
## =============================================================================

library(tidyverse)


## ---------------------------------------------------------------------------
## Paths and knobs — the only block you should normally edit
## ---------------------------------------------------------------------------
base_dir <- "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output"

txi_path    <- file.path(base_dir, "tximport-feb-2025.RData")
counts_path <- file.path(base_dir, "counts.RData")          # cached filtered counts
out_dir     <- file.path(base_dir, "05a-aldex-lean")

# Prevalence filter: keep a transcript if it has >= min_count reads in
# >= min_libraries libraries. Raise min_libraries to cut memory; see the
# threshold table this script prints before it commits.
min_count     <- 10
min_libraries <- 3

mc_samples    <- 128        # Monte Carlo Dirichlet instances for the main pass
gamma_main    <- 0.5        # scale uncertainty for the main pass
gamma_grid    <- c(0, 0.25, 0.5, 1.0)
mc_samples_sweep <- 64      # cheaper sweep; only summaries are retained

mem_budget_gb <- 180        # abort rather than get OOM-killed mid-run

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------------------------------------------
## Params — find config/params.yml wherever this is run from
## ---------------------------------------------------------------------------
# The local repo keeps these in config/params.yml. On HPRC the script may sit
# anywhere, so look in a few plausible places before falling back to the
# values that matter here (seed, fdr_alpha).
param_candidates <- c("params.yml",
                      "config/params.yml",
                      "../../config/params.yml")
param_file <- param_candidates[file.exists(param_candidates)][1]

if (!is.na(param_file)) {
  params <- yaml::read_yaml(param_file)
  message(glue::glue("Params read from: {param_file}"))
} else {
  params <- list(seed = 20260730, fdr_alpha = 0.05)
  message("No params.yml found; using inline defaults (seed 20260730, alpha 0.05).")
}

set.seed(params$seed)


## ---------------------------------------------------------------------------
## Counts — reuse the cached filtered matrix if it exists, else build it
## ---------------------------------------------------------------------------
if (file.exists(counts_path)) {

  message(glue::glue("Loading cached counts: {counts_path}"))
  count_env <- new.env(parent = emptyenv())
  load(counts_path, envir = count_env)
  counts_mat <- count_env$counts

} else {

  message(glue::glue("No cache; building counts from {txi_path}"))
  txi_env <- new.env(parent = emptyenv())
  load(txi_path, envir = txi_env)
  counts_mat <- round(txi_env$txi_frenemies_mcr_gr$counts)
  rm(txi_env); gc()

  message(glue::glue("Raw counts: {nrow(counts_mat)} transcripts x ",
                     "{ncol(counts_mat)} libraries."))

  # Threshold table FIRST, so the filter is a decision and not a default.
  # Each row is what you would be left with, and what it would cost in RAM.
  threshold_report <- tidyr::expand_grid(
      min_count     = c(5, 10, 20),
      min_libraries = c(2, 3, 4, 6)
    ) |>
    dplyr::mutate(
      n_retained = purrr::map2_int(
        min_count, min_libraries,
        \(cnt, lib) sum(rowSums(counts_mat >= cnt) >= lib)),
      pct_retained = round(100 * n_retained / nrow(counts_mat), 1),
      peak_gb_mc128 = round(
        n_retained * ncol(counts_mat) * 128 * 8 / 1e9 * 2.5, 1),
      peak_gb_mc64 = round(
        n_retained * ncol(counts_mat) * 64 * 8 / 1e9 * 2.5, 1))

  print(threshold_report, n = Inf)
  readr::write_csv(threshold_report,
                   file.path(out_dir, "prevalence_threshold_report.csv"))

  keep <- rowSums(counts_mat >= min_count) >= min_libraries
  counts_mat <- counts_mat[keep, , drop = FALSE]

  message(glue::glue("Retained {sum(keep)} of {length(keep)} transcripts ",
                     "(>= {min_count} counts in >= {min_libraries} libraries)."))

  counts <- counts_mat
  save(counts, file = counts_path)
  rm(counts)
}

counts_mat <- as.matrix(counts_mat)
storage.mode(counts_mat) <- "double"
counts_mat <- round(counts_mat)

stopifnot(
  all(counts_mat >= 0),
  all(counts_mat == round(counts_mat)),
  !is.null(colnames(counts_mat)),
  !is.null(rownames(counts_mat)))

message(glue::glue("Working matrix: {nrow(counts_mat)} transcripts x ",
                   "{ncol(counts_mat)} libraries."))


## ---------------------------------------------------------------------------
## Design — parsed from the library names, exactly as in counts-aldex-gamma.R
## ---------------------------------------------------------------------------
design <- tibble(sample_rep = colnames(counts_mat)) |>
  separate(sample_rep,
           into = c("meanfield", "LIBRARY_NUM", "fieldyear", "LOCATION",
                    "SAMPLETYPE", "rest"),
           sep = "_", extra = "merge", remove = FALSE) |>
  separate(rest, into = c("SAMPLEID", "replicate"), sep = "_REP", remove = TRUE) |>
  mutate(
    replicate = paste0("REP", replicate),
    field = case_when(
      str_detect(fieldyear, "Piccard") ~ "Piccard",
      str_detect(fieldyear, "VonDamm") ~ "Von Damm",
      str_detect(fieldyear, "Axial")   ~ "Axial",
      str_detect(fieldyear, "Gorda")   ~ "Gorda Ridge"),
    vent_bin = case_when(
      LOCATION %in% c("Background", "Plume", "BSW") ~ "Non-vent",
      str_detect(LOCATION, "IntlDistrict")          ~ "Non-vent",
      str_detect(LOCATION, "ASHES")                 ~ "Non-vent",
      TRUE                                          ~ "Vent"))

# Verify, don't assume (same defensive pattern as tximport_fixed_alignment.R).
stopifnot(
  identical(design$sample_rep, colnames(counts_mat)),
  !anyNA(design$field), !anyNA(design$vent_bin),
  n_distinct(design$sample_rep) == ncol(counts_mat))

message("Design:")
design |> count(field, vent_bin) |> print(n = Inf)

readr::write_csv(design, file.path(out_dir, "design.csv"))

# One model, every axis. `~ field + vent_bin` already estimates each field
# contrast alongside the vent contrast; fitting separate models per axis would
# cost four Monte Carlo runs to recover coefficients this one fit contains.
model_matrix_design <- model.matrix(~ field + vent_bin, data = design)

message("Model coefficients that will be estimated:")
print(colnames(model_matrix_design))


## ---------------------------------------------------------------------------
## Memory guard — predict peak RAM, refuse rather than get OOM-killed
## ---------------------------------------------------------------------------
predict_peak_gb <- function(n_features, n_libraries, mc) {
  n_features * n_libraries * mc * 8 / 1e9 * 2.5
}

peak_gb <- predict_peak_gb(nrow(counts_mat), ncol(counts_mat), mc_samples)

message(glue::glue(
  "Predicted peak RAM at mc.samples = {mc_samples}: {round(peak_gb, 1)} GB ",
  "(budget {mem_budget_gb} GB)."))

if (peak_gb > mem_budget_gb) {
  stop(glue::glue(
    "Predicted peak RAM ({round(peak_gb, 1)} GB) exceeds mem_budget_gb ",
    "({mem_budget_gb} GB).\n",
    "Options, best first: raise min_libraries (see ",
    "prevalence_threshold_report.csv), lower mc_samples, or request a bigmem ",
    "node. Do not chunk transcripts -- denom = 'all' makes chunks ",
    "non-comparable."), call. = FALSE)
}


## ---------------------------------------------------------------------------
## Helper — pull one statistic for one coefficient, by exact column name
## ---------------------------------------------------------------------------
# aldex.glm() names its columns "<coefficient>:<statistic>", e.g.
# "vent_binVent:Est", "vent_binVent:pval.padj". Those names depend on the model
# matrix, so build them explicitly rather than grepping and hoping. This is the
# helper whose absence silently zeroed the old gamma sweep.
glm_column <- function(glm_table, coefficient, statistic) {
  target <- paste0(coefficient, ":", statistic)
  if (!target %in% names(glm_table)) {
    stop("Column '", target, "' not found in aldex.glm output. Available: ",
         paste(head(names(glm_table), 40), collapse = ", "), call. = FALSE)
  }
  glm_table[[target]]
}


## ---------------------------------------------------------------------------
## Main pass — one aldex.clr, everything harvested, then dropped
## ---------------------------------------------------------------------------
message(glue::glue("Running aldex.clr (mc.samples = {mc_samples}, ",
                   "gamma = {gamma_main}). This is the long step."))

clr_main <- ALDEx2::aldex.clr(
  counts_mat,
  model_matrix_design,
  mc.samples = mc_samples,
  denom      = "all",
  gamma      = gamma_main,
  verbose    = FALSE)

message("aldex.clr done. Fitting GLM across all coefficients.")

glm_main <- ALDEx2::aldex.glm(clr_main, model_matrix_design)


## --- HARVEST 1: tidy table of every coefficient ------------------------------
# Long-then-wide so each row is one transcript x one coefficient, with Est / SE /
# pval / pval.padj side by side. The Intercept is dropped -- it is the library
# geometric mean, not a contrast, and it doubles the file for no information.
coefficient_names <- setdiff(colnames(model_matrix_design), "(Intercept)")

# The pivot below splits column names on ":", which is also how R names
# interaction terms. `~ field + vent_bin` has no interactions, so this is safe
# today -- but assert it, because adding one later would silently corrupt the
# term/statistic split rather than error.
stopifnot(!any(grepl(":", coefficient_names, fixed = TRUE)))

de_all_terms <- glm_main |>
  as_tibble(rownames = "transcript_id") |>
  pivot_longer(
    cols = -transcript_id,
    names_to = c("term", "statistic"),
    names_sep = ":",
    values_to = "value") |>
  filter(term %in% coefficient_names) |>
  pivot_wider(names_from = statistic, values_from = value) |>
  janitor::clean_names() |>
  rename(estimate = est, p_value = pval, p_adj = pval_padj) |>
  arrange(term, p_adj)

message("Significant transcripts per coefficient (BH-adjusted):")
de_all_terms |>
  group_by(term) |>
  summarise(
    n_tested  = n(),
    n_sig     = sum(p_adj < params$fdr_alpha, na.rm = TRUE),
    n_up      = sum(p_adj < params$fdr_alpha & estimate > 0, na.rm = TRUE),
    n_down    = sum(p_adj < params$fdr_alpha & estimate < 0, na.rm = TRUE),
    .groups = "drop") |>
  print(n = Inf)


## --- HARVEST 2: standardised effect sizes for the binary terms ---------------
# ALDEx2's own guidance for small n is to rank on effect size rather than on
# p-values alone: `effect` is the between-group difference standardised by the
# larger of the two within-group dispersions, and |effect| > 1 is the usual
# reporting threshold (Fernandes et al. 2014, Microbiome 2:15).
# Guarded: aldex.glm.effect only handles binary columns, and it is not cheap.
effect_sizes <- tryCatch({
  message("Computing standardised effect sizes.")
  raw_effect <- ALDEx2::aldex.glm.effect(clr_main)
  imap_dfr(raw_effect, \(tbl, nm)
    as_tibble(tbl, rownames = "transcript_id") |>
      mutate(term = nm, .after = transcript_id))
}, error = function(e) {
  message("aldex.glm.effect failed (continuing without it): ", e$message)
  NULL
})


## --- HARVEST 3: CLR point estimates, one number per transcript per library ---
# The Monte Carlo instances collapse to a per-library median. This is the object
# that replaces `clr` for everything downstream -- ordination, per-taxon
# aggregation, the Salazar turnover/expression partition -- at roughly
# 1/mc.samples of the size.
message("Collapsing Monte Carlo instances to per-library CLR medians.")

mc_instances <- ALDEx2::getMonteCarloInstances(clr_main)

row_median <- if (requireNamespace("matrixStats", quietly = TRUE)) {
  matrixStats::rowMedians
} else {
  \(m) apply(m, 1, median)
}

clr_median <- vapply(
  mc_instances,
  \(instance) row_median(as.matrix(instance)),
  FUN.VALUE = numeric(nrow(counts_mat)))

rownames(clr_median) <- rownames(counts_mat)
colnames(clr_median) <- names(mc_instances)

stopifnot(
  identical(colnames(clr_median), colnames(counts_mat)),
  identical(rownames(clr_median), rownames(counts_mat)))

message(glue::glue("CLR median matrix: {nrow(clr_median)} x {ncol(clr_median)} ",
                   "({round(object.size(clr_median) / 1e6, 1)} MB)."))


## --- DROP THE WHALE ----------------------------------------------------------
rm(clr_main, mc_instances, glm_main)
gc()
message("aldex.clr object released.")


## ---------------------------------------------------------------------------
## Scale sensitivity — gamma sweep, summaries only
## ---------------------------------------------------------------------------
# Scale models replace normalisation assumptions with an explicit statement of
# how uncertain the unmeasured system scale is; gamma is that uncertainty. The
# sweep shows whether your call survives the assumption (Gloor, Nixon &
# Silverman 2025, NAR Genom Bioinform 7:lqaf108; Dos Santos et al. 2026,
# PLoS Comput Biol 22:e1014728).
#
# One gamma at a time, with an explicit rm()/gc() between, so peak RAM is one
# run's worth rather than four.
message("Starting gamma sweep.")

gamma_sweep <- map_dfr(gamma_grid, function(this_gamma) {

  set.seed(params$seed)
  message(glue::glue("  gamma = {this_gamma}"))

  clr_sweep <- ALDEx2::aldex.clr(
    counts_mat, model_matrix_design,
    mc.samples = mc_samples_sweep, denom = "all",
    gamma = this_gamma, verbose = FALSE)

  glm_sweep <- ALDEx2::aldex.glm(clr_sweep, model_matrix_design)

  # Same pinned-column helper as the main pass. The old sweep grepped for
  # "BH$", which matches nothing, and so reported 0 for every gamma.
  per_term <- map_dfr(coefficient_names, function(this_term) {
    tibble(
      gamma = this_gamma,
      term  = this_term,
      n_sig = sum(glm_column(glm_sweep, this_term, "pval.padj") <
                    params$fdr_alpha, na.rm = TRUE))
  })

  rm(clr_sweep, glm_sweep)
  gc()

  per_term
})

message("Gamma sweep:")
gamma_sweep |> pivot_wider(names_from = term, values_from = n_sig) |> print(n = Inf)


## ---------------------------------------------------------------------------
## Write small artifacts — .rds and .csv.gz, never .RData with a clr in it
## ---------------------------------------------------------------------------
readr::write_csv(gamma_sweep, file.path(out_dir, "gamma_sweep.csv"))

readr::write_csv(de_all_terms,
                 file.path(out_dir, "de_all_terms.csv.gz"))

saveRDS(clr_median,
        file.path(out_dir, "clr_median.rds"), compress = "xz")

if (!is.null(effect_sizes)) {
  readr::write_csv(effect_sizes,
                   file.path(out_dir, "effect_sizes.csv.gz"))
}

saveRDS(
  list(
    n_transcripts    = nrow(counts_mat),
    n_libraries      = ncol(counts_mat),
    min_count        = min_count,
    min_libraries    = min_libraries,
    mc_samples       = mc_samples,
    gamma_main       = gamma_main,
    gamma_grid       = gamma_grid,
    coefficients     = coefficient_names,
    fdr_alpha        = params$fdr_alpha,
    seed             = params$seed,
    aldex2_version   = as.character(utils::packageVersion("ALDEx2")),
    run_date         = Sys.time()),
  file.path(out_dir, "run_manifest.rds"))

message("Wrote:")
tibble(file = list.files(out_dir, full.names = TRUE)) |>
  mutate(mb = round(file.size(file) / 1e6, 2),
         file = basename(file)) |>
  arrange(desc(mb)) |>
  print(n = Inf)

message(glue::glue("Total: {round(sum(file.size(list.files(out_dir, full.names = TRUE))) / 1e6, 1)} MB ",
                   "(was 17000 MB)."))

sessionInfo()
