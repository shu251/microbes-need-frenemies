## Script information
# ALDEx2 1.42.0


# library set up
library(tidyverse)

params <- yaml::read_yaml("params.yml")
set.seed(params$seed)


# Count-level txi from HPRC (bundle as RData, matching the 02 handoff convention).
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-feb-2025.RData",verbose = TRUE)

names(txi_frenemies_mcr_gr)
dim(txi_frenemies_mcr_gr$counts)
head(colnames(txi_frenemies_mcr_gr$counts))
counts <- round(txi_frenemies_mcr_gr$counts)   # integer counts; cols = sample_id
message(glue::glue("Counts matrix: {nrow(counts)} transcripts x {ncol(counts)} libraries."))

# PREVALENCE FILTER
keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, , drop = FALSE]

message(glue::glue("Retained {sum(keep)} of {length(keep)} transcripts ",
                   "(>=10 counts in >=3 libraries)."))

save(counts, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/counts.RData")

#
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/counts.RData", verbose = TRUE)

# Make matrix
counts_mat <- counts |>
  as.matrix() |>
  round()  

# Get design
design <- tibble(sample_rep = colnames(counts_mat)) |>
  separate(sample_rep,
           into = c("meanfield", "LIBRARY_NUM", "fieldyear", "LOCATION", "SAMPLETYPE", "rest"),
           sep = "_", extra = "merge", remove = FALSE) |>
  separate(rest, into = c("SAMPLEID", "replicate"), sep = "_REP", remove = TRUE) |>
  mutate(
    replicate = paste0("REP", replicate),
    field = case_when(
      str_detect(fieldyear, "Piccard")  ~ "Piccard",
      str_detect(fieldyear, "VonDamm")  ~ "Von Damm",
      str_detect(fieldyear, "Axial")~ "Axial",
      str_detect(fieldyear, "Gorda") ~ "Gorda Ridge"),
    vent_bin = case_when(
      LOCATION %in% c("Background", "Plume", "BSW") ~ "Non-vent",
      str_detect(LOCATION, "IntlDistrict") ~ "Non-vent",
      str_detect(LOCATION, "ASHES") ~ "Non-vent",
      TRUE~ "Vent"))

# same defensive pattern as tximport_fixed_alignment.R: verify, don't assume
stopifnot(
  identical(design$sample_rep, colnames(counts_mat)),
  !anyNA(design$field), !anyNA(design$vent_bin),
  n_distinct(design$sample_rep) == ncol(counts_mat))

# Print design
design |> count(vent_bin, field)

#
# library(ALDEx2)
mm <- model.matrix(~ field + vent_bin, data = design) #stats package

## Run
clr <- ALDEx2::aldex.clr(counts, mm, mc.samples = 128, denom = "all", gamma = 0.5, verbose = FALSE)

glm_res <- ALDEx2::aldex.glm(clr, mm)

glimpse(glm_res)

names(glm_res)   # eyeball once

vent_cols <- grep("vent_binVent", names(glm_res), value = TRUE)
vent_cols 

col_of <- function(suffix) {
  hit <- grep(suffix, vent_cols, value = TRUE)
  if (length(hit) != 1) {
    stop("Expected exactly 1 match for '", suffix, "' in vent_binVent columns, got ",
         length(hit), ": ", paste(hit, collapse = ", "))
  }
  hit
}

res <- glm_res |>
  as_tibble(rownames = "transcript_id") |>
  transmute(transcript_id,
            effect = .data[[col_of("Est$")]],
            se = .data[[col_of("SE$")]],
            pval = .data[[col_of("^vent_binVent:pval$")]],
            pad= .data[[col_of("pval\\.padj$")]]) |>
  arrange(desc(effect))

res |> 
  summarise(n_tested = n(), 
            n_sig = sum(pad < params$fdr_alpha, na.rm = TRUE))


glimpse(res)


# Scale sensitivity
gamma_sweep <- purrr::map_dfr(c(0, 0.25, 0.5, 1.0), function(g) {
  set.seed(params$seed)
  clr_g <- ALDEx2::aldex.clr(counts, mm, mc.samples = 64, denom = "all", gamma = g, verbose = FALSE)
  gr    <- ALDEx2::aldex.glm(clr_g, mm)
  padj  <- gr[[grep("vent_bin.*BH$", names(gr))[1]]]
  tibble::tibble(gamma = g, n_sig = sum(padj < params$fdr_alpha, na.rm = TRUE))
})

gamma_sweep

save(res, clr, glm_res, gamma_sweep, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/filtered-glm-results.RData")


#
sessionInfo()