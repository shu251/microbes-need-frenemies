# QC_fungal_signal_check.R — is the Ascomycota (Zymoseptoria tritici) signal in
# library 70 (Von Damm · Ravelin2) contamination, or biology?
#
# Written to close open item #2 in ANALYSIS_REPORT.md §3.5 cheaply. It is a
# DIAGNOSTIC, not an analysis: it does not rebuild a community matrix, does not
# drop anything, and does not touch the Q1-Q5 outputs. House rule stands —
# dominant features are never silently removed; this script only decides how the
# feature gets DESCRIBED in the Methods.
#
# Four independent lines of evidence, cheapest first. Each one is a test that a
# genuine deep-sea fungus would pass and a contaminant would fail (or vice
# versa), so they can disagree — read them together, not individually.
#
#   1. REPLICATE DISCORDANCE. Ravelin2 was sampled twice (libraries 70 and 80,
#      distinct SAMPLEIDs, distinct measured conditions -> true biological
#      replicates). A fungus living in Ravelin2 fluid should appear in both.
#      A reagent / index-hopping artifact belongs to one library only.
#   2. IDENTITY TO THE REFERENCE (max_pid). A wheat pathogen recovered from a
#      hydrothermal vent should be a DIVERGENT relative of Zymoseptoria (low
#      percent identity to the best database hit). Near-identity to the
#      reference genome means the reference organism itself, i.e. contamination.
#   3. TAXONOMIC BREADTH. A fungal bloom is a community: many genera, many
#      contigs, spread across lineages. Contamination from one source is narrow
#      — one genome, a handful of contigs, one library carrying nearly all TPM.
#   4. ASSEMBLY-GROUP PROVENANCE. SEQ_ID carries the assembly group the contig
#      was built in. If the contigs were assembled in a DIFFERENT sample's
#      assembly group than the library where they are abundant, the TPM is
#      arriving by cross-mapping within the co-assembly, which is a mapping
#      artifact rather than an observation about Ravelin2.
#
# Run:  source("scripts/QC_fungal_signal_check.R")   # from the project root
# Out:  output-tables/qc-fungal-*.csv, figures/qc-fungal-signal.png

# --- setup --------------------------------------------------------------------
# _setup.R finds the project root, loads params + R/ helpers, and resolves
# sample_meta. Paths below are project-root-relative, as everywhere else here.
source("scripts/_setup.R")

load("data/processed/long_frenemies_annot_tax.RData")   # long_frenemies_annot (+ _tax, key_taxa)

# The eukaryote-assigned subset is the denominator every downstream matrix uses,
# so the fungal share is expressed against THAT, not against total library TPM
# (see ANALYSIS_REPORT.md §3.3 — 93.4% of TPM has no domain assignment at all).
euk <- long_frenemies_annot |>
  filter(Domain == "Eukaryota", scaledTPM > 0)

# Genus_spp is free text ("Zymoseptoria tritici IPO323"), so match on the genus
# token rather than an exact string.
is_zymo <- function(x) str_detect(coalesce(x, ""), "Zymoseptoria")

# --- 1. replicate discordance -------------------------------------------------
# One row per library: eukaryote TPM, the fungal and Zymoseptoria shares of it,
# and the biological unit each library belongs to. Replicated units (Ravelin2,
# Lotsoshrimp) are the internal control — two libraries, same vent, same fluid.
lib_fungal <- euk |>
  group_by(LIBRARY_NUM, VENT_FIELD, LOCATION, SAMPLETYPE, SAMPLEID, VENT_BIN) |>
  summarise(
    euk_tpm        = sum(scaledTPM),
    fungi_tpm      = sum(scaledTPM[coalesce(Phylum, "") %in% c("Fungi", "Ascomycota", "Basidiomycota") |
                                     coalesce(Class, "") == "Fungi"]),
    zymo_tpm       = sum(scaledTPM[is_zymo(Genus_spp)]),
    n_zymo_contigs = n_distinct(SEQ_ID[is_zymo(Genus_spp)]),
    .groups = "drop"
  ) |>
  mutate(
    fungi_frac = fungi_tpm / euk_tpm,
    zymo_frac  = zymo_tpm  / euk_tpm
  ) |>
  arrange(desc(zymo_frac))

# The replicate pairs, side by side — the headline comparison.
rep_pairs <- lib_fungal |>
  group_by(VENT_FIELD, LOCATION) |>
  filter(n() > 1) |>
  arrange(VENT_FIELD, LOCATION, LIBRARY_NUM) |>
  ungroup()

# --- 2. identity to the reference --------------------------------------------
# max_pid = percent identity to the best database hit. Divergent (novel, deep-sea)
# lineage -> low; the reference organism itself -> high. Weighted by scaledTPM so
# the dominant contig is not averaged away by rare ones.
zymo_contigs <- euk |>
  filter(is_zymo(Genus_spp)) |>
  group_by(SEQ_ID, Genus_spp, full_classification) |>
  summarise(
    tpm_total   = sum(scaledTPM),
    n_libs      = n_distinct(LIBRARY_NUM),
    max_pid     = max(max_pid, na.rm = TRUE),
    .groups = "drop"
  ) |>
  # assembly group = the sample the contig was BUILT in, parsed out of SEQ_ID
  # (megahit_<GROUP>_k111_<node>_flag=..._multi=<cov>_len=<bp>_<group>.fasta)
  mutate(
    assembly_group = str_replace(SEQ_ID, "^megahit_(.*?)_k\\d+_.*$", "\\1"),
    multi_cov      = as.numeric(str_replace(SEQ_ID, "^.*_multi=([0-9.]+)_.*$", "\\1")),
    contig_len     = as.numeric(str_replace(SEQ_ID, "^.*_len=(\\d+).*$", "\\1"))
  ) |>
  arrange(desc(tpm_total))

# For scale: how divergent is Zymoseptoria compared with everything else
# eukaryotic in this dataset? A contaminant should sit at the top of the pid
# distribution, a novel vent lineage near the middle or below.
pid_context <- euk |>
  mutate(grp = if_else(is_zymo(Genus_spp), "Zymoseptoria", "all other Eukaryota")) |>
  group_by(grp) |>
  summarise(
    n_contigs   = n_distinct(SEQ_ID),
    pid_median  = median(max_pid, na.rm = TRUE),
    pid_q25     = quantile(max_pid, 0.25, na.rm = TRUE),
    pid_q75     = quantile(max_pid, 0.75, na.rm = TRUE),
    pid_max     = max(max_pid, na.rm = TRUE),
    .groups = "drop"
  )

# --- 3. taxonomic breadth of the fungal signal -------------------------------
# In the libraries where fungi are abundant: how many genera, and how lopsided?
# A bloom is diverse; a contaminant is one genome.
fungal_breadth <- euk |>
  filter(coalesce(Phylum, "") %in% c("Fungi", "Ascomycota", "Basidiomycota") |
           coalesce(Class, "") == "Fungi") |>
  group_by(LIBRARY_NUM, LOCATION) |>
  summarise(
    fungi_tpm     = sum(scaledTPM),
    n_genera      = n_distinct(Genus_spp, na.rm = TRUE),
    n_contigs     = n_distinct(SEQ_ID),
    top_genus     = Genus_spp[which.max(scaledTPM)],
    top_genus_frac = max(tapply(scaledTPM, coalesce(Genus_spp, "unassigned"), sum)) / sum(scaledTPM),
    .groups = "drop"
  ) |>
  arrange(desc(fungi_tpm))

# --- 4. assembly-group provenance --------------------------------------------
# Where each Zymoseptoria contig was assembled vs. where its TPM lands. Rows off
# the diagonal are cross-mapping within a co-assembly group.
zymo_provenance <- euk |>
  filter(is_zymo(Genus_spp)) |>
  mutate(assembly_group = str_replace(SEQ_ID, "^megahit_(.*?)_k\\d+_.*$", "\\1")) |>
  group_by(assembly_group, LIBRARY_NUM, LOCATION) |>
  summarise(tpm = sum(scaledTPM), n_contigs = n_distinct(SEQ_ID), .groups = "drop") |>
  arrange(desc(tpm))

# --- 5. does it change the composition result? --------------------------------
# ANALYSIS_REPORT.md §3.6 reports Opisthokonta as the largest vent enrichment
# (+10.4 pp), and flags that the number has to be recomputed once the fungal
# signal is resolved. Per the house rule the contigs are NOT dropped from the
# project — this block only reports the contrast BOTH ways so the Methods can
# state the sensitivity.
#
# Two denominators, because they answer different questions and can disagree:
#   pooled  = sum TPM across libraries within a habitat, then take shares
#             (dominated by the deepest / most expressed libraries)
#   by_lib  = share within each library, then average across libraries
#             (each library weighted equally — the right one for a claim about
#             habitats rather than about sequencing effort)
supergroup_shift <- function(dat, label) {
  pooled <- dat |>
    group_by(VENT_BIN, Supergroup) |>
    summarise(tpm = sum(scaledTPM), .groups = "drop_last") |>
    mutate(frac = tpm / sum(tpm)) |>
    ungroup() |>
    select(VENT_BIN, Supergroup, frac) |>
    pivot_wider(names_from = VENT_BIN, values_from = frac, values_fill = 0)

  by_lib <- dat |>
    group_by(VENT_BIN, LIBRARY_NUM, Supergroup) |>
    summarise(tpm = sum(scaledTPM), .groups = "drop_last") |>
    mutate(frac = tpm / sum(tpm)) |>
    group_by(VENT_BIN, Supergroup) |>
    summarise(frac = mean(frac), .groups = "drop") |>
    pivot_wider(names_from = VENT_BIN, values_from = frac, values_fill = 0)

  full_join(pooled, by_lib, by = "Supergroup", suffix = c("_pooled", "_bylib")) |>
    mutate(scenario = label)
}

comp_with    <- supergroup_shift(euk, "with Zymoseptoria")
comp_without <- supergroup_shift(filter(euk, !is_zymo(Genus_spp)), "without Zymoseptoria")

comp_sensitivity <- bind_rows(comp_with, comp_without) |>
  mutate(
    diff_pooled = Vent_pooled - `Non-vent_pooled`,
    diff_bylib  = Vent_bylib  - `Non-vent_bylib`
  ) |>
  select(scenario, Supergroup, Vent_pooled, `Non-vent_pooled`, diff_pooled,
         Vent_bylib, `Non-vent_bylib`, diff_bylib) |>
  arrange(Supergroup, scenario)

# --- 6. figure ----------------------------------------------------------------
# Panel A: the fungal share of every library, so the reader can see that one
#          library sits ~139-fold (about two orders of magnitude) above its own
#          replicate, while the other replicate pair differs by only ~11-fold.
# Panel B: where the dominant contig's expression lands, versus where it was
#          assembled. Log axes (shares span 1e-5 to 0.9), points not bars.
top_contig <- zymo_contigs$SEQ_ID[1]

panel_a <- lib_fungal |>
  mutate(
    lib_lab  = glue("{LIBRARY_NUM} {LOCATION}"),
    lib_lab  = fct_reorder(lib_lab, zymo_frac),
    replicated = LOCATION %in% c("Ravelin2", "Lotsoshrimp"),
    plot_frac = pmax(zymo_frac, 1e-6)          # floor for the log axis; 0 -> 1e-6
  ) |>
  ggplot(aes(x = plot_frac, y = lib_lab)) +
  geom_point(aes(colour = replicated, shape = VENT_BIN), size = 2.6) +
  scale_x_log10(labels = scales::percent_format(accuracy = 0.001)) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "replicated vent", `FALSE` = "single library"),
                      name = NULL) +
  scale_shape_manual(values = c(Vent = 16, `Non-vent` = 1), name = NULL) +
  labs(
    title = "One library carries the fungal signal; its own replicate does not",
    x = expression(italic("Zymoseptoria")~"share of eukaryote-assigned TPM (log)"),
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

panel_b <- euk |>
  filter(SEQ_ID == top_contig) |>
  group_by(LIBRARY_NUM, LOCATION, VENT_FIELD) |>
  summarise(tpm = sum(scaledTPM), .groups = "drop") |>
  mutate(
    lib_lab = glue("{LIBRARY_NUM} {LOCATION}"),
    lib_lab = fct_reorder(lib_lab, tpm),
    assembled_here = LOCATION == "WhiteCastle"      # the contig's assembly group
  ) |>
  ggplot(aes(x = tpm, y = lib_lab)) +
  geom_point(aes(colour = assembled_here), size = 2.6) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "library it was assembled in",
                                 `FALSE` = "other libraries"),
                      name = NULL) +
  labs(
    title = "Assembled at WhiteCastle,\nbut expressed at Ravelin2",
    x = "scaledTPM of the single dominant contig (log)",
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

qc_fig <- panel_a + panel_b +
  patchwork::plot_annotation(
    tag_levels = "a",
    caption = str_wrap(
      glue("Eukaryote-assigned scaledTPM only (2.5% of total library TPM; see ANALYSIS_REPORT.md sec 3.3). ",
           "n = {nrow(lib_fungal)} libraries in the full transcript table (the eukaryote-restricted ",
           "analysis object carries 21; library 69 MustardStand is absent from it). ",
           "Panel a: zero shares plotted at the 1e-6 axis floor; replicated vents are Ravelin2 (70, 80) ",
           "and Lotsoshrimp (68, 78). Panel b: contig {str_trunc(top_contig, 55)}, assembled in the ",
           "VONDAMM_WhiteCastle_vent group."),
      width = 150)
  ) &
  theme(plot.caption = element_text(size = 6.5, hjust = 0))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/qc-fungal-signal.png", plot = qc_fig,
       width = 10, height = 4.6, dpi = 300)

# --- write out ----------------------------------------------------------------
dir.create("output-tables", showWarnings = FALSE)
write_csv(comp_sensitivity, "output-tables/qc-fungal-composition-sensitivity.csv")
write_csv(lib_fungal,      "output-tables/qc-fungal-by-library.csv")
write_csv(zymo_contigs,    "output-tables/qc-fungal-zymo-contigs.csv")
write_csv(fungal_breadth,  "output-tables/qc-fungal-breadth.csv")
write_csv(zymo_provenance, "output-tables/qc-fungal-provenance.csv")
write_csv(pid_context,     "output-tables/qc-fungal-pid-context.csv")

# --- console summary ----------------------------------------------------------
cat("\n== 1. replicated units (the internal control) ==\n")
print(as.data.frame(rep_pairs |>
  select(LIBRARY_NUM, VENT_FIELD, LOCATION, SAMPLEID, euk_tpm, fungi_frac, zymo_frac, n_zymo_contigs)))

cat("\n== 2. identity to reference ==\n"); print(as.data.frame(pid_context))
cat("\ntop Zymoseptoria contigs:\n")
print(as.data.frame(zymo_contigs |>
  select(assembly_group, contig_len, multi_cov, tpm_total, n_libs, max_pid, Genus_spp) |>
  head(8)))

cat("\n== 3. breadth where fungi are abundant ==\n")
print(as.data.frame(head(fungal_breadth, 8)))

cat("\n== 4. assembled-in vs detected-in ==\n")
print(as.data.frame(head(zymo_provenance, 12)))

cat("\n== 5. Opisthokonta enrichment with vs without the fungal contigs ==\n")
print(as.data.frame(comp_sensitivity |> filter(Supergroup == "Opisthokonta")))
