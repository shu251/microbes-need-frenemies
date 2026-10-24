# ordering schema & color pallettes
yl_red <- c("#ffffcc", "#fed976", "#feb24c", "#fc4e2a", "#e31a1c", "#b10026")

reds <- c("#f7f7f7","#fee5d9", "#fcae91", "#fb6a4a", "#a50f15")

orange_heat <- c("#fff2d4", "#f3d38a", "#ffc94a", "#c76549", "#915e0f", "#62410e")

red_darkred <- c("white","#f2d5d3","#ff7070","#ab333e","darkred", "#410905")

blues <- c("#dde4f0","#b7cae8","#96b3e3",
           "#6b94d6","#3473d9", "#0b4fbd",
           "#0b429c", "#07214a", "#081730")

order_sample_names <- c("Von Damm_Vent_ArrowLoop", "Von Damm_Vent_Bartizan", "Von Damm_Vent_MustardStand", "Von Damm_Vent_OldManTree", "Von Damm_Vent_Ravelin2", "Von Damm_Vent_ShrimpHole", "Von Damm_Vent_WhiteCastle", "Von Damm_Vent_X18",
                        "Piccard_Vent_Lotsoshrimp", "Piccard_Vent_Shrimpocalypse",
                        "Gorda Ridge_Vent_Candelabra", "Gorda Ridge_Vent_MtEdwards", "Gorda Ridge_Vent_SirVentsAlot", "Gorda Ridge_Vent_VentiLatte", 
                        "Von Damm_Non-vent_Plume", "Von Damm_Non-vent_Background",
                        "Piccard_Non-vent_Plume", "Piccard_Non-vent_Background", 
                        "Gorda Ridge_Non-vent_Plume", "Gorda Ridge_Non-vent_BSW")

label_sample_names <- c("Von Damm Arrow Loop", "Von Damm Bartizan", "Von Damm Mustard Stand", "Von Damm Old Man Tree", "Von Damm Ravelin #2", "Von Damm Shrimp Hole", "Von Damm White Castle", "Von Damm X18",
                        "Piccard LotsOshrimp", "Piccard Shrimpocalypse",
                        "Gorda Ridge Candelabra", "Gorda Ridge Mt. Edwards", "Gorda Ridge SirVentsAlot", "Gorda Ridge Venti Latte", 
                        "Von Damm Plume", "Von Damm Background",
                        "Piccard Plume", "Piccard Background", 
                        "Gorda Ridge Plume", "Gorda Ridge Background")

label_sample_names_SHORT <- c("Arrow Loop", "Bartizan", "Mustard Stand", "Old Man Tree", "Ravelin #2", "Shrimp Hole", "White Castle", "X18",
                              "LotsOshrimp", "Shrimpocalypse",
                              "Candelabra", "Mt. Edwards", "SirVentsAlot", "Venti Latte", 
                              "Plume", "Background",
                              "Plume", "Background", 
                              "Plume", "Background")

names(label_sample_names) <- order_sample_names
names(label_sample_names_SHORT) <- order_sample_names

order_field <- c("Von Damm", "Piccard", "Gorda Ridge")
color_field <- c("#973c46", "#b5682d", "#6a89ac")

names(color_field) <- order_field

order_bin <- c("Vent", "Non-vent")
color_bin <- c("#b74a70", "#346D6E")
names(color_bin) <- order_bin

# Curated taxonomic groups (SUPERGROUP_0 from 02's manual curation) and their colors.
# Moved here from 02 on 2026-10-12 so 03, 05 and 06 draw from one palette.
SUPERGROUP_0_ordering <- c("Alveolata-Dinophyta", "Alveolata-Myzozoa", "Alveolata-Syndiniales", "Alveolata-Ciliophora", "Alveolata-Unannotated",
                           "Telonemia", "Stramenopiles", "Rhizaria",
                           "Breviatea", "Excavata", "Hacrobia", "CRuMs",
                           "Apusozoa", "Amoebozoa", "Picozoa",
                           "Ancyromonadida", "Archaeplastida",
                           "Opisthokonta",
                           "Eukaryota incertae sedis", "Eukaryote")

tax_color_metaT <- c("#1abc9c", "#34495e", "#9b59b6", "#fea02f", "#002663","#ff6f61", "#189ab4", "#f5bfb0", "#3f9c2a","#de6600", "#007a7a", "#93c47d","#ee429b", "#d7b6d7", "#e71414", "#ffe8a4", "#753c71","#2986cc", "#bcbcbc", "#444444")

names(tax_color_metaT) <- SUPERGROUP_0_ordering
