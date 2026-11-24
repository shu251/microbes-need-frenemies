# 05_habitat_class_schematic.R -- standalone figure, run from the project root:
#   Rscript scripts/05_habitat_class_schematic.R
# Writes figures/05-habitat-class-schematic.svg. Not part of the Quarto site.
# Schematic of the habitat classes (05/05b rule, 2026-11-23). One example feature per class.
library(tidyverse)

class_lab <- c(vent_only = "Vent-only", vent_skewed = "Vent-skewed", shared = "Shared",
               nonvent_skewed = "Non-vent-skewed", nonvent_only = "Non-vent-only",
               below_filter = "Below filter")
class_rule <- c(vent_only = "≥3 vent, 0 background (plume OK)",
                vent_skewed = "≥3 vent, ≥1 background, <3 non-vent",
                shared = "≥3 vent and ≥3 non-vent",
                nonvent_skewed = "1–2 vent, ≥3 non-vent",
                nonvent_only = "0 vent, ≥3 non-vent",
                below_filter = "<3 vent and <3 non-vent")
# counts from the 2026-11-23 rule (KOs / taxon x KO pairs); update after re-rendering 05/05b
class_n <- c(vent_only = "84 KOs\n279 pairs", vent_skewed = "45\n140", shared = "1,126\n820",
             nonvent_skewed = "1,180\n1,778", nonvent_only = "906\n2,401", below_filter = "")
class_pal <- c(vent_only = "#7A1F45", vent_skewed = "#D98AA8", shared = "grey75",
               nonvent_skewed = "#7FB3B4", nonvent_only = "#1F5556", below_filter = "grey90")
unit_pal <- c(Vent = "#b74a70", Plume = "#7FB3B4", Background = "#346D6E")

# example detections: vent units, plume units, background units
example <- tribble(
  ~class,           ~n_vent, ~n_plume, ~n_bg,
  "vent_only",           5,        1,     0,
  "vent_skewed",         5,        0,     1,
  "shared",              6,        2,     2,
  "nonvent_skewed",      2,        2,     2,
  "nonvent_only",        0,        2,     2,
  "below_filter",        2,        1,     0)

dots <- example %>%
  mutate(units = pmap(list(n_vent, n_plume, n_bg), \(v, p, b) tibble(
    unit_type = rep(c("Vent", "Plume", "Background"), c(13, 3, 3)),
    idx = c(1:13, 1:3, 1:3),
    detected = c(seq_len(13) <= v, seq_len(3) <= p, seq_len(3) <= b)))) %>%
  unnest(units) %>%
  mutate(x = case_when(unit_type == "Vent" ~ idx,
                       unit_type == "Plume" ~ 14.5 + idx,
                       TRUE ~ 18 + idx),
         class = factor(class, levels = names(class_lab)),
         y = as.numeric(fct_rev(class)))

rows <- distinct(dots, class, y) %>%
  mutate(label = class_lab[as.character(class)], rule = class_rule[as.character(class)],
         n = class_n[as.character(class)])
heads <- tibble(x = c(7, 16, 19.5), lab = c("Vent units (13)", "Plume (3)", "Background (3)"))

fig_class_schematic <- ggplot(dots, aes(x, y)) +
  geom_tile(data = rows, aes(x = -10.3, y = y, fill = class), width = 0.45, height = 0.7,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = filter(dots, !detected), shape = 21, size = 4.2, colour = "grey70",
             fill = "white", stroke = 0.5) +
  geom_point(data = filter(dots, detected), aes(colour = unit_type), size = 4.5) +
  geom_vline(xintercept = 14.25, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_text(data = rows, aes(x = -9.8, y = y + 0.13, label = label), hjust = 0,
            fontface = "bold", size = 3.6) +
  geom_text(data = rows, aes(x = -9.8, y = y - 0.2, label = rule), hjust = 0,
            size = 2.9, colour = "grey30") +
  geom_text(data = rows, aes(x = 22.6, y = y, label = n), hjust = 0, size = 3, lineheight = 0.9) +
  geom_text(data = heads, aes(x = x, y = 6.75, label = lab), size = 3.4, fontface = "bold") +
  annotate("text", x = 23.3, y = 6.75, label = "KOs\npairs", size = 3, lineheight = 0.9,
           fontface = "bold") +
  annotate("text", x = 14.5, y = 0.35, hjust = 0, label = "plume detections allowed for vent-only",
           size = 2.7, colour = "grey30") +
  scale_fill_manual(values = class_pal) +
  scale_colour_manual(values = unit_pal, breaks = names(unit_pal), name = "Detected in") +
  scale_x_continuous(limits = c(-10.8, 25), expand = expansion(0)) +
  scale_y_continuous(limits = c(0.2, 7), expand = expansion(0)) +
  labs(caption = paste("Associated = detected in ≥3 units of a habitat",
                       "(3 of 13 vent; 3 of 6 non-vent, plume + background). Open circle = not detected.")) +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, colour = "grey30", size = 8.5),
        plot.margin = margin(8, 8, 8, 8))

ggsave("figures/05-habitat-class-schematic.svg", fig_class_schematic, width = 9.5, height = 4.8)
