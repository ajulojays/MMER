#!/usr/bin/env Rscript

# MMER publication-grade primary figure
# Primary cohorts: Italy + Manitoba only
# Reads frozen outputs from scripts/22_primary_healthy_twoAims.R and
# scripts/23_PC2_ecological_state_interpretation.R.
# No new hypothesis tests are performed here.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

ROOT <- "."
PRIMARY <- file.path(ROOT, "results/primary_healthy_twoAims")
PC2DIR  <- file.path(ROOT, "results/PC2_ecological_state_interpretation")
OUT     <- file.path(ROOT, "results/publication_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

analysis <- read_csv(file.path(PRIMARY, "tables/PRIMARY_analysis_table.csv"), show_col_types = FALSE)
loadings <- read_csv(file.path(PRIMARY, "tables/baseline_PC_family_loadings.csv"), show_col_types = FALSE)
aim1_coef <- read_csv(file.path(PRIMARY, "tables/Aim1_CR2_coefficients.csv"), show_col_types = FALSE)
aim1_omni <- read_csv(file.path(PRIMARY, "tables/Aim1_omnibus_CR2.csv"), show_col_types = FALSE)
family_summary <- read_csv(file.path(PC2DIR, "tables/PC2_family_by_study_tertile_summary.csv"), show_col_types = FALSE)

analysis <- analysis %>%
  mutate(
    PC2_tertile_num = ntile(z_PC2, 3),
    PC2_tertile = factor(PC2_tertile_num, levels = 1:3,
                         labels = c("Low", "Intermediate", "High")),
    study = factor(study, levels = c("Italy", "Manitoba"))
  )

pc2_row <- aim1_coef %>% filter(term == "z_PC2")
pc2_beta <- pc2_row$beta[[1]]
pc2_p <- pc2_row$p_Satt[[1]]
omnibus_p <- aim1_omni$p_value[[1]]

# Families chosen only from the frozen PC2 loading axis.
top_pos <- loadings %>% arrange(desc(PC2_loading)) %>% slice_head(n = 8)
top_neg <- loadings %>% arrange(PC2_loading) %>% slice_head(n = 8)
selected <- bind_rows(top_neg, top_pos) %>% distinct(family, .keep_all = TRUE)
family_order <- selected %>% arrange(PC2_loading) %>% pull(family)

base_theme <- theme_classic(base_size = 10.5) +
  theme(
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9.2),
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    plot.margin = margin(7, 8, 7, 8)
  )

# Panel A: cohort-level resistance distributions.
pA <- ggplot(analysis, aes(study, ecological_resistance)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.5, fill = "white") +
  geom_jitter(aes(shape = study), width = 0.08, height = 0, size = 2.0, alpha = 0.72) +
  labs(
    x = NULL,
    y = "Ecological resistance\n(1 - Bray-Curtis)",
    title = "A  Quarter-level ecological resistance"
  ) +
  base_theme +
  theme(legend.position = "none")

# Panel B: frozen pooled association with study-specific descriptive lines.
pB <- ggplot(analysis, aes(z_PC2, ecological_resistance, shape = study, linetype = study)) +
  geom_point(size = 2.1, alpha = 0.78, colour = "black") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.72, colour = "black") +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("Pooled PC2: beta = %.3f, CR2 p = %.4f\nPC1 + PC2 omnibus p = %.4f",
                    pc2_beta, pc2_p, omnibus_p),
    hjust = -0.03, vjust = 1.12, size = 3.05
  ) +
  labs(
    x = "Baseline ecological state PC2 (SD)",
    y = "Ecological resistance",
    title = "B  Baseline ecological state is associated with resistance",
    shape = NULL,
    linetype = NULL
  ) +
  base_theme +
  theme(
    legend.position = c(0.82, 0.17),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# Panel C: descriptive low -> high PC2 gradient in each cohort.
tertile_summary <- analysis %>%
  group_by(study, PC2_tertile) %>%
  summarise(
    n = n(),
    mean = mean(ecological_resistance),
    se = sd(ecological_resistance) / sqrt(n),
    .groups = "drop"
  )

pC <- ggplot(tertile_summary,
             aes(PC2_tertile, mean, group = study, shape = study, linetype = study)) +
  geom_line(linewidth = 0.72, colour = "black") +
  geom_point(size = 2.4, colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.5, colour = "black") +
  labs(
    x = "Baseline PC2 ecological state",
    y = "Mean ecological resistance +/- SE",
    title = "C  Resistance across baseline PC2 tertiles",
    shape = NULL,
    linetype = NULL
  ) +
  base_theme +
  theme(legend.position = "top")

# Panel D: descriptive compositional interpretation of the frozen PC2 axis.
# Explicit group order fixes the previous High -> Intermediate -> Low display.
group_levels <- c(
  "Italy\nLow", "Italy\nIntermediate", "Italy\nHigh",
  "Manitoba\nLow", "Manitoba\nIntermediate", "Manitoba\nHigh"
)

heat <- family_summary %>%
  filter(family %in% selected$family) %>%
  mutate(
    family = factor(family, levels = family_order),
    PC2_tertile = factor(
      PC2_tertile,
      levels = c("Low_PC2", "Mid_PC2", "High_PC2"),
      labels = c("Low", "Intermediate", "High")
    ),
    group = factor(paste(study, PC2_tertile, sep = "\n"), levels = group_levels),
    mean_RA_percent = 100 * mean_relative_abundance
  )

pD <- ggplot(heat, aes(group, family, fill = mean_RA_percent)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_viridis_c(
    option = "C",
    trans = "sqrt",
    begin = 0.05,
    end = 0.95
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Mean baseline\nrelative abundance (%)",
    title = "D  Baseline composition along the PC2 ecological-state gradient"
  ) +
  theme_minimal(base_size = 9.3) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.2),
    axis.text.y = element_text(size = 8.0),
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 8, 7, 8)
  )

# Publication layout using only base grid, avoiding extra layout dependencies.
draw_four_panel <- function(filename, device = c("pdf", "svg", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(filename, width = 11.4, height = 9.1, useDingbats = FALSE)
  } else if (device == "svg") {
    svg(filename, width = 11.4, height = 9.1)
  } else {
    png(filename, width = 11.4, height = 9.1, units = "in", res = 600)
  }

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = 2,
    widths = unit(c(0.43, 0.57), "npc"),
    heights = unit(c(0.46, 0.54), "npc")
  )))

  print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  dev.off()
}

draw_four_panel(file.path(OUT, "MMER_primary_ecological_resistance_figure.pdf"), "pdf")
draw_four_panel(file.path(OUT, "MMER_primary_ecological_resistance_figure.svg"), "svg")
draw_four_panel(file.path(OUT, "MMER_primary_ecological_resistance_figure_600dpi.png"), "png")

write_csv(tertile_summary, file.path(OUT, "figure_PC2_tertile_summary.csv"))
write_csv(selected, file.path(OUT, "figure_selected_PC2_families.csv"))

cat("\n[PASS] Publication figure generated\n")
cat("PDF: ", file.path(OUT, "MMER_primary_ecological_resistance_figure.pdf"), "\n", sep = "")
cat("SVG: ", file.path(OUT, "MMER_primary_ecological_resistance_figure.svg"), "\n", sep = "")
cat("PNG: ", file.path(OUT, "MMER_primary_ecological_resistance_figure_600dpi.png"), "\n", sep = "")
