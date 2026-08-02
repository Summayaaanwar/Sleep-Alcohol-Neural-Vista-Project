# ============================================================================
# Separate publication figures for the alcohol-use and sleep study
# NHANES 2017-2018
#
# This script DOES NOT rerun any statistical models. It reads the final result
# CSV files and presents the supplied design-based estimates, 95% confidence
# intervals, and p values in publication-ready figures.
#
# Output: ten separate figure/panel files, each saved as a 600-dpi PNG and
# vector PDF.
# The figures can later be combined into multipanel figures in Word, PowerPoint,
# Adobe Illustrator, Inkscape, or with an R package such as patchwork.
# ============================================================================

rm(list = ls())

# ---- 1. Packages ------------------------------------------------------------
required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "readr",
  "scales", "stringr", "tibble", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# ---- 2. Project folders -----------------------------------------------------
# Forward slashes are intentional and work correctly in Windows R.
default_base_dir <- "E:/Neural Vista/projct alcohol and sleep"

# Environment-variable overrides are optional. They make the script portable.
base_dir <- Sys.getenv("SLEEP_ALCOHOL_PROJECT_DIR", unset = default_base_dir)

results_override <- Sys.getenv("SLEEP_ALCOHOL_RESULTS_DIR", unset = "")
figures_dir <- Sys.getenv(
  "SLEEP_ALCOHOL_FIGURES_DIR",
  unset = file.path(base_dir, "outputs", "figures", "separate_panels")
)

# The code checks the most likely result-table locations automatically.
result_candidate_dirs <- unique(c(
  results_override,
  file.path(base_dir, "outputs", "tables"),
  file.path(base_dir, "outputs", "Tables with p value"),
  file.path(base_dir, "outputs", "tables", "For p values"),
  file.path(base_dir, "outputs"),
  base_dir
))
result_candidate_dirs <- result_candidate_dirs[nzchar(result_candidate_dirs)]

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Find either the normal filename or an uploaded-copy name such as "(2).csv".
find_result_file <- function(file_name) {
  exact_paths <- file.path(result_candidate_dirs, file_name)
  exact_match <- exact_paths[file.exists(exact_paths)]

  if (length(exact_match) > 0) {
    message("Using: ", exact_match[[1]])
    return(exact_match[[1]])
  }

  # The required stems contain only letters and underscores, so a simple,
  # transparent pattern is safer than complicated regex escaping.
  stem <- tools::file_path_sans_ext(file_name)
  pattern <- paste0("^", stem, "(\\([0-9]+\\))?\\.csv$")

  existing_dirs <- result_candidate_dirs[dir.exists(result_candidate_dirs)]
  all_csv <- unlist(lapply(
    existing_dirs,
    list.files,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  ), use.names = FALSE)

  copy_match <- all_csv[str_detect(basename(all_csv), regex(pattern, ignore_case = TRUE))]

  if (length(copy_match) == 0) {
    stop(
      "Could not find ", file_name, ". Checked these folders:\n",
      paste0(" - ", result_candidate_dirs, collapse = "\n")
    )
  }

  message("Using: ", copy_match[[1]])
  copy_match[[1]]
}

read_result <- function(file_name) {
  read_csv(find_result_file(file_name), show_col_types = FALSE)
}

# ---- 3. Read the final result files ----------------------------------------
flow_counts <- read_result("flow_counts.csv")
prevalence <- read_result("sleep_trouble_prevalence_by_alcohol.csv")
prevalence_pairwise <- read_result("sleep_trouble_prevalence_pairwise_pvalues.csv")
primary <- read_result("primary_sleep_trouble_models.csv")
primary_overall <- read_result("primary_alcohol_overall_pvalues.csv")
secondary <- read_result("secondary_binary_models.csv")
sleep_hours <- read_result("sleep_hours_model.csv")
dose_response <- read_result("dose_response_models.csv")
sensitivity <- read_result("sensitivity_models.csv")
interaction <- read_result("interaction_models.csv")
interaction_joint <- read_result("interaction_joint_pvalue.csv")
spearman_long <- read_result("spearman_correlation_long.csv")

# ---- 4. Basic validation ----------------------------------------------------
check_columns <- function(data, columns, source_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(source_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

check_columns(flow_counts, c("step", "n"), "flow_counts.csv")
check_columns(prevalence,
              c("alcohol_category", "prevalence", "CI_low", "CI_high", "p_overall"),
              "sleep_trouble_prevalence_by_alcohol.csv")
check_columns(prevalence_pairwise, c("comparison", "OR", "CI_low", "CI_high", "p"),
              "sleep_trouble_prevalence_pairwise_pvalues.csv")
check_columns(primary, c("model", "n", "term", "OR", "CI_low", "CI_high", "p"),
              "primary_sleep_trouble_models.csv")
check_columns(primary_overall, c("model", "n", "p_overall_alcohol"),
              "primary_alcohol_overall_pvalues.csv")
check_columns(secondary, c("outcome", "n", "term", "OR", "CI_low", "CI_high", "p"),
              "secondary_binary_models.csv")
check_columns(sleep_hours, c("n", "term", "beta", "CI_low", "CI_high", "p"),
              "sleep_hours_model.csv")
check_columns(dose_response, c("outcome", "n", "term", "effect", "CI_low", "CI_high", "p"),
              "dose_response_models.csv")
check_columns(sensitivity, c("outcome", "n", "term", "OR", "CI_low", "CI_high", "p"),
              "sensitivity_models.csv")
check_columns(interaction, c("n", "term", "OR", "CI_low", "CI_high", "p"),
              "interaction_models.csv")
check_columns(interaction_joint, c("n", "p"), "interaction_joint_pvalue.csv")
check_columns(spearman_long, c("variable_1", "variable_2", "rho", "p"),
              "spearman_correlation_long.csv")

# ---- 5. Appearance and helper functions ------------------------------------
base_font <- "Times New Roman"

# Colour-blind-friendly palette.
colour_non_binge <- "#0072B2"  # blue
colour_binge <- "#D55E00"      # vermillion
colour_neutral <- "#4D4D4D"    # dark grey
colour_accent <- "#6A3D9A"     # purple
colour_grid <- "#D8D8D8"
colour_band <- "#F2F5F7"

p_text <- function(p) {
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

number_2 <- function(x) formatC(x, format = "f", digits = 2)
number_3 <- function(x) formatC(x, format = "f", digits = 3)

or_result_text <- function(or, low, high, p) {
  paste0(
    "OR ", number_2(or), " (", number_2(low), "-", number_2(high), "); ",
    p_text(p)
  )
}

beta_result_text <- function(beta, low, high, p, digits = 2) {
  fmt <- function(x) formatC(x, format = "f", digits = digits)
  paste0(
    "beta ", fmt(beta), " (", fmt(low), " to ", fmt(high), "); ",
    p_text(p)
  )
}

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = base_font) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(face = "bold", size = base_size, hjust = 0,
                                   margin = margin(b = 9)),
      plot.caption = element_text(size = base_size - 2, hjust = 0,
                                  margin = margin(t = 9)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(face = "bold", colour = "black", size = base_size - 1),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = base_size - 1),
      legend.position = "bottom",
      legend.justification = "left",
      # panel.border frames the graph itself; plot.background frames the
      # complete exported image, including its title, legend, and caption.
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.9),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
}

save_figure <- function(plot_object, file_stem, width, height) {
  png_path <- file.path(figures_dir, paste0(file_stem, ".png"))
  pdf_path <- file.path(figures_dir, paste0(file_stem, ".pdf"))

  ggsave(
    png_path, plot = plot_object,
    width = width, height = height, units = "in",
    dpi = 600, bg = "white"
  )

  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(
    pdf_path, plot = plot_object,
    width = width, height = height, units = "in",
    device = pdf_device, bg = "white"
  )

  message("Saved: ", png_path)
  message("Saved: ", pdf_path)
}

# Add horizontal confidence intervals with small vertical end caps.
add_horizontal_ci <- function(plot_object, colour_column = NULL, cap_height = 0.08) {
  if (is.null(colour_column)) {
    plot_object +
      geom_segment(aes(x = CI_low, xend = CI_high, y = y, yend = y),
                   linewidth = 0.8, colour = colour_neutral) +
      geom_segment(aes(x = CI_low, xend = CI_low,
                       y = y - cap_height, yend = y + cap_height),
                   linewidth = 0.8, colour = colour_neutral) +
      geom_segment(aes(x = CI_high, xend = CI_high,
                       y = y - cap_height, yend = y + cap_height),
                   linewidth = 0.8, colour = colour_neutral)
  } else {
    plot_object +
      geom_segment(aes(x = CI_low, xend = CI_high, y = y, yend = y,
                       colour = .data[[colour_column]]), linewidth = 0.8) +
      geom_segment(aes(x = CI_low, xend = CI_low,
                       y = y - cap_height, yend = y + cap_height,
                       colour = .data[[colour_column]]), linewidth = 0.8) +
      geom_segment(aes(x = CI_high, xend = CI_high,
                       y = y - cap_height, yend = y + cap_height,
                       colour = .data[[colour_column]]), linewidth = 0.8)
  }
}

# ---- 6. Figure 1: methodological and analytical workflow -------------------
# This is a METHODS figure. It summarizes the data-preparation and analysis
# sequence; inferential p values are not applicable to a workflow diagram.

flow_lookup <- setNames(flow_counts$n, flow_counts$step)
required_flow_steps <- c(
  "N_adults", "N_positive_MEC_weight", "N_non_missing_outcome",
  "N_non_missing_exposure", "N_complete_model4"
)

if (!all(required_flow_steps %in% names(flow_lookup))) {
  stop("flow_counts.csv does not contain all rows required for Figure 1.")
}

workflow_nodes <- tibble(
  x = 0.5,
  y = 7:1,
  node_group = c(
    "Data", "Design", "Variables", "Primary",
    "Secondary", "Robustness", "Exploratory"
  ),
  label = c(
    paste0(
      "NHANES 2017-2018 component files merged by SEQN\n",
      "Adults aged 18 years or older (n = ",
      comma(flow_lookup[["N_adults"]]), ")"
    ),
    paste0(
      "Complex survey design specified\n",
      "Positive MEC weight, strata, and primary sampling units (n = ",
      comma(flow_lookup[["N_positive_MEC_weight"]]), ")"
    ),
    paste0(
      "Exposure, sleep outcomes, and covariates constructed\n",
      "Nonmissing primary outcome after weight restriction (n = ",
      comma(flow_lookup[["N_non_missing_outcome"]]), ")"
    ),
    paste0(
      "Primary survey-weighted logistic regression\n",
      "Sequential Models 0-4; n = ",
      comma(flow_lookup[["N_non_missing_exposure"]]), " to ",
      comma(flow_lookup[["N_complete_model4"]])
    ),
    "Secondary sleep outcomes and sleep-duration model\nFully adjusted survey-weighted analyses",
    "Dose-response, sensitivity, and alcohol-pattern-by-sex\ninteraction analyses",
    "Pairwise-complete Spearman correlations\nUnweighted exploratory analysis"
  )
)

workflow_arrows <- tibble(
  x = 0.5, xend = 0.5,
  y = seq(6.65, 1.65, by = -1),
  yend = seq(6.35, 1.35, by = -1)
)

workflow_colours <- c(
  Data = "#D9EAF7",
  Design = "#C7E9E3",
  Variables = "#F2E6C9",
  Primary = "#BFD7EA",
  Secondary = "#D9D2E9",
  Robustness = "#F4CCCC",
  Exploratory = "#E6E6E6"
)

figure1 <- ggplot() +
  geom_segment(
    data = workflow_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.75,
    colour = "#4D4D4D",
    arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "in"))
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x = x, y = y, label = label, fill = node_group),
    family = base_font,
    fontface = "bold",
    size = 3.55,
    lineheight = 1.0,
    label.padding = grid::unit(0.22, "lines"),
    label.r = grid::unit(0.14, "lines"),
    linewidth = 0.7,
    colour = "black"
  ) +
  scale_fill_manual(values = workflow_colours, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.55), clip = "off") +
  labs(
    title = "Methodological and analytical workflow",
    subtitle = "NHANES 2017-2018 alcohol-use and sleep analysis"
  ) +
  theme_void(base_family = base_font) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "bold", size = 11.5, hjust = 0.5,
                                 margin = margin(b = 10)),
    plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.9),
    plot.margin = margin(14, 14, 14, 14)
  )

save_figure(figure1, "Figure_1_methodological_workflow", width = 8.2, height = 10.0)

# ---- 7. Figure 2A: weighted prevalence of trouble sleeping -----------------
prevalence_plot_data <- prevalence %>%
  mutate(
    display_group = case_when(
      alcohol_category == "Lifetime abstainers" ~ "Lifetime abstainers",
      str_detect(alcohol_category, "without") ~
        "Ever drinkers without past-year binge drinking",
      TRUE ~ "Ever drinkers with past-year binge drinking"
    ),
    y = case_when(
      display_group == "Lifetime abstainers" ~ 3,
      str_detect(display_group, "without") ~ 2,
      TRUE ~ 1
    ),
    group_colour = case_when(
      display_group == "Lifetime abstainers" ~ colour_neutral,
      str_detect(display_group, "without") ~ colour_non_binge,
      TRUE ~ colour_binge
    ),
    result_label = paste0(
      formatC(prevalence, format = "f", digits = 1), "% (",
      formatC(CI_low, format = "f", digits = 1), "-",
      formatC(CI_high, format = "f", digits = 1), ")"
    )
  )

if (nrow(prevalence_plot_data) != 3L) {
  stop("Expected three alcohol categories in the prevalence table.")
}

overall_prevalence_p <- unique(prevalence_plot_data$p_overall)
if (length(overall_prevalence_p) != 1L) {
  stop("Expected one overall Rao-Scott p value in the prevalence table.")
}

pair_p_no_binge <- prevalence_pairwise %>%
  filter(str_detect(comparison, "without past-year binge vs lifetime")) %>%
  pull(p)
pair_p_binge <- prevalence_pairwise %>%
  filter(str_detect(comparison, "with past-year binge vs lifetime")) %>%
  pull(p)
pair_p_between <- prevalence_pairwise %>%
  filter(str_detect(comparison, "with vs without")) %>%
  pull(p)

if (length(pair_p_no_binge) != 1L || length(pair_p_binge) != 1L ||
    length(pair_p_between) != 1L) {
  stop("Could not identify all three prespecified pairwise prevalence comparisons.")
}

figure2a <- ggplot(prevalence_plot_data) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, y = y, yend = y, colour = display_group),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(x = CI_low, xend = CI_low, y = y - 0.08, yend = y + 0.08,
        colour = display_group),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(x = CI_high, xend = CI_high, y = y - 0.08, yend = y + 0.08,
        colour = display_group),
    linewidth = 0.9
  ) +
  geom_point(aes(x = prevalence, y = y, colour = display_group), size = 4.0) +
  geom_text(
    aes(x = CI_high + 0.8, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.7,
    colour = "black"
  ) +
  scale_colour_manual(values = setNames(
    prevalence_plot_data$group_colour,
    prevalence_plot_data$display_group
  ), guide = "none") +
  scale_x_continuous(
    limits = c(8, 52),
    breaks = c(10, 20, 30, 40, 50),
    labels = label_percent(scale = 1, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = c(3, 2, 1),
    labels = c(
      "Lifetime abstainers",
      "Ever drinkers without\npast-year binge drinking",
      "Ever drinkers with\npast-year binge drinking"
    ),
    limits = c(0.55, 3.45)
  ) +
  labs(
    title = "Survey-weighted prevalence of reported trouble sleeping",
    subtitle = paste0("Overall Rao-Scott ", p_text(overall_prevalence_p)),
    x = "Survey-weighted prevalence (95% CI)", y = NULL,
    caption = str_wrap(
      paste0(
        "Pairwise design-based tests: without binge vs abstainers, ",
        p_text(pair_p_no_binge), "; with binge vs abstainers, ",
        p_text(pair_p_binge), "; with vs without binge, ",
        p_text(pair_p_between), "."
      ),
      width = 110
    )
  ) +
  theme_publication(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold", lineheight = 0.95),
    panel.grid.major.x = element_line(colour = colour_grid, linewidth = 0.35)
  )

save_figure(figure2a, "Figure_2A_weighted_sleep_trouble_prevalence", width = 9.5, height = 5.4)

# ---- 8. Figure 2B: primary sequential-model forest plot --------------------
# A three-part layout reserves separate columns for row labels, confidence
# intervals, and numerical results. This prevents text/CI overlap.

model_levels <- paste("Model", 0:4)

primary_forest_data <- primary %>%
  filter(str_detect(term, "^alcohol_category")) %>%
  mutate(
    model_index = match(model, model_levels) - 1,
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~ "Without past-year binge",
      str_detect(term, "Binge drinker") ~ "With past-year binge",
      TRUE ~ NA_character_
    ),
    exposure_order = if_else(exposure == "Without past-year binge", 0, 1),
    y = 14 - 3 * model_index - exposure_order,
    result_text = paste0(
      number_2(OR), " (", number_2(CI_low), "-", number_2(CI_high), ")"
    ),
    p_value_text = case_when(
      p < 0.001 ~ "<0.001",
      TRUE ~ formatC(p, format = "f", digits = 3)
    )
  ) %>%
  left_join(
    primary_overall %>% select(model, p_overall_alcohol),
    by = "model"
  ) %>%
  mutate(
    overall_p_text = if_else(
      exposure_order == 0,
      case_when(
        p_overall_alcohol < 0.001 ~ "<0.001",
        TRUE ~ formatC(p_overall_alcohol, format = "f", digits = 3)
      ),
      ""
    )
  ) %>%
  arrange(model_index, exposure_order)

if (nrow(primary_forest_data) != 10L || any(is.na(primary_forest_data$model_index))) {
  stop("Expected ten alcohol-category estimates across primary Models 0-4.")
}

primary_model_rows <- primary_forest_data %>%
  group_by(model, model_index, n) %>%
  summarise(y = mean(y), .groups = "drop") %>%
  mutate(model_label = paste0(model, "\nn = ", comma(n)))

primary_bands <- tibble(
  model_index = 0:4,
  ymin = 12.55 - 3 * model_index,
  ymax = 14.45 - 3 * model_index,
  band_fill = if_else(model_index %% 2 == 0, colour_band, "white")
)

primary_y_limits <- c(0.35, 15.45)

primary_label_panel <- ggplot() +
  geom_rect(
    data = primary_bands,
    aes(xmin = 0, xmax = 3.0, ymin = ymin, ymax = ymax, fill = band_fill),
    colour = NA
  ) +
  scale_fill_identity() +
  geom_text(
    data = primary_model_rows,
    aes(x = 0.02, y = y, label = model_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.15,
    lineheight = 0.92
  ) +
  geom_point(
    data = primary_forest_data,
    aes(x = 1.12, y = y, colour = exposure),
    shape = 15, size = 2.8
  ) +
  geom_text(
    data = primary_forest_data,
    aes(x = 1.28, y = y, label = exposure),
    hjust = 0, family = base_font, fontface = "bold", size = 3.0
  ) +
  annotate("text", x = 0.02, y = 15.12, label = "Model",
           hjust = 0, family = base_font, fontface = "bold", size = 3.25) +
  annotate("text", x = 1.12, y = 15.12, label = "Alcohol-use group",
           hjust = 0, family = base_font, fontface = "bold", size = 3.25) +
  scale_colour_manual(values = c(
    "Without past-year binge" = colour_non_binge,
    "With past-year binge" = colour_binge
  ), guide = "none") +
  scale_x_continuous(limits = c(0, 3.0), expand = c(0, 0)) +
  scale_y_continuous(limits = primary_y_limits, expand = c(0, 0)) +
  theme_void(base_family = base_font) +
  theme(plot.margin = margin(5, 2, 5, 5))

primary_ci_panel <- ggplot(primary_forest_data, aes(y = y)) +
  geom_rect(
    data = primary_bands,
    aes(xmin = 0.48, xmax = 4.5, ymin = ymin, ymax = ymax, fill = band_fill),
    inherit.aes = FALSE, colour = NA
  ) +
  scale_fill_identity() +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555", linewidth = 0.65) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, y = y, yend = y, colour = exposure),
    linewidth = 0.8
  ) +
  geom_segment(
    aes(x = CI_low, xend = CI_low, y = y - 0.08, yend = y + 0.08,
        colour = exposure), linewidth = 0.8
  ) +
  geom_segment(
    aes(x = CI_high, xend = CI_high, y = y - 0.08, yend = y + 0.08,
        colour = exposure), linewidth = 0.8
  ) +
  geom_point(aes(x = OR, colour = exposure), size = 3.1) +
  scale_colour_manual(values = c(
    "Without past-year binge" = colour_non_binge,
    "With past-year binge" = colour_binge
  ), guide = "none") +
  scale_x_log10(
    limits = c(0.48, 4.5), breaks = c(0.5, 1, 2, 4),
    labels = c("0.5", "1", "2", "4"), expand = c(0, 0)
  ) +
  scale_y_continuous(limits = primary_y_limits, expand = c(0, 0)) +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_classic(base_size = 10.5, base_family = base_font) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = colour_grid, linewidth = 0.35),
    plot.margin = margin(5, 3, 5, 3)
  )

primary_numeric_panel <- ggplot() +
  geom_rect(
    data = primary_bands,
    aes(xmin = 0, xmax = 4.7, ymin = ymin, ymax = ymax, fill = band_fill),
    colour = NA
  ) +
  scale_fill_identity() +
  geom_text(
    data = primary_forest_data,
    aes(x = 0.02, y = y, label = result_text),
    hjust = 0, family = base_font, fontface = "bold", size = 3.0
  ) +
  geom_text(
    data = primary_forest_data,
    aes(x = 2.82, y = y, label = p_value_text),
    hjust = 0, family = base_font, fontface = "bold", size = 3.0
  ) +
  geom_text(
    data = primary_forest_data,
    aes(x = 3.78, y = y, label = overall_p_text),
    hjust = 0, family = base_font, fontface = "bold", size = 3.0
  ) +
  annotate("text", x = 0.02, y = 15.12, label = "OR (95% CI)",
           hjust = 0, family = base_font, fontface = "bold", size = 3.2) +
  annotate("text", x = 2.82, y = 15.12, label = "p value",
           hjust = 0, family = base_font, fontface = "bold", size = 3.2) +
  annotate("text", x = 3.78, y = 15.12, label = "Overall p",
           hjust = 0, family = base_font, fontface = "bold", size = 3.2) +
  scale_x_continuous(limits = c(0, 4.7), expand = c(0, 0)) +
  scale_y_continuous(limits = primary_y_limits, expand = c(0, 0)) +
  theme_void(base_family = base_font) +
  theme(plot.margin = margin(5, 5, 5, 2))

figure2b <- (primary_label_panel | primary_ci_panel | primary_numeric_panel) +
  plot_layout(widths = c(1.42, 1.30, 1.78)) +
  plot_annotation(
    title = "Alcohol-use pattern and reported trouble sleeping",
    subtitle = paste0(
      "Sequential survey-weighted logistic regression; ",
      "lifetime abstainers are the reference group"
    ),
    theme = theme(
      text = element_text(family = base_font, colour = "black"),
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 3)),
      plot.subtitle = element_text(face = "bold", size = 11,
                                   margin = margin(b = 8)),
      plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.9),
      plot.margin = margin(10, 10, 8, 10)
    )
  )

save_figure(figure2b, "Figure_2B_primary_sequential_models", width = 12.0, height = 7.4)

# ---- 9. Figure 3A: secondary binary sleep outcomes -------------------------
secondary_plot_data <- secondary %>%
  filter(str_detect(term, "^alcohol_category")) %>%
  mutate(
    outcome_label = recode(
      outcome,
      insufficient_sleep = "Insufficient sleep",
      daytime_sleepy = "Daytime sleepiness",
      snoring_frequent = "Frequent snoring",
      breathing_pauses = "Breathing pauses during sleep"
    ),
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~ "Without past-year binge drinking",
      str_detect(term, "Binge drinker") ~ "With past-year binge drinking",
      TRUE ~ NA_character_
    )
  )

expected_secondary_rows <- 8L
if (nrow(secondary_plot_data) != expected_secondary_rows) {
  stop("Expected 8 alcohol-category rows in secondary_binary_models.csv; found ",
       nrow(secondary_plot_data), ".")
}

outcome_order <- c(
  "Insufficient sleep",
  "Daytime sleepiness",
  "Frequent snoring",
  "Breathing pauses during sleep"
)

secondary_plot_data <- secondary_plot_data %>%
  mutate(
    outcome_index = match(outcome_label, outcome_order),
    outcome_center = rev(seq_along(outcome_order))[outcome_index],
    y = outcome_center + if_else(exposure == "Without past-year binge drinking", 0.15, -0.15),
    result_label = or_result_text(OR, CI_low, CI_high, p)
  )

p3a <- ggplot(secondary_plot_data) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35", linewidth = 0.7)
p3a <- add_horizontal_ci(p3a, "exposure", cap_height = 0.07) +
  geom_point(aes(x = OR, y = y, colour = exposure), size = 3.2) +
  geom_text(
    aes(x = 3.55, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.45,
    colour = "black"
  ) +
  scale_colour_manual(values = c(
    "Without past-year binge drinking" = colour_non_binge,
    "With past-year binge drinking" = colour_binge
  )) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
    labels = label_number(accuracy = 0.01),
    # The unused right-hand portion is reserved for the numerical labels.
    limits = c(0.25, 16)
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(outcome_order)),
    labels = outcome_order,
    expand = expansion(add = c(0.45, 0.45))
  ) +
  labs(
    title = "Secondary binary sleep outcomes",
    subtitle = "Fully adjusted survey-weighted logistic regression",
    x = "Odds ratio (log scale)", y = NULL,
    caption = paste0(
      "Lifetime abstainers are the reference group. Points are odds ratios; lines are design-based 95% CIs. ",
      "Exact p values are printed beside the estimates."
    )
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(p3a, "Figure_3A_secondary_binary_outcomes", width = 11.5, height = 6.2)

# ---- 10. Figure 3B: sleep duration by alcohol-use pattern ------------------
sleep_plot_data <- sleep_hours %>%
  filter(str_detect(term, "^alcohol_category")) %>%
  mutate(
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~ "Without past-year binge drinking",
      str_detect(term, "Binge drinker") ~ "With past-year binge drinking",
      TRUE ~ NA_character_
    ),
    y = if_else(exposure == "Without past-year binge drinking", 2, 1),
    result_label = beta_result_text(beta, CI_low, CI_high, p, digits = 2)
  )

if (nrow(sleep_plot_data) != 2L) {
  stop("Expected 2 alcohol-category rows in sleep_hours_model.csv; found ",
       nrow(sleep_plot_data), ".")
}

p3b <- ggplot(sleep_plot_data) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey35", linewidth = 0.7)
p3b <- add_horizontal_ci(p3b, "exposure", cap_height = 0.07) +
  geom_point(aes(x = beta, y = y, colour = exposure), size = 3.4) +
  geom_text(
    aes(x = 0.43, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.6,
    colour = "black"
  ) +
  scale_colour_manual(values = c(
    "Without past-year binge drinking" = colour_non_binge,
    "With past-year binge drinking" = colour_binge
  ), guide = "none") +
  scale_x_continuous(
    breaks = seq(-0.5, 1.0, by = 0.25),
    limits = c(-0.5, 1.0),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    breaks = c(2, 1),
    labels = c("Without past-year binge drinking", "With past-year binge drinking"),
    expand = expansion(add = c(0.45, 0.45))
  ) +
  labs(
    title = "Alcohol-use pattern and sleep duration",
    subtitle = "Fully adjusted survey-weighted linear regression; n = 3,099",
    x = "Adjusted difference in sleep duration, hours", y = NULL,
    caption = paste0(
      "Lifetime abstainers are the reference group. Points are beta coefficients; lines are design-based 95% CIs. ",
      "Positive values indicate longer reported sleep duration."
    )
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(p3b, "Figure_3B_sleep_duration_by_alcohol_pattern", width = 11.5, height = 4.3)

# ---- 11. Figure 3C: sleep duration per additional drink --------------------
dose_sleep <- dose_response %>%
  filter(outcome == "sleep_hours_by_drinks_per_day", term == "drinks_per_day") %>%
  mutate(
    y = 1,
    result_label = beta_result_text(effect, CI_low, CI_high, p, digits = 3)
  )

if (nrow(dose_sleep) != 1L) {
  stop("Expected one drinks_per_day row for sleep duration; found ", nrow(dose_sleep), ".")
}

p3c <- ggplot(dose_sleep) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey35", linewidth = 0.7)
p3c <- add_horizontal_ci(p3c, cap_height = 0.08) +
  geom_point(aes(x = effect, y = y), size = 3.6, colour = colour_accent) +
  geom_text(
    aes(x = 0.038, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.7,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = seq(-0.06, 0.12, by = 0.03),
    limits = c(-0.06, 0.12),
    labels = label_number(accuracy = 0.001)
  ) +
  scale_y_continuous(
    breaks = 1, labels = "Sleep duration",
    limits = c(0.55, 1.45)
  ) +
  labs(
    title = "Dose-response association with sleep duration",
    subtitle = "Fully adjusted survey-weighted model among ever drinkers; n = 2,750",
    x = "Change in sleep duration per additional drink, hours", y = NULL,
    caption = "The point is the beta coefficient and the horizontal line is the design-based 95% CI."
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(p3c, "Figure_3C_dose_response_sleep_duration", width = 10.5, height = 3.5)

# ---- 12. Figure 3D: trouble sleeping per additional drink ------------------
dose_trouble <- dose_response %>%
  filter(outcome == "sleep_trouble_by_drinks_per_day", term == "drinks_per_day") %>%
  mutate(
    y = 1,
    result_label = or_result_text(effect, CI_low, CI_high, p)
  )

if (nrow(dose_trouble) != 1L) {
  stop("Expected one drinks_per_day row for trouble sleeping; found ", nrow(dose_trouble), ".")
}

p3d <- ggplot(dose_trouble) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35", linewidth = 0.7)
p3d <- add_horizontal_ci(p3d, cap_height = 0.08) +
  geom_point(aes(x = effect, y = y), size = 3.6, colour = colour_accent) +
  geom_text(
    aes(x = 1.11, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.7,
    colour = "black"
  ) +
  scale_x_log10(
    breaks = c(0.75, 0.85, 1.00, 1.15, 1.30, 1.50),
    limits = c(0.75, 1.50),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    breaks = 1, labels = "Reported trouble sleeping",
    limits = c(0.55, 1.45)
  ) +
  labs(
    title = "Dose-response association with reported trouble sleeping",
    subtitle = "Fully adjusted survey-weighted model among ever drinkers; n = 2,758",
    x = "Odds ratio per additional drink (log scale)", y = NULL,
    caption = "The point is the odds ratio and the horizontal line is the design-based 95% CI."
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(p3d, "Figure_3D_dose_response_sleep_trouble", width = 10.5, height = 3.5)

# ---- 13. Figure 4A: sensitivity analyses -----------------------------------
sensitivity_plot_data <- sensitivity %>%
  filter(
    term %in% c(
      "alcohol_userYes",
      "alcohol_riskLow/Moderate",
      "alcohol_riskHigh",
      "alcohol_categoryDrinkers (non-binge)",
      "alcohol_categoryBinge drinker"
    )
  ) %>%
  mutate(
    display_label = case_when(
      term == "alcohol_userYes" ~ "Ever drinker vs lifetime abstainer",
      term == "alcohol_riskLow/Moderate" ~ "Low/moderate risk vs no alcohol use",
      term == "alcohol_riskHigh" ~ "High risk vs no alcohol use",
      outcome == "sleep_trouble_PHQ_full" & str_detect(term, "non-binge") ~
        "Full PHQ: without binge vs abstainer",
      outcome == "sleep_trouble_PHQ_full" & str_detect(term, "Binge drinker") ~
        "Full PHQ: with binge vs abstainer",
      TRUE ~ term
    ),
    plot_group = case_when(
      term == "alcohol_userYes" ~ "Ever alcohol use",
      str_detect(term, "Low/Moderate|non-binge") ~ "Lower-risk/without binge",
      str_detect(term, "High|Binge drinker") ~ "Higher-risk/with binge",
      TRUE ~ "Other"
    )
  )

sensitivity_order <- c(
  "Ever drinker vs lifetime abstainer",
  "Low/moderate risk vs no alcohol use",
  "High risk vs no alcohol use",
  "Full PHQ: without binge vs abstainer",
  "Full PHQ: with binge vs abstainer"
)

sensitivity_plot_data <- sensitivity_plot_data %>%
  mutate(
    y = rev(seq_along(sensitivity_order))[match(display_label, sensitivity_order)],
    result_label = or_result_text(OR, CI_low, CI_high, p)
  )

if (nrow(sensitivity_plot_data) != 5L || any(is.na(sensitivity_plot_data$y))) {
  stop("Expected 5 prespecified sensitivity estimates in sensitivity_models.csv.")
}

p4a <- ggplot(sensitivity_plot_data) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35", linewidth = 0.7)
p4a <- add_horizontal_ci(p4a, "plot_group", cap_height = 0.07) +
  geom_point(aes(x = OR, y = y, colour = plot_group), size = 3.3) +
  geom_text(
    aes(x = 3.65, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.5,
    colour = "black"
  ) +
  scale_colour_manual(values = c(
    "Ever alcohol use" = colour_neutral,
    "Lower-risk/without binge" = colour_non_binge,
    "Higher-risk/with binge" = colour_binge
  )) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4, 8, 16),
    labels = label_number(accuracy = 0.01),
    # The unused right-hand portion is reserved for the numerical labels.
    limits = c(0.5, 16)
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(sensitivity_order)),
    labels = sensitivity_order,
    expand = expansion(add = c(0.45, 0.45))
  ) +
  labs(
    title = "Sensitivity analyses for reported trouble sleeping",
    subtitle = "Fully adjusted survey-weighted logistic regression",
    x = "Odds ratio (log scale)", y = NULL,
    caption = paste0(
      "Reference groups are stated in the row labels. Points are odds ratios; lines are design-based 95% CIs. ",
      "Full PHQ models use the PHQ score including the sleep item."
    )
  ) +
  theme_publication()

save_figure(p4a, "Figure_4A_sensitivity_analyses", width = 12.0, height = 6.5)

# ---- 14. Figure 4B: alcohol-pattern-by-sex interaction ---------------------
# Only the interaction terms are plotted. These are ratios of odds ratios,
# not sex-specific marginal odds ratios.
interaction_plot_data <- interaction %>%
  filter(str_detect(term, ":genderMale$")) %>%
  mutate(
    display_label = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~
        "Without-binge pattern x male sex",
      str_detect(term, "Binge drinker") ~
        "With-binge pattern x male sex",
      TRUE ~ term
    ),
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~ "Without past-year binge drinking",
      str_detect(term, "Binge drinker") ~ "With past-year binge drinking",
      TRUE ~ NA_character_
    ),
    y = if_else(str_detect(term, "non-binge"), 2, 1),
    result_label = paste0(
      "Ratio of ORs ", number_2(OR), " (", number_2(CI_low), "-", number_2(CI_high), "); ",
      p_text(p)
    )
  )

if (nrow(interaction_plot_data) != 2L) {
  stop("Expected 2 alcohol-pattern-by-sex interaction terms; found ",
       nrow(interaction_plot_data), ".")
}

joint_p <- interaction_joint$p[[1]]
joint_n <- interaction_joint$n[[1]]
joint_label <- paste0("Joint alcohol-pattern-by-sex interaction: ", p_text(joint_p))

p4b <- ggplot(interaction_plot_data) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35", linewidth = 0.7)
p4b <- add_horizontal_ci(p4b, "exposure", cap_height = 0.07) +
  geom_point(aes(x = OR, y = y, colour = exposure), size = 3.5) +
  geom_text(
    aes(x = 3.6, y = y, label = result_label),
    hjust = 0, family = base_font, fontface = "bold", size = 3.55,
    colour = "black"
  ) +
  annotate(
    "text", x = 0.22, y = 2.62, label = joint_label,
    hjust = 0, family = base_font, fontface = "bold", size = 4.0
  ) +
  scale_colour_manual(values = c(
    "Without past-year binge drinking" = colour_non_binge,
    "With past-year binge drinking" = colour_binge
  ), guide = "none") +
  scale_x_log10(
    breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
    labels = label_number(accuracy = 0.01),
    # The unused right-hand portion is reserved for the numerical labels.
    limits = c(0.125, 16)
  ) +
  scale_y_continuous(
    breaks = c(2, 1),
    labels = c(
      "Without-binge pattern x male sex",
      "With-binge pattern x male sex"
    ),
    limits = c(0.55, 2.8)
  ) +
  labs(
    title = "Alcohol-use-pattern-by-sex interaction",
    subtitle = paste0("Fully adjusted survey-weighted logistic regression; n = ",
                      format(joint_n, big.mark = ",")),
    x = "Ratio of odds ratios (log scale)", y = NULL,
    caption = paste0(
      "Points are interaction-term ratios of odds ratios; lines are design-based 95% CIs. ",
      "Values are not sex-specific marginal odds ratios."
    )
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(p4b, "Figure_4B_sex_interaction", width = 12.0, height = 5.2)

# ---- 15. Figure 5: Spearman correlation heatmap ----------------------------
# Only one triangle is printed because the correlation matrix is symmetric.
# Blue represents negative correlations, white represents values near zero,
# and vermillion represents positive correlations. This palette is balanced,
# print-readable, and distinguishable for common colour-vision deficiencies.

correlation_variables <- c(
  "drinks_per_day",
  "drinking_frequency_rank",
  "sleep_hours",
  "age_years",
  "bmi",
  "waist_circumference",
  "poverty_ratio",
  "depression_score_no_sleep"
)

correlation_labels <- c(
  drinks_per_day = "Drinks per\ndrinking day",
  drinking_frequency_rank = "Drinking\nfrequency",
  sleep_hours = "Sleep\nduration",
  age_years = "Age",
  bmi = "BMI",
  waist_circumference = "Waist\ncircumference",
  poverty_ratio = "Poverty\nratio",
  depression_score_no_sleep = "PHQ score\n(sleep item excluded)"
)

correlation_plot_data <- spearman_long %>%
  mutate(
    row_index = match(variable_1, correlation_variables),
    column_index = match(variable_2, correlation_variables)
  ) %>%
  filter(!is.na(row_index), !is.na(column_index), row_index > column_index) %>%
  mutate(
    stars = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    cell_label = paste0(formatC(rho, format = "f", digits = 2), stars),
    text_colour = if_else(abs(rho) >= 0.65, "white", "black"),
    row_label = factor(
      unname(correlation_labels[variable_1]),
      levels = unname(correlation_labels[correlation_variables])
    ),
    column_label = factor(
      unname(correlation_labels[variable_2]),
      levels = unname(correlation_labels[correlation_variables])
    )
  )

expected_correlation_cells <- length(correlation_variables) *
  (length(correlation_variables) - 1) / 2
if (nrow(correlation_plot_data) != expected_correlation_cells) {
  stop(
    "Expected ", expected_correlation_cells,
    " unique lower-triangle correlations; found ",
    nrow(correlation_plot_data), "."
  )
}

figure5 <- ggplot(
  correlation_plot_data,
  aes(x = column_label, y = row_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(
    aes(label = cell_label, colour = text_colour),
    family = base_font,
    fontface = "bold",
    size = 4.0,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#FFFFF7",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Spearman rho"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(3.2, "in"),
      barheight = grid::unit(0.20, "in")
    )
  ) +
  labs(
    title = "Spearman correlations among alcohol, sleep, and covariates",
    subtitle = "Pairwise-complete, unweighted exploratory analysis",
    x = NULL, y = NULL,
    caption = str_wrap(
      paste0(
        "Cells show Spearman correlation coefficients. ",
        "*p < 0.05; **p < 0.01; ***p < 0.001. ",
        "Cells without asterisks indicate p >= 0.05; ",
        "p values were not adjusted for multiple comparisons."
      ),
      width = 115
    )
  ) +
  theme_minimal(base_size = 11.5, base_family = base_font) +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 3)),
    plot.subtitle = element_text(face = "bold", size = 11.5,
                                 margin = margin(b = 9)),
    plot.caption = element_text(size = 9.5, hjust = 0, margin = margin(t = 9)),
    axis.text.x = element_text(
      angle = 42, hjust = 1, vjust = 1,
      face = "bold", colour = "black", size = 10,
      lineheight = 0.95
    ),
    axis.text.y = element_text(
      face = "bold", colour = "black", size = 10,
      lineheight = 0.95
    ),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10.5),
    legend.text = element_text(size = 9.5),
    plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.9),
    plot.margin = margin(14, 16, 14, 16)
  )

save_figure(figure5, "Figure_5_spearman_correlation", width = 10.5, height = 9.2)

# ---- 16. Completion message -------------------------------------------------
message("\nAll separate publication figures were created successfully in:")
message(figures_dir)
message("\nFiles created:")
message(" - Figure_1_methodological_workflow")
message(" - Figure_2A_weighted_sleep_trouble_prevalence")
message(" - Figure_2B_primary_sequential_models")
message(" - Figure_3A_secondary_binary_outcomes")
message(" - Figure_3B_sleep_duration_by_alcohol_pattern")
message(" - Figure_3C_dose_response_sleep_duration")
message(" - Figure_3D_dose_response_sleep_trouble")
message(" - Figure_4A_sensitivity_analyses")
message(" - Figure_4B_sex_interaction")
message(" - Figure_5_spearman_correlation")

