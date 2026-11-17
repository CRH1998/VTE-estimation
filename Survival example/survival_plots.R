# Replace the summary_results calculation and plotting code with this script.
# Run after results_df has been created; the simulations do not need to be rerun.
# Requires ggplot2 >= 3.4.0, dplyr, tidyr, and scales.

library(dplyr)
library(tidyr)
library(ggplot2)

# Shared theme: figure titles can be supplied in the article's captions.
theme_article <- function(base_size = 12, base_family = "serif") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        face = "bold", size = base_size + 1, hjust = 0,
        margin = margin(b = 10)
      ),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(1.1, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

# Keep the original beta values for grouping; round only their displayed labels.
# Named mappings ensure the same beta has the same colour in every figure.
beta_levels <- sort(unique(summary_results$beta))
article_palette <- c(
  "#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9"
)
article_shapes <- c(16, 17, 15, 18, 3, 7, 8)
stopifnot(length(beta_levels) <= length(article_palette))

beta_colours <- setNames(
  article_palette[seq_along(beta_levels)], as.character(beta_levels)
)
beta_shapes <- setNames(
  article_shapes[seq_along(beta_levels)], as.character(beta_levels)
)
beta_labels <- function(x) sprintf("%.2f", as.numeric(as.character(x)))

beta_colour_scale <- function() {
  scale_colour_manual(
    name = expression(beta), values = beta_colours, labels = beta_labels
  )
}

sample_size_scale <- function(data) {
  scale_x_continuous(
    breaks = sort(unique(data$n)),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.03, 0.03))
  )
}

test_scales <- function() {
  list(
    scale_linetype_manual(
      name = "Test",
      values = c("Log-rank" = "solid", "Two-Step" = "dashed")
    ),
    scale_shape_manual(
      name = "Test", values = c("Log-rank" = 16, "Two-Step" = 17)
    )
  )
}

# 1. Empirical level: only beta = 0, so a beta colour legend is unnecessary.
ggplot_empirical_levels <- function(summary_results, alpha = 0.05) {
  plot_data <- summary_results %>%
    filter(beta == 0) %>%
    pivot_longer(
      cols = c(LogRankLevel, TwoStepLevel),
      names_to = "test", values_to = "empirical_level"
    ) %>%
    mutate(test = factor(
      test, levels = c("LogRankLevel", "TwoStepLevel"),
      labels = c("Log-rank", "Two-Step")
    ))

  ggplot(plot_data, aes(
    x = n, y = empirical_level, linetype = test, shape = test
  )) +
    geom_hline(yintercept = alpha, linetype = "dotted",
               colour = "grey50", linewidth = 0.45) +
    geom_line(linewidth = 0.65, colour = "black") +
    geom_point(size = 2.2, colour = "black") +
    test_scales() +
    sample_size_scale(plot_data) +
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.08))
    ) +
    expand_limits(y = 2 * alpha) +
    labs(
      title = "Empirical type I error",
      x = expression("Sample size (" * n * ")"), y = "Type I error"
    ) +
    theme_article()
}

# 2. Empirical power: colour identifies beta; line and symbol identify the test.
ggplot_empirical_power <- function(summary_results) {
  plot_data <- summary_results %>%
    filter(beta > 0) %>%
    mutate(beta_plot = factor(beta, levels = beta_levels)) %>%
    pivot_longer(
      cols = c(LogRankPower, TwoStepPower),
      names_to = "test", values_to = "empirical_power"
    ) %>%
    mutate(test = factor(
      test, levels = c("LogRankPower", "TwoStepPower"),
      labels = c("Log-rank", "Two-Step")
    ))

  ggplot(plot_data, aes(
    x = n, y = empirical_power, colour = beta_plot,
    linetype = test, shape = test, group = interaction(beta_plot, test)
  )) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2) +
    beta_colour_scale() +
    test_scales() +
    sample_size_scale(plot_data) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = "Empirical power",
      x = expression("Sample size (" * n * ")"), y = "Power"
    ) +
    guides(
      colour = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(linetype = "solid", shape = 16)
      ),
      linetype = guide_legend(order = 2, nrow = 1),
      shape = guide_legend(order = 2, nrow = 1)
    ) +
    theme_article()
}

# 3. RMSE: symbols supplement colour so curves are easier to distinguish in print.
ggplot_empirical_RMSE <- function(summary_results) {
  plot_data <- summary_results %>%
    mutate(beta_plot = factor(beta, levels = beta_levels))

  ggplot(plot_data, aes(
    x = n, y = RMSE, colour = beta_plot, shape = beta_plot, group = beta_plot
  )) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.2) +
    beta_colour_scale() +
    scale_shape_manual(
      name = expression(beta), values = beta_shapes, labels = beta_labels
    ) +
    sample_size_scale(plot_data) +
    scale_y_continuous(
      limits = c(0, NA), labels = scales::label_number(),
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      title = "RMSE of the reduced two-step estimator",
      x = expression("Sample size (" * n * ")"), y = "RMSE"
    ) +
    theme_article()
}


# 4. Build one shared legend independently of the plotted power data.
# Every beta is represented in this legend-only data, including beta = 0.
# These rows are used only to draw legend keys; they never enter any data panel.
make_shared_legend <- function() {
  legend_data <- expand.grid(
    key_x = c(0, 1),
    beta_plot = factor(beta_levels, levels = beta_levels),
    test = factor(c("Log-rank", "Two-Step"),
                  levels = c("Log-rank", "Two-Step")),
    KEEP.OUT.ATTRS = FALSE
  )
  
  legend_plot <- ggplot(legend_data, aes(
    x = key_x, y = 0, colour = beta_plot, linetype = test, shape = test,
    group = interaction(beta_plot, test)
  )) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2) +
    beta_colour_scale() +
    test_scales() +
    guides(
      colour = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(shape = NA, linetype = "solid", linewidth = 0.65)
      ),
      linetype = guide_legend(
        order = 2, nrow = 1,
        override.aes = list(colour = "black", linewidth = 0.65, size = 2)
      ),
      shape = guide_legend(order = 2, nrow = 1)
    ) +
    theme_article() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(0.7, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
  # Find the populated legend box; its name differs across ggplot2 versions.
  legend_table <- ggplotGrob(legend_plot)
  legend_index <- which(
    grepl("^guide-box", legend_table$layout$name) &
      vapply(legend_table$grobs, inherits, logical(1), what = "gtable")
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the shared legend from the legend-only plot.")
  }
  legend_table$grobs[[legend_index]]
}

# 5. Combined figure: all three panel legends are hidden, and the complete
# shared legend is placed in its own row beneath the panels.
combine_article_plots <- function(p_levels, p_power, p_rmse, alpha = 0.05) {
  level_panel <- p_levels +
    labs(title = "Type I error") +
    guides(colour = "none", linetype = "none", shape = "none")
  
  power_panel <- p_power +
    labs(title = "Power") +
    guides(colour = "none", linetype = "none", shape = "none")
  
  rmse_panel <- p_rmse +
    labs(title = "RMSE of reduced two-step estimator") +
    guides(colour = "none", shape = "none", linetype = "none")
  
  shared_legend <- patchwork::wrap_elements(
    full = make_shared_legend(), ignore_tag = TRUE
  )
  
  combined <- patchwork::wrap_plots(
    A = level_panel, B = power_panel, C = rmse_panel, D = shared_legend,
    design = "ABC\nDDD", widths = c(1, 1, 1), heights = c(1, 0.16)
  ) & theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold"),
    plot.tag = element_text(size = 11, face = "bold", family = "serif"),
    plot.margin = margin(6, 6, 6, 6)
  )
  
  combined + patchwork::plot_annotation(
    theme = theme(
      text = element_text(family = "serif"),
      plot.title = element_text(size = 13, face = "bold", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 6))
    )
  )
}





################################################################################
#                           Non-proportional hazards plots                     #
################################################################################

theme_article_nonprop <- function(base_size = 12, base_family = "serif") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        size = base_size + 1, face = "bold", hjust = 0, margin = margin(b = 10)
      ),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(0.7, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# Preserve exact parameter values for grouping; round labels only.
# Names are distinct from those in article_survival_plots.R so both scripts
# can be used in the same R session.
np_beta_levels <- sort(unique(summary_results_nonprop$beta))
np_delta_levels <- sort(unique(summary_results_nonprop$delta), decreasing = TRUE)
np_test_levels <- c("Log-rank", "Two-Step")
np_palette <- c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#E69F00", "#56B4E9")
np_shapes <- c(16, 17, 15, 18, 3, 7, 8)
stopifnot(length(np_delta_levels) <= length(np_palette),
          length(np_beta_levels) <= length(np_shapes))

np_delta_colours <- setNames(np_palette[seq_along(np_delta_levels)],
                             as.character(np_delta_levels))
np_beta_shapes <- setNames(np_shapes[seq_along(np_beta_levels)],
                           as.character(np_beta_levels))
np_test_linetypes <- c("Log-rank" = "solid", "Two-Step" = "dashed")
np_coefficient_labels <- function(x) sprintf("%.2f", as.numeric(as.character(x)))

nonprop_scales <- function() {
  list(
    scale_colour_manual(name = expression(delta), values = np_delta_colours,
                        labels = np_coefficient_labels),
    scale_shape_manual(name = expression(beta), values = np_beta_shapes,
                       labels = np_coefficient_labels),
    scale_linetype_manual(name = "Test", values = np_test_linetypes)
  )
}

nonprop_guides <- function() {
  guides(
    colour = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(shape = NA, linetype = "solid", linewidth = 0.65)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(colour = "black", linetype = "blank", size = 2.2)
    ),
    linetype = guide_legend(
      order = 3, nrow = 1,
      override.aes = list(colour = "black", shape = NA, linewidth = 0.65)
    )
  )
}

nonprop_sample_size_scale <- function(data) {
  scale_x_continuous(
    breaks = sort(unique(data$n)),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.03, 0.03))
  )
}

nonprop_plot_data <- function(summary_results_nonprop) {
  summary_results_nonprop %>%
    pivot_longer(
      cols = c(LogRankReject, TwoStepReject),
      names_to = "test", values_to = "rejection_rate"
    ) %>%
    mutate(
      test = factor(test, levels = c("LogRankReject", "TwoStepReject"),
                    labels = np_test_levels),
      beta_plot = factor(beta, levels = np_beta_levels),
      delta_plot = factor(delta, levels = np_delta_levels)
    )
}

ggplot_empirical_levels_nonprop <- function(summary_results_nonprop,
                                            alpha = 0.05) {
  plot_data <- nonprop_plot_data(summary_results_nonprop) %>% filter(is_null)
  if (nrow(plot_data) == 0L) stop("The level plot requires beta = delta = 0.")
  
  ggplot(plot_data, aes(
    x = n, y = rejection_rate, colour = delta_plot,
    shape = beta_plot, linetype = test, group = test
  )) +
    geom_hline(yintercept = alpha, linetype = "dotted",
               colour = "grey50", linewidth = 0.45) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.2) +
    nonprop_scales() +
    nonprop_guides() +
    nonprop_sample_size_scale(plot_data) +
    scale_y_continuous(
      limits = c(0, NA), labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.08))
    ) +
    expand_limits(y = 2 * alpha) +
    labs(
      title = "Empirical type I error",
      x = expression("Sample size (" * n * ")"), y = "Type I error"
    ) +
    theme_article_nonprop()
}

ggplot_empirical_power_nonprop <- function(summary_results_nonprop) {
  plot_data <- nonprop_plot_data(summary_results_nonprop) %>% filter(!is_null)
  if (nrow(plot_data) == 0L) stop("The power plot requires a non-null setting.")
  
  ggplot(plot_data, aes(
    x = n, y = rejection_rate, colour = delta_plot,
    shape = beta_plot, linetype = test,
    group = interaction(beta_plot, delta_plot, test)
  )) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.1) +
    nonprop_scales() +
    nonprop_guides() +
    nonprop_sample_size_scale(plot_data) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = "Empirical power",
      x = expression("Sample size (" * n * ")"), y = "Power"
    ) +
    theme_article_nonprop()
}

# Build the complete shared legend independently of either panel's data.
# The rows below draw only legend keys; they never enter a data panel.
make_shared_legend_nonprop <- function() {
  legend_data <- expand.grid(
    key_x = c(0, 1),
    beta_plot = factor(np_beta_levels, levels = np_beta_levels),
    delta_plot = factor(np_delta_levels, levels = np_delta_levels),
    test = factor(np_test_levels, levels = np_test_levels),
    KEEP.OUT.ATTRS = FALSE
  )
  
  legend_plot <- ggplot(legend_data, aes(
    x = key_x, y = 0, colour = delta_plot, shape = beta_plot, linetype = test,
    group = interaction(beta_plot, delta_plot, test)
  )) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.1) +
    nonprop_scales() +
    nonprop_guides() +
    theme_article_nonprop() +
    theme(legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0))
  
  legend_table <- ggplotGrob(legend_plot)
  legend_index <- which(
    grepl("^guide-box", legend_table$layout$name) &
      vapply(legend_table$grobs, inherits, logical(1), what = "gtable")
  )
  if (length(legend_index) != 1L) stop("Could not extract the shared legend.")
  legend_table$grobs[[legend_index]]
}

combine_article_plots_nonprop <- function(p_levels_nonprop, p_power_nonprop,
                                          alpha = 0.05) {
  no_guides <- guides(colour = "none", shape = "none", linetype = "none")
  level_panel <- p_levels_nonprop + no_guides + labs(title = "Type I error")
  power_panel <- p_power_nonprop + no_guides + labs(title = "Power")
  shared_legend <- patchwork::wrap_elements(
    full = make_shared_legend_nonprop(), ignore_tag = TRUE
  )
  
  # The power panel has seven alternative settings (14 curves for two tests)
  # with the supplied grid, so it receives more horizontal space.
  combined <- patchwork::wrap_plots(
    A = level_panel, B = power_panel, C = shared_legend,
    design = "AB\nCC", widths = c(1, 1.7), heights = c(1, 0.18)
  ) & theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold"),
    plot.tag = element_text(size = 11, face = "bold", family = "serif"),
    plot.margin = margin(6, 6, 6, 6)
  )
  
  combined + patchwork::plot_annotation(
    theme = theme(
      text = element_text(family = "serif"),
      plot.title = element_text(size = 13, face = "bold", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 6))
    )
  )
}