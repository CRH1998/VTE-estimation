# Empirical level and power for the PM heterogeneity simulation
# Source immediately after assigning the column names of results_df in section 1.4.
# Uses results_df, settings, setting_index, and analysis_type from that section.
# Creates pm_summary plus pm_plots$level, pm_plots$power, and pm_plots$combined.
# A single sigma in the results is selected automatically. For older runs with
# several sigma values, set pm_sigma to the desired value before sourcing.

library(dplyr)
library(ggplot2)

pm_rejection_summary <- function(results_df, settings, setting_index, analysis_type,
                                 alpha = 0.05) {
  analyses <- unique(as.character(analysis_type))
  outcome <- c(SB = "ci_low_multifold", TwoStep = "CRHP", PM = "PM_test_reject",
               Dukes = "DukesP")
  # stopifnot(length(analyses) > 0L, all(analyses %in% names(outcome)),
  #           all(c("n", "Sigma", unname(outcome[analyses])) %in% names(results_df)))
  
  number <- function(x) as.numeric(if (is.factor(x)) as.character(x) else x)
  n_obs <- number(results_df$n)
  sigma <- number(results_df$Sigma)
  
  # The posted results_df stacks one block per settings row, in settings order.
  # Attach the design beta before sorting. An explicit degree column also lets
  # this function work on a results_df that has already been reordered.
  if ("heterogeneity_degree" %in% names(results_df)) {
    beta <- number(results_df$heterogeneity_degree)
  } else {
    stopifnot(length(setting_index) == nrow(results_df),
              all(setting_index %in% seq_len(nrow(settings))))
    design <- settings[sort(setting_index), , drop = FALSE]
    if (!isTRUE(all(n_obs == design$n & sigma == design$sigma))) {
      stop("The result rows do not match the settings blocks. Attach the original heterogeneity_degree to each replicate before sorting.")
    }
    beta <- number(design$heterogeneity_degree)
  }
  stopifnot(all(is.finite(n_obs)), all(is.finite(sigma)), all(is.finite(beta)))
  
  # Use the same decision under the null and under alternatives. For SB this
  # is the positive-lower-bound rule used in the supplied power calculation.
  decisions <- lapply(analyses, function(method) {
    x <- results_df[[outcome[[method]]]]
    if (method == "PM") {
      flags <- c("TRUE" = TRUE, "FALSE" = FALSE, "1" = TRUE, "0" = FALSE)
      reject <- unname(flags[toupper(trimws(as.character(x)))])
      if (any(!is.na(x) & is.na(reject))) stop("PM_test_reject must contain TRUE/FALSE or 1/0.")
    } else if (method == "SB") {
      reject <- number(x) > 0
    } else {
      p_value <- number(x)
      if (any(p_value < 0 | p_value > 1, na.rm = TRUE)) stop("Invalid p-value for ", method)
      reject <- p_value <= alpha
    }
    data.frame(n = n_obs, beta = beta, sigma = sigma, analysis = method, reject = reject)
  })
  
  rates <- bind_rows(decisions) %>%
    mutate(analysis = factor(analysis, levels = analyses)) %>%
    group_by(n, analysis, beta, sigma) %>%
    summarise(n_total = dplyr::n(), n_used = sum(!is.na(reject)),
              n_reject = sum(reject, na.rm = TRUE), .groups = "drop") %>%
    mutate(rejection_rate = if_else(n_used > 0, n_reject / n_used, NA_real_),
           panel = factor(if_else(beta == 0, "Type I error", "Power"),
                          levels = c("Type I error", "Power"))) %>%
    group_by(analysis, beta, sigma) %>%
    mutate(curve_id = cur_group_id()) %>%
    ungroup() %>%
    arrange(curve_id, n)
  
  if (any(rates$n_used < rates$n_total)) {
    warning("Some test decisions are missing. Rates use the available decisions. Check n_used and n_total in pm_summary.", call. = FALSE)
  }
  rates
}

pm_rejection_plots <- function(rates, alpha = 0.05, sigma_value = NULL) {
  sigmas <- sort(unique(rates$sigma))
  if (is.null(sigma_value) && length(sigmas) == 1L) sigma_value <- sigmas[[1L]]
  if (length(sigma_value) != 1L || !is.numeric(sigma_value) || !(sigma_value %in% sigmas)) {
    stop("Set pm_sigma before sourcing, or pass sigma_value to pm_rejection_plots(). ",
         "Available sigma values: ", paste(sigmas, collapse = ", "), call. = FALSE)
  }
  rates <- rates[rates$sigma == sigma_value, , drop = FALSE]
  
  method_lines <- c(SB = "solid", TwoStep = "dashed", PM = "dotdash", Dukes = "longdash")
  method_shapes <- c(SB = 16, TwoStep = 17, PM = 15, Dukes = 18)
  methods <- unique(as.character(rates$analysis))
  betas <- sort(unique(rates$beta))
  alternatives <- betas[betas != 0]
  alt_colours <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  stopifnot(length(alternatives) <= length(alt_colours))
  beta_colours <- c("0" = "#222222",
                    setNames(alt_colours[seq_along(alternatives)], as.character(alternatives)))
  coefficient_labels <- function(x) {
    formatC(as.numeric(as.character(x)), digits = 2, format = "f", drop0trailing = TRUE)
  }
  percent_labels <- function(x) {
    paste0(formatC(100 * x, digits = 1, format = "f", drop0trailing = TRUE), "%")
  }
  rates <- mutate(rates, beta_key = factor(beta, levels = betas))
  sample_sizes <- sort(unique(rates$n))
  
  draw <- function(data, title, combined = FALSE) {
    # Power spans [0, 1]. The level axis includes zero, the 5% target (at the
    # default alpha), and every observed level, even if there is substantial inflation.
    anchors <- distinct(data, panel) %>%
      mutate(n = min(sample_sizes),
             upper = if_else(panel == "Power", 1, min(1, 2 * alpha)))
    reference <- filter(anchors, panel == "Type I error") %>% mutate(target = alpha)
    #reference_note <- if (nrow(reference) > 0L) {
    #  paste0(" The dotted line marks the ", percent_labels(alpha), " reference level.")
    #} else ""
    #caption <- bquote(paste(sigma, " = ", .(sigma_value), ".", .(reference_note)))
    
    g <- ggplot(data, aes(x = n, y = rejection_rate)) +
      geom_blank(data = anchors, aes(x = n, y = upper), inherit.aes = FALSE) +
      geom_hline(data = reference, aes(yintercept = target), linetype = "dotted",
                 colour = "grey40", linewidth = 0.45, show.legend = FALSE) +
      geom_line(aes(group = curve_id, colour = beta_key, linetype = analysis),
                linewidth = 0.8) +
      geom_point(aes(colour = beta_key, shape = analysis), size = 2.4) +
      scale_colour_manual(name = expression(beta), values = beta_colours, labels = coefficient_labels) +
      scale_linetype_manual(name = "Test", values = method_lines, breaks = methods) +
      scale_shape_manual(name = "Test", values = method_shapes, breaks = methods) +
      scale_x_continuous(breaks = sample_sizes, labels = scales::label_comma(accuracy = 1),
                         expand = expansion(mult = 0.07)) +
      scale_y_continuous(limits = c(0, NA), labels = percent_labels,
                         breaks = function(limits) {
                           ticks <- pretty(limits, n = 4)
                           ticks[ticks >= 0 & ticks <= 1]
                         }, expand = expansion(mult = c(0.025, 0.045))) +
      guides(
        colour = guide_legend(order = 2, nrow = 1,
                              override.aes = list(shape = NA, linetype = "solid")),
        linetype = guide_legend(order = 1, nrow = 1,
                                override.aes = list(colour = "grey20")),
        shape = guide_legend(order = 1, nrow = 1)
      ) +
      labs(title = title, x = expression("Sample size (" * italic(n) * ")"),
           y = "Rejection rate") + #, caption = caption
      theme_classic(base_size = 11, base_family = "serif") +
      theme(
        panel.grid.major.x = element_line(colour = "grey89", linewidth = 0.3),
        panel.grid.major.y = element_line(colour = "grey89", linewidth = 0.3),
        panel.grid.minor = element_blank(), axis.text = element_text(colour = "grey15"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold", hjust = 0),
        #plot.caption = element_text(hjust = 0, size = 9),
        legend.position = "bottom", legend.box = "horizontal",
        legend.key.width = grid::unit(0.85, "cm"),
        panel.spacing.x = grid::unit(6, "mm")
      )
    # Faceting the full data builds one shared legend, including beta = 0.
    if (combined) g <- g + facet_wrap(vars(panel), nrow = 1, scales = "free_y")
    g
  }
  
  list(
    level = draw(filter(rates, beta == 0), "Type I error"),
    power = draw(filter(rates, beta != 0), "Power"),
    combined = draw(rates, "", combined = TRUE) #Testing for treatment effect heterogeneity
  )
}

# pm_alpha <- 0.05  # Match the level used to construct the SB interval and the PM decision.
# pm_summary <- pm_rejection_summary(results_df, settings, setting_index, analysis_type, pm_alpha)
# pm_plots <- pm_rejection_plots(pm_summary, pm_alpha, sigma_value = get0("pm_sigma", inherits = FALSE))
# print(pm_plots$combined)

# Use print(pm_plots$level) or print(pm_plots$power) to display a single panel.
# To count the curves drawn at the selected sigma, run
# count(distinct(pm_plots$combined$data, panel, curve_id), panel, name = "n_curves")
# Optional vector PDF export
# ggsave("pm_level_power.pdf", plot = pm_plots$combined, width = 190, height = 100,
#        units = "mm", device = grDevices::pdf, useDingbats = FALSE, bg = "white")
