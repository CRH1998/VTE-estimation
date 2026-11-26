#-------------------------------------------------------------------------------
# Process simulation results
#-------------------------------------------------------------------------------

# Load simulations
load("C:/Users/brf337/Desktop/VTE/VTE-estimation/MainFiles/resl.rda")


# One list element per setting, with replicates in columns
results_df <- bind_rows(resl)


source("heterogeneity_plots.R")


# Summarize simulations into dataframe and plot
pm_alpha <- 0.05  # Match the level used to construct the SB interval and the PM decision.
pm_summary <- pm_rejection_summary(results_df = results_df, 
                                   settings = settings, 
                                   analysis_type = c("SB","TwoStep", "PM"), 
                                   alpha = pm_alpha)
pm_plots <- pm_rejection_plots(pm_summary, pm_alpha)
print(pm_plots$combined)


# Save plot
# ggsave(
#   "smooth_sim_results21112026.pdf",
#   plot = pm_plots$combined,
#   width = 200, height = 100, units = "mm",
#   device = grDevices::pdf,
#   useDingbats = FALSE, bg = "white"
# )


# write.table(results_df, file = "rough_sim_results21112026.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
# ggsave(
#   "rough_sim_results21112026.pdf",
#   plot = pm_plots$combined,
#   width = 200, height = 100, units = "mm",
#   device = grDevices::pdf,
#   useDingbats = FALSE, bg = "white"
# )










