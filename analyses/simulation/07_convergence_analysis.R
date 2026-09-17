#===============================================================================
# Fit convergence diagnostics
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(scales)

source(here("analyses", "simulation", "00_design.R"))

candidate_forms <- c("ou", "mse", "rate_state")

form_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)

required_columns <- c("truth_form",
                      "fitted_form",
                      "rep",
                      "catalogue_seed",
                      "fit_seed",
                      "n_fit",
                      "runtime_minutes",
                      "fit_status",
                      "converged",
                      "hit_max",
                      "inla_failure",
                      "nan_inf_logl",
                      "vb_aborted",
                      "n_iter",
                      "min_sd",
                      "degenerate",
                      "usable",
                      "status",
                      "error",
                      "file")

missing_columns <- setdiff(required_columns, names(fit_manifest))
stopifnot(length(missing_columns) == 0)

#-------------------------------------------------------------------------------
# Validate production design
#-------------------------------------------------------------------------------

stopifnot(nrow(fit_manifest) == 3 * n_rep * 3,
          !anyDuplicated(fit_manifest[, c("truth_form", "rep", "fitted_form")]),
          setequal(unique(fit_manifest$truth_form), candidate_forms),
          setequal(unique(fit_manifest$fitted_form), candidate_forms))

cell_counts <- table(fit_manifest$truth_form, fit_manifest$fitted_form)
stopifnot(all(cell_counts[candidate_forms, candidate_forms] == n_rep))

fit_manifest$truth_form <- factor(fit_manifest$truth_form,
                                    levels = candidate_forms)
fit_manifest$fitted_form <- factor(fit_manifest$fitted_form,
                                     levels = candidate_forms)

#-------------------------------------------------------------------------------
# Overall convergence
#-------------------------------------------------------------------------------

overall_summary <- fit_manifest %>%
  summarise(n_total = n(),
            n_converged = sum(converged %in% TRUE),
            convergence_rate = mean(converged %in% TRUE),
            n_usable = sum(usable %in% TRUE),
            usable_rate = mean(usable %in% TRUE),
            n_hit_max = sum(hit_max %in% TRUE))

convergence_3x3 <- fit_manifest %>%
  group_by(truth_form, fitted_form) %>%
  summarise(n_total = n(),
            n_converged = sum(converged %in% TRUE),
            convergence_rate = mean(converged %in% TRUE),
            n_usable = sum(usable %in% TRUE),
            usable_rate = mean(usable %in% TRUE),
            n_hit_max = sum(hit_max %in% TRUE),
            .groups = "drop")

write.csv(overall_summary,
          file.path(convergence_table_dir, "convergence_overall.csv"),
          row.names = FALSE)

write.csv(convergence_3x3,
          file.path(convergence_table_dir, "convergence_3x3.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Convergence heatmap
#-------------------------------------------------------------------------------

p_convergence <- convergence_3x3 %>%
  ggplot(aes(x = fitted_form, y = truth_form, fill = convergence_rate)) +
  geom_tile(colour = "white", linewidth = 1) + 
  geom_text(aes(label = n_converged, colour = convergence_rate >= 0.75), 
            size = 8, fontface = "bold") +
  scale_x_discrete(labels = form_labels, expand = expansion(add = 0)) +
  scale_y_discrete(labels = form_labels, expand = expansion(add = 0)) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), 
                      labels = scales::percent) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), 
                      guide = "none") +
  labs(x = "Fitted temporal form", y = "Generating temporal form",
       fill = "Convergence rate") +
  coord_equal() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.title = element_text(size = 17), 
        axis.text = element_text(size = 15), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))

ggsave(file.path(convergence_figure_dir, "convergence_3x3.pdf"),
       p_convergence, width = 7, height = 5, dpi = 300)

#-------------------------------------------------------------------------------
# Iterations among converged fits
#-------------------------------------------------------------------------------

converged_fits <- fit_manifest %>%
  filter(converged, !is.na(n_iter))

iterations_3x3 <- converged_fits %>%
  group_by(truth_form, fitted_form) %>%
  summarise(n_converged = n(), median_iter = median(n_iter),
            q25_iter = quantile(n_iter, 0.25), q75_iter = quantile(n_iter, 0.75),
            p90_iter = quantile(n_iter, 0.90), .groups = "drop")

write.csv(iterations_3x3, file.path(convergence_table_dir, "iterations_3x3.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Iteration distributions
#-------------------------------------------------------------------------------

p_iterations <- converged_fits %>%
  ggplot(aes(x = fitted_form, y = n_iter)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ truth_form, nrow = 1, labeller = as_labeller(form_labels)) +
  scale_x_discrete(labels = form_labels) +
  labs(x = "Fitted temporal form", y = "Outer iterations to convergence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(convergence_figure_dir, "iterations_converged.png"),
       p_iterations, width = 10, height = 4, dpi = 300)

#-------------------------------------------------------------------------------
# Runtime
#-------------------------------------------------------------------------------

runtime_summary <- fit_manifest %>%
  group_by(truth_form, fitted_form) %>%
  summarise(n_fits = n(),
            median_runtime = median(runtime_minutes, na.rm = TRUE),
            q25_runtime = quantile(runtime_minutes, 0.25, na.rm = TRUE),
            q75_runtime = quantile(runtime_minutes, 0.75, na.rm = TRUE),
            p90_runtime = quantile(runtime_minutes, 0.90, na.rm = TRUE),
            .groups = "drop")

runtime_by_convergence <- fit_manifest %>%
  group_by(truth_form, fitted_form, converged) %>%
  summarise(n = n(),
            median_runtime = median(runtime_minutes, na.rm = TRUE),
            q25_runtime = quantile(runtime_minutes, 0.25, na.rm = TRUE),
            q75_runtime = quantile(runtime_minutes, 0.75, na.rm = TRUE),
            .groups = "drop")

write.csv(runtime_summary,
          file.path(convergence_table_dir, "runtime_3x3.csv"),
          row.names = FALSE)

write.csv(runtime_by_convergence,
          file.path(convergence_table_dir, "runtime_by_convergence.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Final summary
#-------------------------------------------------------------------------------

message("Converged: ", overall_summary$n_converged, "/",
        overall_summary$n_total,
        sprintf(" (%.1f%%)", 100 * overall_summary$convergence_rate))

message("Usable: ", overall_summary$n_usable, "/",
        overall_summary$n_total,
        sprintf(" (%.1f%%)", 100 * overall_summary$usable_rate))

message("Hit maximum iterations: ", overall_summary$n_hit_max)

message("Saved convergence tables to: ", convergence_table_dir)
message("Saved convergence figures to: ", convergence_figure_dir)
message("Finished fit convergence diagnostics.")