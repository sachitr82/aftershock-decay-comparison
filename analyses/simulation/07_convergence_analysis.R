#===============================================================================
# Fit convergence diagnostics
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(tidyr)

source(here("analyses", "simulation", "00_design.R"))

candidate_kernels <- c("ou", "mse", "rate_state")

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

#-------------------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------------------

main_fit_dir <- file.path(fit_dir, "main_simulation")
diagnostic_dir <- file.path(main_fit_dir, "diagnostics")
diagnostic_table_dir <- file.path(diagnostic_dir, "tables")
diagnostic_figure_dir <- file.path(diagnostic_dir, "figures")

dir.create(diagnostic_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_figure_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)

required_columns <- c("truth_kernel", "fitted_kernel", "rep", "catalogue_seed", 
                      "fit_seed", "n_fit", "runtime_minutes", "converged", 
                      "hit_max", "n_iter", "status", "error", "file")

missing_columns <- setdiff(required_columns, names(fit_manifest))
stopifnot(length(missing_columns) == 0)

#-------------------------------------------------------------------------------
# Validate production design
#-------------------------------------------------------------------------------

stopifnot(nrow(fit_manifest) == 3 * n_rep * 3)
stopifnot(!anyDuplicated(fit_manifest[, c("truth_kernel", "rep", "fitted_kernel")]))
stopifnot(setequal(unique(fit_manifest$truth_kernel), candidate_kernels))
stopifnot(setequal(unique(fit_manifest$fitted_kernel), candidate_kernels))

cell_counts <- table(fit_manifest$truth_kernel, fit_manifest$fitted_kernel)
stopifnot(all(cell_counts[candidate_kernels, candidate_kernels] == n_rep))

fit_manifest$truth_kernel <- factor(fit_manifest$truth_kernel,
                                    levels = candidate_kernels)
fit_manifest$fitted_kernel <- factor(fit_manifest$fitted_kernel,
                                     levels = candidate_kernels)

fit_manifest <- fit_manifest %>%
  mutate(specification = if_else(truth_kernel == fitted_kernel,
                                 "Correctly specified", "Misspecified"))

#-------------------------------------------------------------------------------
# Overall convergence
#-------------------------------------------------------------------------------

overall_summary <- fit_manifest %>%
  summarise(n_total = n(), n_converged = sum(converged),
            convergence_rate = mean(converged), n_hit_max = sum(hit_max),
            median_iter = median(n_iter), median_runtime = median(runtime_minutes))

print(overall_summary)

write.csv(overall_summary, file.path(diagnostic_table_dir, "convergence_overall.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Convergence across generating--fitted combinations
#-------------------------------------------------------------------------------

convergence_3x3 <- fit_manifest %>%
  group_by(truth_kernel, fitted_kernel) %>%
  summarise(n_total = n(), n_converged = sum(converged),
            convergence_rate = mean(converged), n_hit_max = sum(hit_max),
            .groups = "drop")

print(convergence_3x3)

write.csv(convergence_3x3, file.path(diagnostic_table_dir, "convergence_3x3.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Compact convergence matrix
#-------------------------------------------------------------------------------

convergence_matrix <- convergence_3x3 %>%
  mutate(convergence = sprintf("%d/%d (%.1f%%)", n_converged, n_total,
                               100 * convergence_rate)) %>%
  select(truth_kernel, fitted_kernel, convergence) %>%
  pivot_wider(names_from = fitted_kernel, values_from = convergence)

print(convergence_matrix)

write.csv(convergence_matrix,
          file.path(diagnostic_table_dir, "convergence_matrix.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Convergence heatmap
#-------------------------------------------------------------------------------

p_convergence <- convergence_3x3 %>%
  ggplot(aes(x = fitted_kernel, y = truth_kernel, fill = convergence_rate)) +
  geom_tile(colour = "white", linewidth = 1) + 
  geom_text(aes(label = n_converged, colour = convergence_rate >= 0.75), 
            size = 8, fontface = "bold") +
  scale_x_discrete(labels = kernel_labels, expand = expansion(add = 0)) +
  scale_y_discrete(labels = kernel_labels, expand = expansion(add = 0)) +
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

print(p_convergence)

ggsave(file.path(diagnostic_figure_dir, "convergence_3x3.pdf"),
       p_convergence, width = 7, height = 5, dpi = 300)

#-------------------------------------------------------------------------------
# Marginal convergence summaries
#-------------------------------------------------------------------------------

convergence_by_fit <- fit_manifest %>%
  group_by(fitted_kernel) %>%
  summarise(n_total = n(), n_converged = sum(converged),
            convergence_rate = mean(converged), n_hit_max = sum(hit_max),
            .groups = "drop")

convergence_by_truth <- fit_manifest %>%
  group_by(truth_kernel) %>%
  summarise(n_total = n(), n_converged = sum(converged),
            convergence_rate = mean(converged), n_hit_max = sum(hit_max),
            .groups = "drop")

convergence_by_specification <- fit_manifest %>%
  group_by(specification) %>%
  summarise(n_total = n(), n_converged = sum(converged),
            convergence_rate = mean(converged), n_hit_max = sum(hit_max),
            .groups = "drop")

print(convergence_by_fit)
print(convergence_by_truth)
print(convergence_by_specification)

write.csv(convergence_by_fit,
          file.path(diagnostic_table_dir, "convergence_by_fitted_model.csv"),
          row.names = FALSE)
write.csv(convergence_by_truth,
          file.path(diagnostic_table_dir, "convergence_by_generating_model.csv"),
          row.names = FALSE)
write.csv(convergence_by_specification,
          file.path(diagnostic_table_dir, "convergence_by_specification.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Iterations among converged fits
#-------------------------------------------------------------------------------

converged_fits <- fit_manifest %>%
  filter(converged, !is.na(n_iter))

iterations_3x3 <- converged_fits %>%
  group_by(truth_kernel, fitted_kernel) %>%
  summarise(n_converged = n(), median_iter = median(n_iter),
            q25_iter = quantile(n_iter, 0.25), q75_iter = quantile(n_iter, 0.75),
            p90_iter = quantile(n_iter, 0.90), .groups = "drop")

print(iterations_3x3)

write.csv(iterations_3x3, file.path(diagnostic_table_dir, "iterations_3x3.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Iteration distributions
#-------------------------------------------------------------------------------

p_iterations <- converged_fits %>%
  ggplot(aes(x = fitted_kernel, y = n_iter)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ truth_kernel, nrow = 1, labeller = as_labeller(kernel_labels)) +
  scale_x_discrete(labels = kernel_labels) +
  labs(x = "Fitted temporal form", y = "Outer iterations to convergence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_iterations)

ggsave(file.path(diagnostic_figure_dir, "iterations_converged.png"),
       p_iterations, width = 10, height = 4, dpi = 300)

#-------------------------------------------------------------------------------
# Catalogue size and convergence
#-------------------------------------------------------------------------------

catalogue_size_summary <- fit_manifest %>%
  mutate(convergence_class = if_else(converged, "Converged", "Did not converge")) %>%
  group_by(fitted_kernel, convergence_class) %>%
  summarise(n_fits = n(), median_n_fit = median(n_fit),
            q25_n_fit = quantile(n_fit, 0.25), q75_n_fit = quantile(n_fit, 0.75),
            .groups = "drop")

print(catalogue_size_summary)

write.csv(catalogue_size_summary,
          file.path(diagnostic_table_dir, "catalogue_size_by_convergence.csv"),
          row.names = FALSE)

p_size_iterations <- converged_fits %>%
  ggplot(aes(x = n_fit, y = n_iter)) +
  geom_point(alpha = 0.35) +
  facet_grid(truth_kernel ~ fitted_kernel,
             labeller = labeller(truth_kernel = as_labeller(kernel_labels),
                                 fitted_kernel = as_labeller(kernel_labels))) +
  labs(x = "Number of catalogue events", y = "Outer iterations to convergence") +
  theme_bw()

print(p_size_iterations)

ggsave(file.path(diagnostic_figure_dir, "catalogue_size_vs_iterations.png"),
       p_size_iterations, width = 9, height = 8, dpi = 300)

#-------------------------------------------------------------------------------
# Runtime
#-------------------------------------------------------------------------------

runtime_summary <- fit_manifest %>%
  group_by(truth_kernel, fitted_kernel) %>%
  summarise(n_fits = n(), median_runtime = median(runtime_minutes, na.rm = TRUE),
            q25_runtime = quantile(runtime_minutes, 0.25, na.rm = TRUE),
            q75_runtime = quantile(runtime_minutes, 0.75, na.rm = TRUE),
            p90_runtime = quantile(runtime_minutes, 0.90, na.rm = TRUE),
            .groups = "drop")

print(runtime_summary)

runtime_by_convergence <- fit_manifest %>%
  group_by(truth_kernel, fitted_kernel, converged) %>%
  summarise(n = n(), median_runtime = median(runtime_minutes), 
            q25_runtime = quantile(runtime_minutes, 0.25), 
            q75_runtime = quantile(runtime_minutes, 0.75), 
            .groups = "drop")

print(runtime_by_convergence)

write.csv(runtime_summary, file.path(diagnostic_table_dir, "runtime_3x3.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Final summary
#-------------------------------------------------------------------------------

cat("Converged:", overall_summary$n_converged, "/", overall_summary$n_total,
    sprintf("(%.1f%%)\n", 100 * overall_summary$convergence_rate))
cat("Hit maximum iterations:", overall_summary$n_hit_max, "\n")