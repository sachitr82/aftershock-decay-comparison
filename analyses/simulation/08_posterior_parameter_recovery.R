#===============================================================================
# Posterior parameter recovery
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(ETAS.inlabru)
source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------------------

main_fit_dir <- file.path(fit_dir, "main_simulation")
recovery_dir <- file.path(main_fit_dir, "parameter_recovery")
figure_dir <- file.path(recovery_dir, "figures")

dir.create(recovery_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"), 
                         stringsAsFactors = FALSE)

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")
shared_parameters <- c("mu", "K", "alpha")
model_parameters <- list(ou = c("mu", "K", "alpha", "c", "p"), 
                         mse = c("mu", "K", "alpha", "d", "rho", "gamma"), 
                         rate_state = c("mu", "K", "alpha", "B", "ta"))

fit_manifest$fit_file <- file.path(main_fit_dir, 
                                   paste0("truth_", fit_manifest$truth_kernel),
                                   fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)
stopifnot(all(c("truth_kernel", "fitted_kernel", "rep", "converged", "file") %in% names(fit_manifest)))
stopifnot(all(file.exists(fit_manifest$fit_file[fit_manifest$converged])))

#-------------------------------------------------------------------------------
# Extract posterior parameter summaries from converged fits
#-------------------------------------------------------------------------------

manifest_converged <- fit_manifest %>% filter(converged)

posterior_rows <- vector("list", nrow(manifest_converged))

for (i in seq_len(nrow(manifest_converged))) {
  
  row_i <- manifest_converged[i, ]
  truth_i <- row_i$truth_kernel
  fitted_i <- row_i$fitted_kernel
  params_i <- if (truth_i == fitted_i) model_parameters[[fitted_i]] else shared_parameters
  
  message(i, "/", nrow(manifest_converged), " | ", truth_i, " -> ", fitted_i,
          " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  post_i <- ETAS.inlabru::get_posterior_param(
    list(
      model.fit = obj_i$fit, link.functions = obj_i$link.functions, 
      kernel = fitted_i)
    )$post.summary
  
  post_i <- post_i %>% filter(param %in% params_i) %>%
    rename(parameter = param, q025 = q0.025, q975 = q0.975)
  
  stopifnot(nrow(post_i) == length(params_i))
  stopifnot(setequal(post_i$parameter, params_i))
  
  post_i$truth <- vapply(post_i$parameter, function(p) 
    as.numeric(obj_i$truth_parameters[[p]]), numeric(1))
  post_i$truth_kernel <- truth_i
  post_i$fitted_kernel <- fitted_i
  post_i$rep <- row_i$rep
  
  posterior_rows[[i]] <- post_i %>% 
    select(truth_kernel, fitted_kernel, rep, parameter, truth, mean, q025, median, q975)
  
  rm(obj_i, post_i)
  gc(verbose = FALSE)
}

posterior_fit_summary <- bind_rows(posterior_rows)
posterior_fit_summary <- posterior_fit_summary %>%
  mutate(error = median - truth, squared_error = error^2, ci_width = q975 - q025,
         covered = q025 <= truth & q975 >= truth)

write.csv(posterior_fit_summary, file.path(recovery_dir, "posterior_fit_summary.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Parameter recovery under correct specification
#-------------------------------------------------------------------------------

recovery_summary <- posterior_fit_summary %>%
  filter(truth_kernel == fitted_kernel) %>%
  group_by(truth_kernel, parameter) %>%
  summarise(n = n(), truth = first(truth), mean_posterior_median = mean(median),
            bias = mean(error), rmse = sqrt(mean(squared_error)),
            coverage = mean(covered),
            coverage_mcse = sqrt(coverage * (1 - coverage) / n),
            median_ci_width = median(ci_width), .groups = "drop")

print(recovery_summary)
write.csv(recovery_summary,
          file.path(recovery_dir, "correct_specification_recovery.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Check number of correctly specified converged fits
#-------------------------------------------------------------------------------

recovery_counts <- posterior_fit_summary %>%
  filter(truth_kernel == fitted_kernel) %>%
  distinct(truth_kernel, rep) %>%
  count(truth_kernel, name = "n_converged")

print(recovery_counts)

#-------------------------------------------------------------------------------
# Diagnostic plot of posterior medians under correct specification
#-------------------------------------------------------------------------------

recovery_plot_data <- posterior_fit_summary %>%
  filter(truth_kernel == fitted_kernel)

p_recovery <- ggplot(recovery_plot_data, aes(x = factor(rep), y = median)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(aes(yintercept = truth), linetype = "dashed") +
  facet_wrap(~ truth_kernel + parameter, scales = "free_y") +
  labs(x = "Synthetic catalogue", y = "Posterior median") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())

print(p_recovery)
ggsave(file.path(figure_dir, "parameter_recovery_diagnostic.pdf"),
       p_recovery, width = 10, height = 7)

#-------------------------------------------------------------------------------
# Shared parameters under temporal-form misspecification
#-------------------------------------------------------------------------------

shared_parameter_summary <- posterior_fit_summary %>%
  filter(parameter %in% shared_parameters) %>%
  group_by(truth_kernel, fitted_kernel, parameter) %>%
  summarise(n = n(), truth = first(truth), mean_posterior_median = mean(median), 
            bias = mean(error), rmse = sqrt(mean(squared_error)), 
            truth_in_95ci = mean(covered), median_ci_width = median(ci_width), 
            .groups = "drop")

print(shared_parameter_summary)
write.csv(shared_parameter_summary, 
          file.path(recovery_dir, "shared_parameter_misspecification.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Functional recovery under correct specification
#-------------------------------------------------------------------------------

n_samp <- 1000
t_grid <- 10^seq(-3, log10(T_fit_end), length.out = 300)

diagonal_manifest <- fit_manifest %>% filter(converged, truth_kernel == fitted_kernel)

functional_rows <- vector("list", nrow(diagonal_manifest))

for (i in seq_len(nrow(diagonal_manifest))) {
  
  row_i <- diagonal_manifest[i, ]
  kernel_i <- row_i$fitted_kernel
  
  message(i, "/", nrow(diagonal_manifest), " | functional recovery | ",
          kernel_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  set.seed(2000000 + i)
  
  functional_i <- ETAS.inlabru::posterior_temporal_summary(
    list(model.fit = obj_i$fit, link.functions = obj_i$link.functions, 
         kernel = kernel_i),
    t.eval = t_grid, n.samp = n_samp)$summary
  
  functional_i <- functional_i %>% rename(q025 = q0.025, q975 = q0.975)
  
  true_g_i <- ETAS.inlabru::temporal_kernel(
    dt = t_grid, theta = truths[[kernel_i]], kernel = kernel_i)
  time_index_i <- match(functional_i$time, t_grid)
  
  functional_i$truth <- ifelse(functional_i$quantity == "g(t)",
                               true_g_i[time_index_i], 
                               truths[[kernel_i]]$K * true_g_i[time_index_i])
  
  functional_i$truth_kernel <- kernel_i
  functional_i$rep <- row_i$rep
  
  functional_rows[[i]] <- functional_i %>% select(truth_kernel, rep, quantity,
                                                  time, truth, q025, median, q975)
  
  rm(obj_i, functional_i)
  gc(verbose = FALSE)
}

functional_fit_summary <- bind_rows(functional_rows)
functional_fit_summary <- functional_fit_summary %>%
  mutate(covered = q025 <= truth & q975 >= truth)

saveRDS(functional_fit_summary, file.path(recovery_dir, "functional_fit_summary.rds"))

#-------------------------------------------------------------------------------
# Summarise functional recovery across synthetic catalogues
#-------------------------------------------------------------------------------

functional_summary <- functional_fit_summary %>%
  group_by(truth_kernel, quantity, time) %>%
  summarise(n = n(), truth = first(truth),
            median_recovered = median(median), 
            q10_recovered = quantile(median, 0.10),
            q90_recovered = quantile(median, 0.90),
            coverage = mean(covered), .groups = "drop")

saveRDS(functional_summary, file.path(recovery_dir, "functional_recovery_summary.rds"))

#-------------------------------------------------------------------------------
# Plot recovery of temporal triggering functions
#-------------------------------------------------------------------------------

p_functional <- functional_summary %>%
  filter(quantity == "g(t)") %>%
  ggplot(aes(x = time, y = median_recovered)) +
  geom_ribbon(aes(ymin = q10_recovered, ymax = q90_recovered), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_line(aes(y = truth), linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ truth_kernel, labeller = as_labeller(kernel_labels)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Time since parent event (days)", y = "Temporal triggering function") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

print(p_functional)
ggsave(file.path(figure_dir, "functional_recovery.pdf"),
       p_functional, width = 8, height = 4)

#-------------------------------------------------------------------------------
# Pointwise credible-interval coverage of temporal functions
#-------------------------------------------------------------------------------

functional_coverage <- functional_summary %>%
  group_by(truth_kernel, quantity) %>%
  summarise(mean_coverage = mean(coverage), min_coverage = min(coverage), 
            max_coverage = max(coverage), .groups = "drop")

print(functional_coverage)
write.csv(functional_coverage, 
          file.path(recovery_dir, "functional_coverage.csv"),
          row.names = FALSE)