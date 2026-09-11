#===============================================================================
# Ridgecrest posterior temporal functions
#===============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(ETAS.inlabru)

#-------------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

ridgecrest_fit_dir <- here("outputs", "ridgecrest", "baseline")

ridgecrest_summary_dir <- here(ridgecrest_fit_dir, "summaries")

ridgecrest_figure_dir <- here(ridgecrest_fit_dir, "figures")

dir.create(ridgecrest_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ridgecrest_figure_dir, recursive = TRUE, showWarnings = FALSE)

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

#-------------------------------------------------------------------------------
# Settings
#-------------------------------------------------------------------------------

n_samp <- 10000
T_fit_end <- 2371

t_grid <- c(0, 10^seq(-3, log10(T_fit_end), length.out = 300))

#-------------------------------------------------------------------------------
# Trapezoidal integration weights
#-------------------------------------------------------------------------------

trapz_weights <- function(x) {
  dx <- diff(x)
  c(dx[1] / 2, (head(dx, -1) + tail(dx, -1)) / 2, tail(dx, 1) / 2)
}

trap_w <- trapz_weights(t_grid)

stopifnot(length(trap_w) == length(t_grid))

#-------------------------------------------------------------------------------
# Load baseline fits
#-------------------------------------------------------------------------------

fit_objects <- setNames(
  lapply(candidate_kernels, function(k) {
    readRDS(file.path(ridgecrest_fit_dir, paste0("fit_", k, ".rds")))
  }),
  candidate_kernels)

#-------------------------------------------------------------------------------
# Joint posterior propagation through g(t) and h(t)
#-------------------------------------------------------------------------------

temporal_rows <- vector("list", length(candidate_kernels))

for (i in seq_along(candidate_kernels)) {
  
  kernel_i <- candidate_kernels[i]
  obj_i <- fit_objects[[kernel_i]]
  
  message("Posterior temporal functions | ", kernel_i)
  
  set.seed(5000000 + i)
  
  post_samp_i <- ETAS.inlabru::post_sampling(
    input.list = list(model.fit = obj_i$fit,
                      link.functions = obj_i$link.functions,
                      kernel = kernel_i), 
    n.samp = n_samp, max.batch = 1000)
  
  g_draws_i <- vapply(
    seq_len(nrow(post_samp_i)),
    function(s) {

      theta_s <- as.list(post_samp_i[s, , drop = FALSE])
      
      ETAS.inlabru::temporal_kernel(dt = t_grid, theta = theta_s, kernel = kernel_i)
    },
    numeric(length(t_grid)))
  
  stopifnot(all(is.finite(g_draws_i)), all(g_draws_i >= 0))
  
  masses_i <- colSums(g_draws_i * trap_w)
  
  stopifnot(all(is.finite(masses_i)), all(masses_i > 0))
  
  h_draws_i <- sweep(g_draws_i, MARGIN = 2, STATS = masses_i, FUN = "/")
  
  # Pointwise posterior summaries
  g_summary_i <- data.frame(model = kernel_i,
                            quantity = "g(t)",
                            time = t_grid,
                            q025 = apply(g_draws_i, 1, quantile, probs = 0.025),
                            median = apply(g_draws_i, 1, median),
                            q975 = apply(g_draws_i, 1, quantile, probs = 0.975))
  
  h_summary_i <- data.frame(model = kernel_i,
                            quantity = "h(t)",
                            time = t_grid,
                            q025 = apply(h_draws_i, 1, quantile, probs = 0.025),
                            median = apply(h_draws_i, 1, median),
                            q975 = apply(h_draws_i, 1, quantile, probs = 0.975))
  
  temporal_rows[[i]] <- bind_rows(g_summary_i, h_summary_i)
  
  rm(post_samp_i, g_draws_i, h_draws_i)
  gc()
}

temporal_summary <- bind_rows(temporal_rows) %>%
  mutate(model = factor(model, levels = candidate_kernels,
                        labels = c("OU", "MSE", "RS")),
         quantity = factor(quantity, levels = c("g(t)", "h(t)"),
                      labels = c("Unnormalised decay, g(t)",
                                 "Window-normalised decay, h(t)")))

#-------------------------------------------------------------------------------
# Save numerical summaries
#-------------------------------------------------------------------------------

saveRDS(temporal_summary,
        file.path(ridgecrest_summary_dir, "baseline_temporal_function_summary.rds"))

write.csv(temporal_summary,
          file.path(ridgecrest_summary_dir, "baseline_temporal_function_summary.csv"), 
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

model_levels <- c("OU", "MSE", "RS")

# Match colours and line types used in the earlier kernel comparison
model_cols <- c("OU"  = "steelblue", "MSE" = "orange", "RS"  = "springgreen4")

model_lty <- c("OU"  = "solid", "MSE" = "dashed", "RS"  = "dotted")

# Tick labels as 10^x rather than 1e+x
log_lab <- scales::trans_format( "log10", scales::math_format(10^.x))

plot_data <- temporal_summary %>% filter(time > 0)

#-------------------------------------------------------------------------------
# Plot g(t) and h(t)
#-------------------------------------------------------------------------------

p_temporal <- plot_data %>%
  ggplot(aes(x = time, y = median, colour = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.22, colour = NA,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  facet_wrap(~ quantity, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = 10^seq(-3, 3, by = 1), labels = log_lab,
                expand = expansion(mult = 0)) +
  scale_y_log10(labels = log_lab,  expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(1e-3, T_fit_end), ylim = c(1e-7, NA)) +
  scale_colour_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  scale_linetype_manual(values = model_lty) +
  labs(x = "Time since parent event (days)", y = NULL, colour = "Temporal form",
       linetype = "Temporal form") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_line(linewidth = 0.15),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(1.1, "lines"))

print(p_temporal)
ggsave(here(ridgecrest_figure_dir, "ridgecrest-functional-posteriors.pdf"), 
       p_temporal, width = 7, height = 3.5)
