#===============================================================================
# Ridgecrest posterior temporal functions
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(ETAS.inlabru)
library(scales)
library(grid)
library(patchwork)

source(here("analyses", "ridgecrest", "00_design.R"))

candidate_forms <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Settings
#-------------------------------------------------------------------------------

n_samp <- 10000

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
  lapply(candidate_forms, function(k) {
    readRDS(file.path(baseline_fit_dir, paste0("fit_", k, ".rds")))
  }),
  candidate_forms)

#-------------------------------------------------------------------------------
# Joint posterior propagation through g(t) and h(t)
#-------------------------------------------------------------------------------

temporal_rows <- vector("list", length(candidate_forms))

for (i in seq_along(candidate_forms)) {
  
  form_i <- candidate_forms[i]
  obj_i <- fit_objects[[form_i]]
  
  message("Posterior temporal functions | ", form_i)
  
  set.seed(5000000 + i)
  
  # Joint posterior samples
  post_samp_i <- ETAS.inlabru::post_sampling(
    input.list = list(model.fit = obj_i$fit,
                      link.functions = obj_i$link.functions,
                      form = form_i), 
    n.samp = n_samp, max.batch = 1000)
  
  # Pointwise definition of posterior functions
  g_draws_i <- vapply(
    seq_len(nrow(post_samp_i)),
    function(s) {

      theta_s <- as.list(post_samp_i[s, , drop = FALSE])
      
      ETAS.inlabru::temporal_decay(dt = t_grid, theta = theta_s, form = form_i)
    },
    numeric(length(t_grid)))
  
  stopifnot(all(is.finite(g_draws_i)), all(g_draws_i >= 0))
  
  masses_i <- colSums(g_draws_i * trap_w)
  
  stopifnot(all(is.finite(masses_i)), all(masses_i > 0))
  
  # Normalised over fitting window
  h_draws_i <- sweep(g_draws_i, MARGIN = 2, STATS = masses_i, FUN = "/")
  
  # Pointwise posterior summaries
  g_summary_i <- data.frame(model = form_i,
                            quantity = "g(t)",
                            time = t_grid,
                            q025 = apply(g_draws_i, 1, quantile, probs = 0.025),
                            median = apply(g_draws_i, 1, median),
                            q975 = apply(g_draws_i, 1, quantile, probs = 0.975))
  
  h_summary_i <- data.frame(model = form_i,
                            quantity = "h(t)",
                            time = t_grid,
                            q025 = apply(h_draws_i, 1, quantile, probs = 0.025),
                            median = apply(h_draws_i, 1, median),
                            q975 = apply(h_draws_i, 1, quantile, probs = 0.975))
  
  temporal_rows[[i]] <- bind_rows(g_summary_i, h_summary_i)
  
  rm(post_samp_i, g_draws_i, h_draws_i)
  gc()
}

decay_levels <- c("Omori\u2013Utsu",
                  "Modified stretched exponential",
                  "Rate-state")

temporal_summary <- bind_rows(temporal_rows) %>%
  mutate(model = factor(model, levels = candidate_forms, labels = decay_levels),
         quantity = factor(quantity, levels = c("g(t)", "h(t)"),
                           labels = c("Unnormalised decay, g(t)",
                                      "Window-normalised decay, h(t)")))

#-------------------------------------------------------------------------------
# Save numerical summaries
#-------------------------------------------------------------------------------

write.csv(temporal_summary,
          file.path(temporal_function_table_dir,
                    "baseline_temporal_function_summary.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

# Match colours and line types used in the earlier form comparison
model_cols <- setNames(c("steelblue", "orange", "springgreen4"), decay_levels)
model_lty <- setNames(c("solid", "dashed", "dotted"), decay_levels)

# Tick labels as 10^x rather than 1e+x
log_lab <- scales::trans_format( "log10", scales::math_format(10^.x))

plot_data <- temporal_summary %>% filter(time > 0)

#-------------------------------------------------------------------------------
# Plot g(t) and h(t)
#-------------------------------------------------------------------------------

base_plot <- ggplot(plot_data,
                    aes(x = time, y = median, colour = model,
                        fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.22, colour = NA,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = 10^seq(-3, 3, by = 1), labels = log_lab,
                expand = expansion(mult = 0)) +
  scale_y_log10(labels = log_lab, expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(1e-3, T_fit_end), ylim = c(1e-7, NA)) +
  scale_colour_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  scale_linetype_manual(values = model_lty) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_line(linewidth = 0.15),
        legend.title = element_blank(),
        legend.key.width = unit(0.9, "lines"),
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.text = element_text(size = 11))

p_g <- base_plot %+% filter(plot_data, quantity == "Unnormalised decay, g(t)") +
  labs(x = "Time since parent event (days)", y = expression(g[k](t)))

p_h <- base_plot %+% filter(plot_data, quantity == "Window-normalised decay, h(t)") +
  labs(x = "Time since parent event (days)", y = expression(h[k](t)))

p_temporal <- (p_g | p_h) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 10, face = "plain"))

ggsave(file.path(temporal_function_figure_dir,
                 "ridgecrest-functional-posteriors.pdf"), 
       p_temporal, width = 7, height = 3.5)

message("Saved temporal-function summary to: ", temporal_function_table_dir)
message("Saved temporal-function figure to: ", temporal_function_figure_dir)
message("Finished Ridgecrest posterior temporal functions.")
