#===============================================================================
# Parameter and functional recovery under correct specification.
# NOTE: This script draws 1,000 joint posterior samples for each of 300 fitted
# catalogues. On an Apple M4 Mac, the full run takes approximately 30–60 minutes.
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(grid)
library(here)
library(ETAS.inlabru)

source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Plotting labels
#-------------------------------------------------------------------------------

form_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

parameter_labels <- c(mu = "mu", K = "K", alpha = "alpha", c = "c", p = "p",
                      d = "d", rho = "rho", gamma = "gamma", B = "B", ta = "t[a]")

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(fit_dir, "fit_manifest.csv"), 
                         stringsAsFactors = FALSE)

fit_manifest$fit_file <- file.path(fit_dir, 
                                   paste0("truth_", fit_manifest$truth_form),
                                   fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)

# Retain the 300 usable correctly specified fits
diagonal_manifest <- fit_manifest %>%
  filter(usable, truth_form == fitted_form)

stopifnot(nrow(diagonal_manifest) == 300)

#-------------------------------------------------------------------------------
# Extract posterior parameter distributions and summaries under correct specification
#-------------------------------------------------------------------------------

posterior_rows <- vector("list", nrow(diagonal_manifest))
marginal_rows <- vector("list", nrow(diagonal_manifest))

for (i in seq_len(nrow(diagonal_manifest))) {
  
  #-------------------------------------------------------------------------
  # Identify form and load saved object
  #-------------------------------------------------------------------------
  row_i <- diagonal_manifest[i, ]
  form_i <- row_i$fitted_form
  
  message(i, "/", nrow(diagonal_manifest), " | ", form_i,
          " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  #-------------------------------------------------------------------------
  # Transform from latent Gaussian to ETAS scale, summarise posterior
  #-------------------------------------------------------------------------
  post_out_i <- ETAS.inlabru::get_posterior_param(
    list(model.fit = obj_i$fit,
         link.functions = obj_i$link.functions,
         form = form_i))
  
  post_i <- post_out_i$post.summary %>%
    rename(parameter = param, q025 = q0.025, q975 = q0.975)
  
  #-------------------------------------------------------------------------
  # Attach known truth to each parameter
  #-------------------------------------------------------------------------
  post_i$truth <- vapply(post_i$parameter,
                         function(p) as.numeric(obj_i$truth_parameters[[p]]),
                         numeric(1))
  
  post_i$truth_form <- form_i
  post_i$rep <- row_i$rep
  
  posterior_rows[[i]] <- post_i %>%
    select(truth_form, rep, parameter, truth, q025, median, q975)
    
    marg_i <- post_out_i$post.df
    marg_i$truth_form <- form_i
    marg_i$rep <- row_i$rep
    marg_i$truth <- vapply(marg_i$param, function(p)
      as.numeric(obj_i$truth_parameters[[p]]), numeric(1))
    
    marginal_rows[[i]] <- marg_i
}

# Catalogue-level parameter recovery quantities
posterior_fit_summary <- bind_rows(posterior_rows) %>%
  mutate(error = median - truth, squared_error = error^2, ci_width = q975 - q025,
         covered = q025 <= truth & q975 >= truth)
posterior_marginals <- bind_rows(marginal_rows)

#-------------------------------------------------------------------------------
# Posterior marginal recovery
#-------------------------------------------------------------------------------

# Visualisation aids
facet_breaks <- function(x) {
  b <- pretty(x, n = 4)
  b <- b[b >= min(x) & b <= max(x)]
  
  if (length(b) >= 2) {
    c(first(b), last(b))
  } else {
    seq(min(x), max(x), length.out = 4)[2:3]
  }
}
facet_labels <- function(x) {
  if (all(x > 0.99 & x < 1.001)) {
    formatC(x, format = "f", digits = 5)
  } else if (diff(range(x)) < 0.1) {
    scales::label_number(accuracy = 0.001, trim = TRUE)(x)
  } else if (diff(range(x)) < 1) {
    scales::label_number(accuracy = 0.01, trim = TRUE)(x)
  } else {
    scales::label_number(accuracy = 1, trim = TRUE)(x)
  }
}

# Generating parameter value
truth_lines_marginal <- posterior_marginals %>%
  distinct(truth_form, param, truth)

# Overlay catalogue-level posterior marginals with the generating value
p_marginals <- ggplot(posterior_marginals, aes(x = x, y = y, group = rep)) +
  geom_line(alpha = 0.08, linewidth = 0.3) +
  geom_vline(data = truth_lines_marginal, aes(xintercept = truth), linetype = "dashed") +
  facet_wrap(~ truth_form + param, scales = "free", labeller = labeller(
    truth_form = as_labeller(form_labels), 
    param = as_labeller(parameter_labels, label_parsed)))  +
  scale_x_continuous(breaks = facet_breaks, labels = facet_labels) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
  labs(x = "Parameter value", y = "Posterior density") +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        panel.grid.minor = element_blank())

ggsave(file.path(recovery_figure_dir, "parameter_posterior_overlays.pdf"), p_marginals,
       width = 10, height = 7)

#-------------------------------------------------------------------------------
# Numerical parameter recovery summary across catalogues
#-------------------------------------------------------------------------------

posterior_fit_summary_table <- posterior_fit_summary %>%
  mutate(
    
    # Report RS recovery as 1-B because B lies close to one
    parameter = ifelse(truth_form == "rate_state" & parameter == "B",
                       "1-B", parameter),
    truth = ifelse(parameter == "1-B", 1 - truth, truth),
    median_old = median,
    q025_old = q025,
    q975_old = q975,
    median = ifelse(parameter == "1-B", 1 - median_old, median_old),
    q025 = ifelse(parameter == "1-B", 1 - q975_old, q025_old),
    q975 = ifelse(parameter == "1-B", 1 - q025_old, q975_old)) %>%
  mutate(error = median - truth,
         squared_error = error^2,
         ci_width = q975 - q025,
         covered = q025 <= truth & q975 >= truth)

recovery_summary_table <- posterior_fit_summary_table %>%
  group_by(truth_form, parameter) %>%
  summarise(n = n(),
            truth = first(truth),
            bias = mean(error),
            rmse = sqrt(mean(squared_error)),
            relative_bias = bias / truth,
            relative_rmse = rmse / truth,
            coverage = mean(covered),
            median_ci_width = median(ci_width),
            relative_ci_width = median_ci_width/ truth,
            .groups = "drop")

write.csv(recovery_summary_table,
          file.path(recovery_table_dir, "correct_specification_recovery.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Functional recovery under correct specification
#-------------------------------------------------------------------------------

# Draw 1,000 joint posterior samples and evaluate g(t) over the fitting window
n_samp <- 1000
t_grid <- c(0, 10^seq(-3, log10(T_fit_end), length.out = 300))

functional_rows <- vector("list", nrow(diagonal_manifest))

for (i in seq_len(nrow(diagonal_manifest))) {
  
  row_i <- diagonal_manifest[i, ]
  form_i <- row_i$fitted_form
  
  message(i, "/", nrow(diagonal_manifest), " | functional recovery | ",
          form_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  set.seed(2000000 + i)
  
  functional_i <- ETAS.inlabru::posterior_temporal_summary(
    list(model.fit = obj_i$fit,
         link.functions = obj_i$link.functions, 
         form = form_i),
    t.eval = t_grid, n.samp = n_samp)$summary %>%
    rename(q025 = q0.025, q975 = q0.975) %>%
    filter(quantity == "g(t)")
  
  true_g_i <- ETAS.inlabru::temporal_decay(dt = t_grid,
                                           theta = truths[[form_i]],
                                           form = form_i)
  time_index_i <- match(functional_i$time, t_grid)
  
  functional_i$truth <- true_g_i[time_index_i]
  functional_i$truth_form <- form_i
  functional_i$rep <- row_i$rep
  
  functional_rows[[i]] <- functional_i %>% select(truth_form, rep, time, 
                                                  truth, q025, median, q975)
}

functional_fit_summary <- bind_rows(functional_rows)
functional_fit_summary <- functional_fit_summary %>%
  mutate(covered = q025 <= truth & q975 >= truth)

#-------------------------------------------------------------------------------
# Plot recovery of temporal triggering functions
#-------------------------------------------------------------------------------

g_plot_data <- functional_fit_summary %>%
  filter(time > 0)

x_range_g <- range(g_plot_data$time)

p_functional_g <- g_plot_data%>%
  ggplot(aes(x = time, y = median, group = rep)) +
  geom_line(alpha = 0.08, linewidth = 0.25) + 
  geom_line(data = g_plot_data %>%
      distinct(truth_form, time, truth),
    aes(x = time, y = truth), inherit.aes = FALSE, colour = "firebrick",
    linewidth = 1.1,linetype = "dashed") +
  facet_wrap(~truth_form,
             labeller = as_labeller(c(mse = "MSE", ou = "OU", rate_state = "RS"))) +
  scale_x_log10(labels = scales::label_log(), expand = expansion(mult = c(0, 0))) +
  scale_y_log10(breaks = 10^seq(-10, 0, by = 2), labels = scales::label_log()) +
  coord_cartesian(xlim = x_range_g, ylim = c(1e-10, NA)) +
  labs(x = "Time since parent event (days)", y = expression(g(t))) +
  theme_bw(base_size = 16) + theme(panel.spacing.x = unit(1.5, "lines"))

ggsave(file.path(recovery_figure_dir, "functional_recovery.pdf"),
       p_functional_g, width = 8, height = 4)

#-------------------------------------------------------------------------------
# Total-variation recovery of normalised temporal shape
#-------------------------------------------------------------------------------

trapz_num <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

tv_recovery <- functional_fit_summary %>%
  group_by(truth_form, rep) %>%
  arrange(time, .by_group = TRUE) %>%
  group_modify(~ {
    
    # Normalise fitted and generating functions over the observation window
    # so TV compares temporal shape rather than total triggering mass.
    recovered_mass <- trapz_num(.x$time, .x$median)
    true_mass <- trapz_num(.x$time, .x$truth)
    
    h_recovered <- .x$median / recovered_mass
    h_true <- .x$truth / true_mass
    
    data.frame(tv = 0.5 * trapz_num(.x$time, abs(h_recovered - h_true)))
  }) %>%
  ungroup()

stopifnot(all(is.finite(tv_recovery$tv)),
          all(tv_recovery$tv >= 0),
          all(tv_recovery$tv <= 1))

tv_summary <- tv_recovery %>%
  group_by(truth_form) %>%
  summarise(n = n(),
            median_tv = median(tv),
            q10_tv = quantile(tv, 0.10),
            q90_tv = quantile(tv, 0.90),
            .groups = "drop")

write.csv(tv_recovery, file.path(recovery_table_dir, "functional_tv_by_catalogue.csv"),
          row.names = FALSE)

write.csv(tv_summary, file.path(recovery_table_dir, "functional_tv_recovery.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Pointwise credible-interval coverage of g(t)
#-------------------------------------------------------------------------------

# Calculate coverage across catalogues at each time, then average over the grid.
coverage_by_time <- functional_fit_summary %>%
  group_by(truth_form, time) %>%
  summarise(pointwise_coverage = mean(covered),.groups = "drop")

functional_coverage <- coverage_by_time %>%
  group_by(truth_form) %>%
  summarise(mean_coverage = mean(pointwise_coverage), .groups = "drop")

write.csv(functional_coverage, 
          file.path(recovery_table_dir, "functional_coverage.csv"),
          row.names = FALSE)

message("Saved parameter recovery tables to: ", recovery_table_dir)
message("Saved parameter recovery figures to: ", recovery_figure_dir)
message("Finished parameter and functional recovery analysis.")
