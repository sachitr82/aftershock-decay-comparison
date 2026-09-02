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

# Remove outputs from previous run
if (dir.exists(recovery_dir)) {
  unlink(recovery_dir, recursive = TRUE, force = TRUE)
}

dir.create(recovery_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------------------------------------
# Plotting labels
#-------------------------------------------------------------------------------

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

parameter_labels <- c(mu = "mu", K = "K", alpha = "alpha", c = "c", p = "p",
                      d = "d", rho = "rho", gamma = "gamma", B = "B", ta = "t[a]")

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"), 
                         stringsAsFactors = FALSE)

fit_manifest$fit_file <- file.path(main_fit_dir, 
                                   paste0("truth_", fit_manifest$truth_kernel),
                                   fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)
fit_manifest <- fit_manifest %>%
  mutate(usable = converged & !hit_max & !degenerate & !vb_aborted &
           !inla_failure & !nan_inf_logl)

stopifnot(all(fit_manifest$converged == fit_manifest$usable, na.rm = TRUE))

model_parameters <- list(
  ou = c("mu", "K", "alpha", "c", "p"),
  mse = c("mu", "K", "alpha", "d", "rho", "gamma"),
  rate_state = c("mu", "K", "alpha", "B", "ta"))

diagonal_manifest <- fit_manifest %>%
  filter(usable, truth_kernel == fitted_kernel)

stopifnot(nrow(diagonal_manifest) == 300)

#-------------------------------------------------------------------------------
# Extract posterior parameter distributions and summaries under correct specification
#-------------------------------------------------------------------------------

posterior_rows <- vector("list", nrow(diagonal_manifest))
marginal_rows <- vector("list", nrow(diagonal_manifest))

for (i in seq_len(nrow(diagonal_manifest))) {
  
  row_i <- diagonal_manifest[i, ]
  kernel_i <- row_i$fitted_kernel
  
  message(i, "/", nrow(diagonal_manifest), " | ", kernel_i,
          " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  post_out_i <- ETAS.inlabru::get_posterior_param(
    list(
      model.fit = obj_i$fit,
      link.functions = obj_i$link.functions,
      kernel = kernel_i))
  
  post_i <- post_out_i$post.summary %>%
    rename(parameter = param, q025 = q0.025, q975 = q0.975)
  
  post_i$truth <- vapply(
    post_i$parameter,
    function(p) as.numeric(obj_i$truth_parameters[[p]]),
    numeric(1))
  
  post_i$truth_kernel <- kernel_i
  post_i$rep <- row_i$rep
  
  posterior_rows[[i]] <- post_i %>%
    select(
      truth_kernel, rep, parameter, truth,
      mean, q025, median, q975)
    
    marg_i <- post_out_i$post.df
    marg_i$truth_kernel <- kernel_i
    marg_i$rep <- row_i$rep
    marg_i$truth <- vapply(marg_i$param, function(p)
      as.numeric(obj_i$truth_parameters[[p]]), numeric(1))
    
    marginal_rows[[i]] <- marg_i
}

posterior_fit_summary <- bind_rows(posterior_rows) %>%
  mutate(error = median - truth, squared_error = error^2, ci_width = q975 - q025,
         covered = q025 <= truth & q975 >= truth)
posterior_marginals <- bind_rows(marginal_rows)

print(posterior_fit_summary)

write.csv(posterior_fit_summary, file.path(recovery_dir, 
                                           "posterior_fit_summary_correct.csv"),
          row.names = FALSE)

saveRDS(posterior_marginals, file.path(recovery_dir,
                                       "posterior_marginals_correct_spec.rds"))

#-------------------------------------------------------------------------------
# Prior marginal densities
#-------------------------------------------------------------------------------

set.seed(123)
n_prior <- 50000

prior_rows <- lapply(names(model_parameters), function(k) {
  
  links_k <- make_links_P0(k)
  
  bind_rows(lapply(model_parameters[[k]], function(p) {
    
    link_name <- if (p == "c") "c_" else p
    prior_samp <- links_k[[link_name]](rnorm(n_prior))
    prior_density <- density(prior_samp, n = 512)
    
    data.frame(truth_kernel = k, param = p, x = prior_density$x,
               y = prior_density$y)
  }))
})

prior_marginals <- bind_rows(prior_rows)

#-------------------------------------------------------------------------------
# Posterior marginal recovery
#-------------------------------------------------------------------------------
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

truth_lines_marginal <- posterior_marginals %>%
  distinct(truth_kernel, param, truth)

posterior_ranges <- posterior_marginals %>%
  group_by(truth_kernel, param) %>%
  summarise(xmin = min(x), xmax = max(x), .groups = "drop")

prior_marginals_plot <- prior_marginals %>%
  left_join(posterior_ranges, by = c("truth_kernel", "param")) %>%
  filter(x >= xmin, x <= xmax)

p_marginals <- ggplot(posterior_marginals, aes(x = x, y = y, group = rep)) +
  geom_line(alpha = 0.08, linewidth = 0.3) +
  geom_vline(data = truth_lines_marginal, aes(xintercept = truth), linetype = "dashed") +
  facet_wrap(~ truth_kernel + param, scales = "free", labeller = labeller(
    truth_kernel = as_labeller(kernel_labels), 
    param = as_labeller(parameter_labels, label_parsed)))  +
  scale_x_continuous(breaks = facet_breaks, labels = facet_labels) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
  labs(x = "Parameter value", y = "Posterior density") +
theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 15),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    panel.grid.minor = element_blank())

p_marginals_with_priors <- ggplot(posterior_marginals, aes(x = x, y = y, group = rep)) +
  geom_line(alpha = 0.08, linewidth = 0.3) +
  geom_line(data = prior_marginals_plot, aes(x = x, y = y), inherit.aes = FALSE,
            linewidth = 0.9) +
  geom_vline(data = truth_lines_marginal, aes(xintercept = truth), linetype = "dashed") +
  facet_wrap(~ truth_kernel + param, scales = "free", labeller = labeller(
    truth_kernel = as_labeller(kernel_labels), 
    param = as_labeller(parameter_labels, label_parsed)))  +
  labs(x = "Parameter value", y = "Density") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

print(p_marginals)
print(p_marginals_with_priors)

ggsave(file.path(figure_dir, "parameter_posterior_overlays.pdf"), p_marginals,
       width = 10, height = 7)
ggsave(file.path(figure_dir, "parameter_posterior_overlays_with_priors.pdf"), 
       p_marginals_with_priors, width = 10, height = 7)

#-------------------------------------------------------------------------------
# Numerical parameter recovery summary
#-------------------------------------------------------------------------------

posterior_fit_summary_table <- posterior_fit_summary %>%
  mutate(
    parameter = ifelse(truth_kernel == "rate_state" & parameter == "B",
                       "1-B", parameter),
    truth = ifelse(parameter == "1-B", 1 - truth, truth),
    median_old = median,
    q025_old = q025,
    q975_old = q975,
    median = ifelse(parameter == "1-B", 1 - median_old, median_old),
    q025 = ifelse(parameter == "1-B", 1 - q975_old, q025_old),
    q975 = ifelse(parameter == "1-B", 1 - q025_old, q975_old)) %>%
  mutate(
    error = median - truth,
    squared_error = error^2,
    ci_width = q975 - q025,
    covered = q025 <= truth & q975 >= truth)

recovery_summary_table <- posterior_fit_summary_table %>%
  group_by(truth_kernel, parameter) %>%
  summarise(
    n = n(),
    truth = first(truth),
    mean_posterior_median = mean(median),
    bias = mean(error),
    rmse = sqrt(mean(squared_error)),
    relative_bias = bias / truth,
    relative_rmse = rmse / truth,
    coverage = mean(covered),
    median_ci_width = median(ci_width),
    relative_ci_width = median_ci_width/ truth,
    .groups = "drop")

print(recovery_summary_table)
write.csv(recovery_summary_table,
          file.path(recovery_dir, "correct_specification_recovery.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Functional recovery under correct specification
#-------------------------------------------------------------------------------

n_samp <- 1000
t_grid <- c(0, 10^seq(-3, log10(T_fit_end), length.out = 300))

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
    t.eval = t_grid, n.samp = n_samp)$summary %>%
    rename(q025 = q0.025, q975 = q0.975) %>%
    filter(quantity == "g(t)")
  
  true_g_i <- ETAS.inlabru::temporal_kernel(
    dt = t_grid, theta = truths[[kernel_i]], kernel = kernel_i)
  time_index_i <- match(functional_i$time, t_grid)
  
  functional_i$truth <- true_g_i[time_index_i]
  functional_i$truth_kernel <- kernel_i
  functional_i$rep <- row_i$rep
  
  functional_rows[[i]] <- functional_i %>% select(truth_kernel, rep, time, 
                                                  truth, q025, median, q975)
}

functional_fit_summary <- bind_rows(functional_rows)
functional_fit_summary <- functional_fit_summary %>%
  mutate(covered = q025 <= truth & q975 >= truth)

saveRDS(functional_fit_summary, file.path(recovery_dir, "functional_fit_summary.rds"))

#-------------------------------------------------------------------------------
# Plot recovery of temporal triggering functions
#-------------------------------------------------------------------------------

g_plot_data <- functional_fit_summary %>%
  filter(time > 0)

x_range_g <- range(g_plot_data$time)

p_functional_g <- g_plot_data%>%
  ggplot(aes(x = time, y = median, group = rep), linetype = "solid") +
  geom_line(alpha = 0.08, linewidth = 0.25) + 
  geom_line(data = g_plot_data %>%
      distinct(truth_kernel, time, truth),
    aes(x = time, y = truth), inherit.aes = FALSE, colour = "firebrick",
    linewidth = 1.1,linetype = "dashed") +
  facet_wrap(~truth_kernel,
             labeller = as_labeller(c(mse = "MSE", ou = "OU", rate_state = "RS"))) +
  scale_x_log10(labels = scales::label_log(), expand = expansion(mult = c(0, 0))) +
  scale_y_log10(breaks = 10^seq(-10, 0, by = 2), labels = scales::label_log()) +
  coord_cartesian(xlim = x_range_g, ylim = c(1e-10, NA)) +
  labs(x = "Time since parent event (days)", y = expression(g(t))) +
  theme_bw(base_size = 16) + theme(panel.spacing.x = unit(1.5, "lines"))

print(p_functional_g)

ggsave(file.path(figure_dir, "functional_recovery.pdf"),
       p_functional_g, width = 8, height = 4)

#-------------------------------------------------------------------------------
# Total-variation recovery of normalised temporal shape
#-------------------------------------------------------------------------------

trapz_num <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

tv_recovery <- functional_fit_summary %>%
  group_by(truth_kernel, rep) %>%
  arrange(time, .by_group = TRUE) %>%
  group_modify(~ {
    
    recovered_mass <- trapz_num(.x$time, .x$median)
    true_mass <- trapz_num(.x$time, .x$truth)
    
    h_recovered <- .x$median / recovered_mass
    h_true <- .x$truth / true_mass
    
    data.frame(tv = 0.5 * trapz_num(.x$time, abs(h_recovered - h_true)))
  }) %>%
  ungroup()

stopifnot(
  all(is.finite(tv_recovery$tv)),
  all(tv_recovery$tv >= 0),
  all(tv_recovery$tv <= 1))

tv_summary <- tv_recovery %>%
  group_by(truth_kernel) %>%
  summarise(n = n(),
    median_tv = median(tv),
    q10_tv = quantile(tv, 0.10),
    q90_tv = quantile(tv, 0.90),
    .groups = "drop")

print(tv_summary)

write.csv(tv_recovery, file.path(recovery_dir, "functional_tv_by_catalogue.csv"),
          row.names = FALSE)

write.csv(tv_summary, file.path(recovery_dir, "functional_tv_recovery.csv"),
          row.names = FALSE)
#-------------------------------------------------------------------------------
# Pointwise credible-interval coverage of g(t)
#-------------------------------------------------------------------------------

coverage_by_time <- functional_fit_summary %>%
  group_by(truth_kernel, time) %>%
  summarise(pointwise_coverage = mean(covered),.groups = "drop")

functional_coverage <- coverage_by_time %>%
  group_by(truth_kernel) %>%
  summarise(mean_coverage = mean(pointwise_coverage), .groups = "drop")

print(functional_coverage)
write.csv(functional_coverage, 
          file.path(recovery_dir, "functional_coverage.csv"),
          row.names = FALSE)