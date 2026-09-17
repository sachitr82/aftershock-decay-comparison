#===============================================================================
# Temporal-form misspecification.
# NOTE: This script draws 1,000 joint posterior samples across all usable fitted
# models. On an Apple M4 Mac, the full run takes approximately 30 minutes.
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
# Load fit manifest and restrict to usable fits
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)

fit_manifest$fit_file <- file.path(fit_dir,
                                   paste0("truth_", fit_manifest$truth_form),
                                   fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)

misspec_manifest <- fit_manifest %>%
  filter(usable, truth_form != fitted_form)

stopifnot(nrow(misspec_manifest) == 481)

#-------------------------------------------------------------------------------
# Functional behaviour under misspecification
#-------------------------------------------------------------------------------

# Number of joint posterior draws
n_samp <- 1000
t_grid <- c(0, 10^seq(-3, log10(T_fit_end), length.out = 300))

functional_rows <- vector("list", nrow(misspec_manifest))

for (i in seq_len(nrow(misspec_manifest))) {
  
  # Loop through every misspecified fit
  row_i <- misspec_manifest[i, ]
  truth_form_i <- row_i$truth_form
  fitted_form_i <- row_i$fitted_form
  
  message(i, "/", nrow(misspec_manifest), " | truth ", truth_form_i,
          " | fit ", fitted_form_i, " | rep ", row_i$rep)
  
  # Load fitted model
  obj_i <- readRDS(row_i$fit_file)
  
  set.seed(3000000 + i)
  
  # Obtain pointwise posterior summary of g(t)
  functional_i <- ETAS.inlabru::posterior_temporal_summary(
    list(model.fit = obj_i$fit,
         link.functions = obj_i$link.functions, 
         form = fitted_form_i),
    t.eval = t_grid, n.samp = n_samp)$summary %>%
    filter(quantity == "g(t)")
  
  # True generating function
  true_g_i <- ETAS.inlabru::temporal_decay(
    dt = t_grid, theta = truths[[truth_form_i]], form = truth_form_i)
  
  # Truth and catalogue identifiers
  time_index_i <- match(functional_i$time, t_grid)
  
  functional_i$truth <- true_g_i[time_index_i]
  functional_i$truth_form <- truth_form_i
  functional_i$fitted_form <- fitted_form_i
  functional_i$rep <- row_i$rep
  
  functional_rows[[i]] <- functional_i %>% select(truth_form, fitted_form,
                                                  rep, time, truth, median)
}

functional_misspec_summary <- bind_rows(functional_rows)

#-------------------------------------------------------------------------------
# Total-variation distance under misspecification
#-------------------------------------------------------------------------------

trapz_num <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

# Window-normalise posterior-median and true curves, then compute TV distance
tv_misspec <- functional_misspec_summary %>%
  group_by(truth_form, fitted_form, rep) %>%
  arrange(time, .by_group = TRUE) %>%
  group_modify(~ {
    
    recovered_mass <- trapz_num(.x$time, .x$median)
    true_mass <- trapz_num(.x$time, .x$truth)
    
    h_recovered <- .x$median / recovered_mass
    h_true <- .x$truth / true_mass
    
    data.frame(tv_wrong = 0.5 * trapz_num(.x$time, abs(h_recovered - h_true)))
  }) %>%
  ungroup()

stopifnot(nrow(tv_misspec) == nrow(misspec_manifest),
          all(is.finite(tv_misspec$tv_wrong)),
          all(tv_misspec$tv_wrong >= 0),
          all(tv_misspec$tv_wrong <= 1))

# Pair with correct-specification recovery (within catalogue pairing)
correct_tv <- read.csv(
  file.path(recovery_table_dir, "functional_tv_by_catalogue.csv"),
  stringsAsFactors = FALSE) %>%
  rename(tv_correct = tv)

tv_paired <- tv_misspec %>%
  left_join(correct_tv, by = c("truth_form", "rep")) %>%
  mutate(delta_tv = tv_wrong - tv_correct)

stopifnot(nrow(tv_paired) == nrow(misspec_manifest),
          all(!is.na(tv_paired$tv_correct)),
          all(is.finite(tv_paired$delta_tv)))

#-------------------------------------------------------------------------------
# Functional misspecification summary
#-------------------------------------------------------------------------------

tv_misspec_summary <- tv_paired %>%
  group_by(truth_form, fitted_form) %>%
  summarise(
    n = n(),
    median_delta_tv = median(delta_tv),
    q10_delta_tv = quantile(delta_tv, 0.10),
    q90_delta_tv = quantile(delta_tv, 0.90),
    prop_wrong_better = mean(delta_tv < 0),
    .groups = "drop")

write.csv(tv_misspec_summary, file.path(misspecification_table_dir,
                                        "functional_tv_summary.csv"),
          row.names = FALSE)

#===============================================================================
# Posterior distribution of functional TV
#===============================================================================

#-------------------------------------------------------------------------------
# All usable fits
#-------------------------------------------------------------------------------

n_post_tv <- 1000

usable_manifest <- fit_manifest %>% filter(usable)

stopifnot(nrow(usable_manifest) == 300 + nrow(misspec_manifest))

#-------------------------------------------------------------------------------
# Trapezoidal integration weights
#-------------------------------------------------------------------------------

trapz_weights <- function(x) {
  
  dx <- diff(x)
  
  c(dx[1] / 2, (head(dx, -1) + tail(dx, -1)) / 2, tail(dx, 1) / 2)
}

trap_w <- trapz_weights(t_grid)

stopifnot(length(trap_w) == length(t_grid),
          abs(sum(trap_w) - (max(t_grid) - min(t_grid))) < 1e-8)

#-------------------------------------------------------------------------------
# True window-normalised temporal shapes
#-------------------------------------------------------------------------------

truth_h <- lapply(
  c("ou", "mse", "rate_state"),
  function(k) {
    
    g_true <- ETAS.inlabru::temporal_decay(
      dt = t_grid, theta = truths[[k]], form = k)
    
    mass_true <- sum(trap_w * g_true)
    g_true / mass_true
  }
)

names(truth_h) <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Posterior TV draws: normalise each joint draw and calculate TV from truth
#-------------------------------------------------------------------------------

posterior_tv_rows <- vector("list", nrow(usable_manifest))

for (i in seq_len(nrow(usable_manifest))) {
  
  row_i <- usable_manifest[i, ]
  
  truth_form_i <- row_i$truth_form
  fitted_form_i <- row_i$fitted_form
  
  message(i, "/", nrow(usable_manifest), " | posterior TV | truth ", 
          truth_form_i, " | fit ", fitted_form_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  # Reproducible joint posterior draws
  set.seed(4000000 + i)
  
  post_samp_i <- ETAS.inlabru::post_sampling(
    input.list = list(model.fit = obj_i$fit, link.functions = obj_i$link.functions,
      form = fitted_form_i), n.samp = n_post_tv, max.batch = n_post_tv)
  
  # Evaluate g(t) for each joint posterior draw
  g_draws_i <- vapply(
    seq_len(nrow(post_samp_i)),
    function(s) {
      
      theta_s <- as.list(post_samp_i[s, , drop = FALSE])
      
      ETAS.inlabru::temporal_decay(dt = t_grid, theta = theta_s,
                                    form = fitted_form_i)
    },
    numeric(length(t_grid)))
  
  # Normalise each posterior draw separately over the fitting window
  masses_i <- colSums(g_draws_i * trap_w)
  
  stopifnot(all(is.finite(masses_i)), all(masses_i > 0))
  
  h_draws_i <- sweep(g_draws_i, MARGIN = 2, STATS = masses_i, FUN = "/")
  
  # Known generating h(t)
  h_true_i <- truth_h[[truth_form_i]]
  
  # Calculate TV distance from the known generating shape for each draw
  
  abs_diff_i <- abs(sweep(h_draws_i, MARGIN = 1, STATS = h_true_i, FUN = "-"))
  
  tv_draws_i <- 0.5 * colSums(abs_diff_i * trap_w)
  
  stopifnot(length(tv_draws_i) == n_post_tv, all(is.finite(tv_draws_i)),
    all(tv_draws_i >= 0), all(tv_draws_i <= 1))
  
  posterior_tv_rows[[i]] <- data.frame(
    truth_form = truth_form_i,
    fitted_form = fitted_form_i,
    rep = row_i$rep,
    draw = seq_len(n_post_tv),
    tv = tv_draws_i)
  
  # Occasional memory cleanup
  rm(obj_i, post_samp_i, g_draws_i, h_draws_i, abs_diff_i)
  
  if (i %% 25 == 0) {
    gc()
  }
}

posterior_tv_draws <- bind_rows(posterior_tv_rows)

stopifnot(nrow(posterior_tv_draws) == nrow(usable_manifest) * n_post_tv)

# Summarise 1,000 TV values to one posterior interval per fit
posterior_tv_by_fit <- posterior_tv_draws %>%
  group_by(truth_form, fitted_form, rep) %>%
  summarise(
    posterior_q025_tv = quantile(tv, 0.025),
    posterior_median_tv = median(tv),
    posterior_q975_tv = quantile(tv, 0.975),
    .groups = "drop")

stopifnot(nrow(posterior_tv_by_fit) == nrow(usable_manifest))

#-------------------------------------------------------------------------------
# Posterior-TV figure
#-------------------------------------------------------------------------------

posterior_tv_plot_data <- posterior_tv_by_fit %>%
  group_by(truth_form, fitted_form) %>%
  arrange(posterior_median_tv, .by_group = TRUE) %>%
  mutate(catalogue_rank = row_number()) %>%
  ungroup() %>%
  mutate(truth_label = factor(truth_form, levels = c("ou", "mse", "rate_state"),
                              labels = c("OU truth", "MSE truth", "RS truth")),
    fitted_label = factor(fitted_form, levels = c("ou", "mse", "rate_state"),
                          labels = c("OU fit", "MSE fit", "RS fit")))

# Plot posterior interval
p_posterior_tv <- ggplot(posterior_tv_plot_data, aes(x = catalogue_rank,
                                                     y = posterior_median_tv)) +
  geom_linerange(aes(ymin = posterior_q025_tv, ymax = posterior_q975_tv),
                 alpha = 0.35, linewidth = 0.35) +
  geom_point(size = 0.7, alpha = 0.8) +
  facet_grid(truth_label ~ fitted_label) +
  labs(x = "Synthetic catalogue (ordered by posterior median TV)",
       y = "Posterior TV distance") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 20))

ggsave(file.path(misspecification_figure_dir, "posterior_functional_tv.pdf"),
       p_posterior_tv, width = 10, height = 8)

message("Saved misspecification tables to: ", misspecification_table_dir)
message("Saved misspecification figures to: ", misspecification_figure_dir)
message("Finished temporal-form misspecification analysis.")
