#===============================================================================
# Temporal-form misspecification
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
misspec_dir <- file.path(main_fit_dir, "misspecification")

dir.create(misspec_dir, recursive = TRUE, showWarnings = FALSE)

fit_manifest <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)
dir.create(misspec_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load fit manifest and restrict to usable fits (converged)
#-------------------------------------------------------------------------------

fit_manifest$fit_file <- file.path(main_fit_dir, 
                                   paste0("truth_", fit_manifest$truth_kernel),
                                   fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)

fit_manifest <- fit_manifest %>%
  mutate(usable = converged & !hit_max & !degenerate & !vb_aborted &
           !inla_failure & !nan_inf_logl)

misspec_manifest <- fit_manifest %>%
  filter(usable, truth_kernel != fitted_kernel)

stopifnot(nrow(misspec_manifest) == 473)

print(misspec_manifest %>% count(truth_kernel, fitted_kernel))

#-------------------------------------------------------------------------------
# Functional behaviour under misspecification
#-------------------------------------------------------------------------------

n_samp <- 1000
t_grid <- c(0, 10^seq(-3, log10(T_fit_end), length.out = 300))

functional_rows <- vector("list", nrow(misspec_manifest))

for (i in seq_len(nrow(misspec_manifest))) {
  
  row_i <- misspec_manifest[i, ]
  truth_kernel_i <- row_i$truth_kernel
  fitted_kernel_i <- row_i$fitted_kernel
  
  message(i, "/", nrow(misspec_manifest), " | truth ", truth_kernel_i,
          " | fit ", fitted_kernel_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  set.seed(3000000 + i)
  
  functional_i <- ETAS.inlabru::posterior_temporal_summary(
    list(model.fit = obj_i$fit, link.functions = obj_i$link.functions, 
         kernel = fitted_kernel_i),
    t.eval = t_grid, n.samp = n_samp)$summary %>%
    rename(q025 = q0.025, q975 = q0.975) %>%
    filter(quantity == "g(t)")
  
  true_g_i <- ETAS.inlabru::temporal_kernel(
    dt = t_grid, theta = truths[[truth_kernel_i]], kernel = truth_kernel_i)
  time_index_i <- match(functional_i$time, t_grid)
  
  functional_i$truth <- true_g_i[time_index_i]
  functional_i$truth_kernel <- truth_kernel_i
  functional_i$fitted_kernel <- fitted_kernel_i
  functional_i$rep <- row_i$rep
  
  functional_rows[[i]] <- functional_i %>% select(truth_kernel, fitted_kernel,
                                                  rep, time, truth, q025, median, q975)
}

functional_misspec_summary <- bind_rows(functional_rows)

saveRDS(functional_misspec_summary, file.path(misspec_dir, 
                                              "functional_misspec_summary.rds"))

#-------------------------------------------------------------------------------
# Total-variation distance under misspecification
#-------------------------------------------------------------------------------

trapz_num <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

tv_misspec <- functional_misspec_summary %>%
  group_by(truth_kernel, fitted_kernel, rep) %>%
  arrange(time, .by_group = TRUE) %>%
  group_modify(~ {
    
    recovered_mass <- trapz_num(.x$time, .x$median)
    true_mass <- trapz_num(.x$time, .x$truth)
    
    h_recovered <- .x$median / recovered_mass
    h_true <- .x$truth / true_mass
    
    data.frame(tv_wrong = 0.5 * trapz_num(.x$time, abs(h_recovered - h_true)))
  }) %>%
  ungroup()

stopifnot(nrow(tv_misspec) == 473, all(is.finite(tv_misspec$tv_wrong)),
          all(tv_misspec$tv_wrong >= 0), all(tv_misspec$tv_wrong <= 1))

# Pair with correct-specification recovery
correct_tv <- read.csv(
  file.path(recovery_dir, "functional_tv_by_catalogue.csv"),
  stringsAsFactors = FALSE) %>%
  rename(tv_correct = tv)

tv_paired <- tv_misspec %>%
  left_join(correct_tv, by = c("truth_kernel", "rep")) %>%
  mutate(delta_tv = tv_wrong - tv_correct)

stopifnot(nrow(tv_paired) == 473, all(!is.na(tv_paired$tv_correct)),
          all(is.finite(tv_paired$delta_tv)))

#-------------------------------------------------------------------------------
# Functional misspecification summary
#-------------------------------------------------------------------------------

tv_misspec_summary <- tv_paired %>%
  group_by(truth_kernel, fitted_kernel) %>%
  summarise(
    n = n(),
    
    median_tv_wrong = median(tv_wrong),
    q10_tv_wrong = quantile(tv_wrong, 0.10),
    q90_tv_wrong = quantile(tv_wrong, 0.90),
    
    median_tv_correct = median(tv_correct),
    
    median_delta_tv = median(delta_tv),
    q10_delta_tv = quantile(delta_tv, 0.10),
    q90_delta_tv = quantile(delta_tv, 0.90),
    
    prop_wrong_better = mean(delta_tv < 0),
    
    .groups = "drop")

print(tv_misspec_summary)

write.csv(tv_paired, file.path(misspec_dir, "functional_tv_by_catalogue.csv"),
          row.names = FALSE)

write.csv(tv_misspec_summary, file.path(misspec_dir, "functional_tv_summary.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Example where misspecified MSE recovers OU truth better than correct OU fit
#-------------------------------------------------------------------------------

example_case <- tv_paired %>%
  filter(truth_kernel == "ou", fitted_kernel == "mse", delta_tv < 0) %>%
  arrange(delta_tv) %>%
  slice(1)

print(example_case)

functional_correct <- readRDS(file.path(recovery_dir, "functional_fit_summary.rds"))

rep_example <- example_case$rep

correct_example <- functional_correct %>%
  filter(truth_kernel == "ou",rep == rep_example) %>%
  arrange(time)

misspec_example <- functional_misspec_summary %>%
  filter(truth_kernel == "ou", fitted_kernel == "mse", rep == rep_example) %>%
  arrange(time)

true_mass <- trapz_num(correct_example$time, correct_example$truth)

correct_mass <- trapz_num(correct_example$time, correct_example$median)

mse_mass <- trapz_num(misspec_example$time, misspec_example$median)

comparison_data <- bind_rows(
  data.frame(time = correct_example$time, h = correct_example$truth / true_mass,
             curve = "Generating OU"),
  data.frame(time = correct_example$time, h = correct_example$median / correct_mass,
             curve = "Correct OU fit"),
  data.frame(time = misspec_example$time, h = misspec_example$median / mse_mass,
             curve = "Misspecified MSE fit")
)

p_example <- comparison_data %>%
  filter(time > 0) %>%
  ggplot(aes(x = time, y = h, colour = curve, linetype = curve)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Generating OU" = "firebrick",
                                 "Correct OU fit" = "black",
                                 "Misspecified MSE fit" = "steelblue")) +
  scale_linetype_manual(values = c("Generating OU" = "dashed",
                                   "Correct OU fit" = "solid", 
                                   "Misspecified MSE fit" = "solid")) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10(labels = scales::label_log()) +
  labs(x = "Time since parent event (days)", y = expression(h(t)),
       colour = NULL, linetype = NULL) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

print(p_example)

#-------------------------------------------------------------------------------
# Posterior distribution of functional TV
#-------------------------------------------------------------------------------

n_post_tv <- 1000

usable_manifest <- fit_manifest %>% filter(usable)

stopifnot(nrow(usable_manifest) == 773)

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
    
    g_true <- ETAS.inlabru::temporal_kernel(
      dt = t_grid, theta = truths[[k]], kernel = k)
    
    mass_true <- sum(trap_w * g_true)
    g_true / mass_true
  }
)

names(truth_h) <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Pre-flight check on one fit
#-------------------------------------------------------------------------------

test_row <- usable_manifest[1, ]

test_obj <- readRDS(test_row$fit_file)

set.seed(4000001)

test_samp <- ETAS.inlabru::post_sampling(
  input.list = list(model.fit = test_obj$fit, 
                    link.functions = test_obj$link.functions,
                    kernel = test_row$fitted_kernel),
  n.samp = 5, max.batch = 5)

print(test_samp)
print(names(test_samp))

#-------------------------------------------------------------------------------
# Posterior TV draws
#-------------------------------------------------------------------------------

posterior_tv_rows <- vector("list", nrow(usable_manifest))

for (i in seq_len(nrow(usable_manifest))) {
  
  row_i <- usable_manifest[i, ]
  
  truth_kernel_i <- row_i$truth_kernel
  fitted_kernel_i <- row_i$fitted_kernel
  
  message(i, "/", nrow(usable_manifest), " | posterior TV | truth ", 
          truth_kernel_i, " | fit ", fitted_kernel_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  # Reproducible joint posterior draws
  set.seed(4000000 + i)
  
  post_samp_i <- ETAS.inlabru::post_sampling(
    input.list = list(model.fit = obj_i$fit, link.functions = obj_i$link.functions,
      kernel = fitted_kernel_i), n.samp = n_post_tv, max.batch = n_post_tv)
  
  # ----------------------------------------------------------
  # Propagate every joint posterior draw through g(t)
  #
  # Matrix dimensions:
  # rows    = times
  # columns = posterior draws
  # ----------------------------------------------------------
  
  g_draws_i <- vapply(
    seq_len(nrow(post_samp_i)),
    function(s) {
      
      theta_s <- as.list(post_samp_i[s, , drop = FALSE])
      
      ETAS.inlabru::temporal_kernel(dt = t_grid, theta = theta_s,
                                    kernel = fitted_kernel_i)
    },
    numeric(length(t_grid)))
  
  # ----------------------------------------------------------
  # Normalise each posterior draw separately over [0, T]
  # ----------------------------------------------------------
  
  masses_i <- colSums(g_draws_i * trap_w)
  
  stopifnot(all(is.finite(masses_i)), all(masses_i > 0))
  
  h_draws_i <- sweep(g_draws_i, MARGIN = 2, STATS = masses_i, FUN = "/")
  
  # Known generating h(t)
  h_true_i <- truth_h[[truth_kernel_i]]
  
  # ----------------------------------------------------------
  # TV(draw, truth)
  # ----------------------------------------------------------
  
  abs_diff_i <- abs(
    sweep(h_draws_i, MARGIN = 1, STATS = h_true_i, FUN = "-"))
  
  tv_draws_i <- 0.5 * colSums(abs_diff_i * trap_w)
  
  stopifnot(length(tv_draws_i) == n_post_tv, all(is.finite(tv_draws_i)),
    all(tv_draws_i >= 0), all(tv_draws_i <= 1))
  
  posterior_tv_rows[[i]] <- data.frame(
    truth_kernel = truth_kernel_i,
    fitted_kernel = fitted_kernel_i,
    rep = row_i$rep,
    draw = seq_len(n_post_tv),
    tv = tv_draws_i)
  
  # Occasional memory cleanup
  rm(obj_i, post_samp_i, g_draws_i, h_draws_i,  abs_diff_i)
  
  if (i %% 25 == 0) {
    gc()
  }
}

posterior_tv_draws <- bind_rows(posterior_tv_rows)

stopifnot(nrow(posterior_tv_draws) == 773 * n_post_tv)
saveRDS(posterior_tv_draws, file.path(misspec_dir,
                                      "posterior_functional_tv_draws.rds"))

posterior_tv_by_fit <- posterior_tv_draws %>%
  group_by(truth_kernel, fitted_kernel, rep) %>%
  summarise(
    posterior_mean_tv = mean(tv),
    
    posterior_q025_tv = quantile(tv, 0.025),
    posterior_q10_tv = quantile(tv, 0.10),
    
    posterior_median_tv = median(tv),
    
    posterior_q90_tv = quantile(tv, 0.90),
    posterior_q975_tv = quantile(tv, 0.975),
    
    posterior_iqr_tv = IQR(tv),
    
    .groups = "drop")

stopifnot(nrow(posterior_tv_by_fit) == 773)

write.csv(posterior_tv_by_fit, 
          file.path(misspec_dir, "posterior_functional_tv_by_fit.csv"),
  row.names = FALSE)

posterior_tv_summary <- posterior_tv_by_fit %>%
  group_by(
    truth_kernel,
    fitted_kernel
  ) %>%
  summarise(
    n = n(),
    
    median_posterior_tv =
      median(posterior_median_tv),
    
    q10_posterior_median_tv =
      quantile(posterior_median_tv, 0.10),
    
    q90_posterior_median_tv =
      quantile(posterior_median_tv, 0.90),
    
    median_80_width =
      median(posterior_q90_tv - posterior_q10_tv),
    
    median_width =
      median(posterior_q975_tv - posterior_q025_tv),
    
    .groups = "drop"
  )

print(posterior_tv_summary)

write.csv(posterior_tv_summary, 
          file.path(misspec_dir, "posterior_functional_tv_summary.csv"), 
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Posterior-TV figure
#-------------------------------------------------------------------------------

posterior_tv_plot_data <- posterior_tv_by_fit %>%
  group_by(truth_kernel, fitted_kernel) %>%
  arrange(posterior_median_tv, .by_group = TRUE) %>%
  mutate(catalogue_rank = row_number()) %>%
  ungroup() %>%
  mutate(truth_label = factor(truth_kernel, levels = c("ou", "mse", "rate_state"),
                              labels = c("OU truth", "MSE truth", "RS truth")),
    fitted_label = factor(fitted_kernel, levels = c("ou", "mse", "rate_state"),
                          labels = c("OU fit", "MSE fit", "RS fit")))

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

print(p_posterior_tv)

ggsave(file.path(misspec_dir, "posterior_functional_tv.pdf"),
       p_posterior_tv, width = 10, height = 8)
