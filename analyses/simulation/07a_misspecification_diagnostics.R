#===============================================================================
# Misspecification diagnostics (why no RS fits to MSE truth converged)
# NOTE: This script performs 20 targeted RS refits and writes approximately
# 200 MB of fitted model objects to disk. Takes 2-4 hours to run
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(here)
library(ggplot2)
library(dplyr)
library(patchwork)
library(inlabru)

source(here("analyses", "simulation", "00_design.R"))
source(here("src", "fit_helpers", "fit_diagnostics.R"))
source(here("src", "fit_helpers", "fit_temporal_etas.R"))

#-------------------------------------------------------------------------------
# Load production fits and catalogues
#-------------------------------------------------------------------------------

fm <- read.csv(file.path(fit_dir, "fit_manifest.csv"),
               stringsAsFactors = FALSE)

simulation_manifest <- read.csv(file.path(catalogue_dir, "simulation_manifest.csv"),
                                stringsAsFactors = FALSE)

fm$fit_file <- file.path(fit_dir, paste0("truth_", fm$truth_form), fm$file)

rs_link <- make_links_P0("rate_state")

#-------------------------------------------------------------------------------
# Diagnostic 1: terminal rate-state parameter geometry
# RS fits to OU-, MSE-, and RS-generated catalogues. 
#-------------------------------------------------------------------------------

# Retain all saved RS fits, including non-usable fits, to examine where
# optimisation terminates in RS parameter space.
rs_rows <- fm %>%
  filter(fitted_form == "rate_state",
         status != "error",
         !is.na(file))

# Terminal values from non-converged fits are used only as numerical diagnostics,
# not as valid posterior parameter estimates.
rs_diag <- do.call(rbind, lapply(seq_len(nrow(rs_rows)), function(i) {
  r <- rs_rows[i, ]
  o <- readRDS(r$fit_file)
  
  # Transform terminal latent-scale posterior means to the physical RS scale.
  B <- rs_link$B(o$fit$summary.fixed["th.B", "mean"])
  ta <- rs_link$ta(o$fit$summary.fixed["th.ta", "mean"])
  
  data.frame(truth_form = r$truth_form,
             rep = r$rep,
             n_fit = r$n_fit,
             fit_status = r$fit_status,
             converged = r$converged,
             hit_max = r$hit_max,
             inla_failure = r$inla_failure,
             nan_inf_logl = r$nan_inf_logl,
             vb_aborted = r$vb_aborted,
             degenerate = r$degenerate,
             usable = r$usable,
             min_sd = r$min_sd,
             B = B,
             one_minus_B = 1 - B,
             ta = ta)
}))


rs_diag$truth_form <- factor(rs_diag$truth_form, 
                               levels = c("ou", "mse", "rate_state"), 
                               labels = c("OU truth", "MSE truth", "RS truth"))
rs_diag$fit_result <- case_when(rs_diag$degenerate %in% TRUE ~ "Degenerate",
                                rs_diag$usable %in% TRUE ~ "Usable",
                                TRUE ~ "Non-usable")

rs_diag$fit_result <- factor(rs_diag$fit_result,
                             levels = c("Usable", "Non-usable", "Degenerate"))


rs_parameter_summary <- rs_diag %>%
  group_by(truth_form, fit_result) %>%
  summarise(n = n(),
            median_B = median(B, na.rm = TRUE),
            median_1mB = median(one_minus_B, na.rm = TRUE),
            median_ta = median(ta, na.rm = TRUE),
            q25_ta = quantile(ta, 0.25, na.rm = TRUE),
            q75_ta = quantile(ta, 0.75, na.rm = TRUE),
            .groups = "drop")

# Save per-fit terminal RS parameter values and fit diagnostics
write.csv(rs_diag,
          file.path(misspec_diagnostic_table_dir, "rs_parameter_geometry.csv"),
          row.names = FALSE)

# Save grouped summaries of terminal RS parameter values
write.csv(rs_parameter_summary,
          file.path(misspec_diagnostic_table_dir, "rs_parameter_summary.csv"),
          row.names = FALSE)

p_ta <- ggplot(rs_diag %>% filter(is.finite(ta), ta > 0), 
               aes(truth_form, ta, shape = fit_result)) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 2) +
  scale_y_log10() +
  labs(x = NULL, y = expression(t[a]~"(days)"), title = expression("Terminal "*t[a])) +
  theme_bw() +
  theme(legend.position = "bottom")

p_B <- ggplot(rs_diag %>% filter(is.finite(one_minus_B), one_minus_B > 0), 
              aes(truth_form, one_minus_B, shape = fit_result)) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 2) +
  scale_y_log10() +
  labs(x = NULL, y = expression(1-B), title = expression("Terminal "*(1-B))) +
  theme_bw() +
  theme(legend.position = "bottom")

# Plot per-fit terminal t_a and 1-B values by truth and fit outcome
p_rs_geometry <- (p_ta | p_B ) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(file.path(misspec_diagnostic_figure_dir, "rs_parameter_geometry.pdf"), 
       p_rs_geometry, width = 8, height = 4)

#-------------------------------------------------------------------------------
# Diagnostic 2: Bidirectional RS--MSE functional approximation
# Can RS approximate the MSE truth at all?
#-------------------------------------------------------------------------------

# Log-spaced time grid
T <- T_fit_end
tg <- c(0, 10^seq(-5, log10(T), length.out = 2000))

# Total variation function
trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
tv <- function(h1, h2) 0.5 * trapz(tg, abs(h1 - h2))

# Window-normalise each temporal decay function for shape-only comparison
norm_shape <- function(theta, form) {
  g <- temporal_decay(dt = tg, theta = theta, form = form)
  G <- temporal_decay_integral(a = 0, b = T, theta = theta, form = form)
  g / G
}

# True generating functions
h_mse_true <- norm_shape(list(d = d_true, rho = rho_true, gamma = gamma_true), "mse")
h_rs_true <- norm_shape(list(B = B_true, ta = ta_true), "rate_state")


# Minimise finite-window TV between MSE truth and an RS approximation
obj_rs <- function(z) {
  B <- plogis(z[1]); ta <- exp(z[2])
  h <- norm_shape(list(B = B, ta = ta), "rate_state")
  tv(h_mse_true, h)
}

# Use multiple starts to reduce sensitivity to local optima
rs_starts <- rbind(c(qlogis(B_true), log(ta_true)), c(qlogis(0.999), log(50)), 
                   c(qlogis(0.9999), log(100)), c(qlogis(0.99999),log(500)),
                   c(qlogis(0.999999), log(2000)))

# Optimise and find optimising values
rs_opts <- lapply(seq_len(nrow(rs_starts)), function(i) 
  optim(rs_starts[i, ], obj_rs, method = "Nelder-Mead", control = list(maxit = 5000)))

rs_best <- rs_opts[[which.min(vapply(rs_opts, function(x) x$value, numeric(1)))]]
stopifnot(rs_best$convergence == 0)

# RS parameter values which minimised discrepancy with true normalised MSE
B_best <- plogis(rs_best$par[1])
ta_best <- exp(rs_best$par[2])
tv_rs_to_mse <- rs_best$value

# Minimise finite-window TV between RS truth and an MSE approximation
obj_mse <- function(z) {
  d <- exp(z[1]); rho <- exp(z[2]); gamma <- plogis(z[3])
  h <- norm_shape(list(d = d, rho = rho, gamma = gamma), "mse")
  tv(h_rs_true, h)
}

# Same as before
mse_starts <- rbind(c(log(d_true), log(rho_true), qlogis(gamma_true)), 
                    c(log(0.03), log(1), qlogis(0.2)), 
                    c(log(0.1), log(2), qlogis(0.4)), 
                    c(log(0.01), log(0.5), qlogis(0.1)), 
                    c(log(0.3), log(3), qlogis(0.7)))
mse_opts <- lapply(seq_len(nrow(mse_starts)), function(i) 
  optim(mse_starts[i, ], obj_mse, method = "Nelder-Mead", control = list(maxit = 5000)))
mse_best <- mse_opts[[which.min(vapply(mse_opts, function(x) x$value, numeric(1)))]]
stopifnot(mse_best$convergence == 0)

d_best <- exp(mse_best$par[1])
rho_best <- exp(mse_best$par[2])
gamma_best <- plogis(mse_best$par[3])
tv_mse_to_rs <- mse_best$value

# Final data frame containing parameter values that provided the closest match
rs_mse_approximation <- data.frame(
  direction = c("MSE truth -> RS", "RS truth -> MSE"),
  tv = c(tv_rs_to_mse, tv_mse_to_rs),
  B = c(B_best, NA),
  one_minus_B = c(1 - B_best, NA),
  ta = c(ta_best, NA),
  d = c(NA, d_best),
  rho = c(NA, rho_best),
  gamma = c(NA, gamma_best)
)

write.csv(rs_mse_approximation,
          file.path(misspec_diagnostic_table_dir,
                    "rs_mse_functional_approximation.csv"),
          row.names = FALSE)

h_rs_best <- norm_shape(list(B = B_best, ta = ta_best), "rate_state")
h_mse_best <- norm_shape(list(d = d_best, rho = rho_best, gamma = gamma_best), "mse")

plot_dat <- rbind(data.frame(t = tg[-1], h = h_mse_true[-1], 
                             comparison = "MSE truth -> RS fit", form = "Truth"),
                  data.frame(t = tg[-1], h = h_rs_best[-1], 
                             comparison = "MSE truth -> RS fit", 
                             form = "Best approximation"), 
                  data.frame(t = tg[-1], h = h_rs_true[-1],
                             comparison = "RS truth -> MSE fit", form = "Truth"),
                  data.frame(t = tg[-1], h = h_mse_best[-1], 
                             comparison = "RS truth -> MSE fit", form = "Best approximation"))

# Plot true shapes against their best alternative-form approximations
p_rs_mse <- ggplot(plot_dat, aes(t, h, linetype = form)) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~comparison) +
  labs(x = "Time since triggering event (days)", 
       y = "Normalised temporal triggering", linetype = NULL) +
  theme_bw()

ggsave(file.path(misspec_diagnostic_figure_dir,
                 "rs_mse_functional_approximation.pdf"),
       p_rs_mse, width = 9, height = 4)

#-------------------------------------------------------------------------------
# Diagnostic 2b: OU-to-RS control (repeat above but with OU truth and RS alternatives)
#-------------------------------------------------------------------------------

h_ou_true <- norm_shape(list(c = c_true, p = p_true), "ou")

obj_rs_ou <- function(z) {
  B <- plogis(z[1]); ta <- exp(z[2])
  h <- norm_shape(list(B = B, ta = ta), "rate_state")
  tv(h_ou_true, h)
}

rs_ou_opts <- lapply(seq_len(nrow(rs_starts)), function(i)
  optim(rs_starts[i, ], obj_rs_ou, method = "Nelder-Mead",
        control = list(maxit = 5000)))

rs_ou_best <- rs_ou_opts[[which.min(vapply(rs_ou_opts, function(x) x$value, numeric(1)))]]
stopifnot(rs_ou_best$convergence == 0)

B_ou_best <- plogis(rs_ou_best$par[1])
ta_ou_best <- exp(rs_ou_best$par[2])

tv_rs_to_ou <- rs_ou_best$value

rs_control_summary <- data.frame(
  truth_form = c("ou", "mse"),
  tv = c(tv_rs_to_ou, tv_rs_to_mse),
  B = c(B_ou_best, B_best),
  one_minus_B = c(1 - B_ou_best, 1 - B_best),
  ta = c(ta_ou_best, ta_best),
  B_prior_z = c(qlogis(B_ou_best) - 7.5,
                qlogis(B_best) - 7.5),
  ta_prior_z = c((log(ta_ou_best) - log(200)) / 0.5,
                 (log(ta_best) - log(200)) / 0.5)
)

write.csv(rs_control_summary,
          file.path(misspec_diagnostic_table_dir,
                    "rs_functional_control.csv"),
          row.names = FALSE)
h_rs_ou_best <- norm_shape(list(B = B_ou_best, ta = ta_ou_best), "rate_state")

plot_ou_rs <- rbind(data.frame(t = tg[-1], h = h_ou_true[-1], form = "OU truth"),
                    data.frame(t = tg[-1], h = h_rs_ou_best[-1], form = "Best RS approximation"))

p_ou_rs <- ggplot(plot_ou_rs, aes(t, h, linetype = form)) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "OU truth -> RS approximation",
       x = "Time since triggering event (days)",
       y = "Normalised temporal triggering", linetype = NULL) +
  theme_bw()

ggsave(file.path(misspec_diagnostic_figure_dir,
                 "ou_rs_functional_approximation.pdf"),
       p_ou_rs, width = 6, height = 4)

#-------------------------------------------------------------------------------
# Diagnostic 3: Initialisation sensitivity (does this impact convergence?)
#-------------------------------------------------------------------------------

# Extract all MSE truth RS fits
rs_mse_rows <- fm %>%
  filter(truth_form == "mse",
         fitted_form == "rate_state")

# Check none converged
stopifnot(nrow(rs_mse_rows) == n_rep,
          !any(rs_mse_rows$converged %in% TRUE))

# Extract 5
test_rows <- rs_mse_rows %>%
  filter(fit_status == "hit_max", !degenerate) %>%
  arrange(rep) %>%
  slice_head(n = 5)

stopifnot(nrow(test_rows) == 5)

# Reruns fits under same conditions except with initial value set to the
# values obtained from Diagnostic 2
warm_results <- lapply(seq_len(nrow(test_rows)), function(i) {
  r <- test_rows[i, ]
  
  sim_i <- simulation_manifest %>%
    filter(form == "mse", rep == r$rep)
  
  catalogue_i <- readRDS(file.path(catalogue_dir, sim_i$file))$catalogue
  
  catalogue_i <- prepare_temporal_catalogue(catalogue_i,
                                            T1 = T_fit_start,
                                            T2 = T_fit_end)
  
  set.seed(r$fit_seed)
  
  link_i <- make_links_P0("rate_state")
  inv_i <- make_inverse_links_P0("rate_state")
  bru_i <- make_bru_options_P0("rate_state", rel_tol = fit_control$rel_tol, 
                               max_iter = fit_control$max_iter)
  
  bru_i$bru_initial$th.B <- inv_i$B(B_best)
  bru_i$bru_initial$th.ta <- inv_i$ta(ta_best)
  bru_i$control.compute <- list(config = TRUE)
  
  start_time <- Sys.time()
  
  fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue_i,
                                       M0 = M0, 
                                       T1 = T_fit_start,
                                       T2 = T_fit_end, 
                                       link.functions = link_i, 
                                       coef.t. = temporal_binning$coef.t, 
                                       delta.t. = temporal_binning$delta.t, 
                                       N.max. = temporal_binning$N.max, 
                                       bru.opt = bru_i,
                                       form = "rate_state")
  
  runtime_i <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  diag_i <- fit_diagnostics(fit_i)
  usable_i <- fit_is_usable(diag_i)
  
  saveRDS(fit_i,
          file.path(misspec_diagnostic_fit_dir,
                    sprintf("mse_rep_%04d_rs_warm_start.rds", r$rep)))
  
  data.frame(rep = r$rep,
             n_fit = nrow(catalogue_i),
             baseline_status = r$fit_status,
             baseline_usable = r$usable,
             baseline_iter = r$n_iter,
             baseline_min_sd = r$min_sd,
             warm_status = diag_i$fit_status,
             warm_usable = usable_i,
             warm_iter = diag_i$n_iter,
             warm_min_sd = diag_i$min_sd,
             warm_degenerate = diag_i$degenerate,
             runtime_minutes = runtime_i)
})

warm_summary <- bind_rows(warm_results)

write.csv(warm_summary,
          file.path(misspec_diagnostic_table_dir,
                    "warm_start_comparison.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Diagnostic 4: RS prior sensitivity. Unif(1, 200), wider lognormal, 
# centred lognormal
#-------------------------------------------------------------------------------

# Priors considered: wider lognormal, recentred log normal, uniform.
ta_priors <- list(
  wide_lognormal = list(type = "lognormal", meanlog = log(200), sdlog = 2),
  mse_centered = list(type = "lognormal", meanlog = log(ta_best), sdlog = 0.5),
  uniform_1_200 = list(type = "uniform", min = 1, max = 200))

# New links to latent Gaussian scale
make_rs_ta_links <- function(pr) {
  link <- make_links_P0("rate_state")
  
  if (pr$type == "lognormal") {
    link$ta <- \(x) loggaus_t(x, pr$meanlog, pr$sdlog)
    inv_ta <- \(x) inv_loggaus_t(x, pr$meanlog, pr$sdlog)
  }
  
  if (pr$type == "uniform") {
    link$ta <- \(x) unif_t(x, pr$min, pr$max)
    inv_ta <- \(x) inv_unif_t(x, pr$min, pr$max)
  }
  
  list(link = link, inv_ta = inv_ta)
}

prior_results <- list()

# Refit with new prior specifications
for (prior_name in names(ta_priors)) {
  pr <- ta_priors[[prior_name]]
  transforms <- make_rs_ta_links(pr)
  
  for (i in seq_len(nrow(test_rows))) {
    r <- test_rows[i, ]
    
    sim_i <- simulation_manifest %>%
      filter(form == "mse", rep == r$rep)
    
    catalogue_i <- readRDS(file.path(catalogue_dir, sim_i$file))$catalogue
    
    catalogue_i <- prepare_temporal_catalogue(catalogue_i,
                                              T1 = T_fit_start,
                                              T2 = T_fit_end)
    
    set.seed(r$fit_seed)
    
    link_i <- transforms$link
    bru_i <- make_bru_options_P0("rate_state", rel_tol = fit_control$rel_tol, 
                                 max_iter = fit_control$max_iter)
    bru_i$bru_initial$th.ta <- transforms$inv_ta(initials[["rate_state"]]["ta"])
    bru_i$control.compute <- list(config = TRUE)
    start_time <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue_i,
                                         M0 = M0, 
                                         T1 = T_fit_start,
                                         T2 = T_fit_end, 
                                         link.functions = link_i, 
                                         coef.t. = temporal_binning$coef.t, 
                                         delta.t. = temporal_binning$delta.t, 
                                         N.max. = temporal_binning$N.max, 
                                         bru.opt = bru_i,
                                         form = "rate_state")
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    d <- fit_diagnostics(fit_i)
    usable_i <- fit_is_usable(d)
    
    prior_results[[paste(prior_name, r$rep, sep = "_")]] <- data.frame(
      prior = prior_name,
      rep = r$rep,
      fit_status = d$fit_status,
      usable = usable_i,
      converged = d$converged,
      n_iter = d$n_iter,
      min_sd = d$min_sd,
      degenerate = d$degenerate,
      runtime_minutes = runtime_i)
    saveRDS(fit_i,
            file.path(misspec_diagnostic_fit_dir,
                      sprintf("rs_ta_%s_rep_%04d.rds", prior_name, r$rep)))
  }
}

prior_summary <- bind_rows(prior_results)

prior_summary_by_prior <- prior_summary %>%
  group_by(prior) %>%
  summarise(n = n(),
            usable = sum(usable),
            converged = sum(converged),
            degenerate = sum(degenerate),
            median_iter = median(n_iter),
            .groups = "drop")

write.csv(prior_summary,
          file.path(misspec_diagnostic_table_dir,
                    "prior_sensitivity.csv"),
          row.names = FALSE)

write.csv(prior_summary_by_prior,
          file.path(misspec_diagnostic_table_dir,
                    "prior_sensitivity_summary.csv"),
          row.names = FALSE)

message("Saved misspecification diagnostic fits to: ",
        misspec_diagnostic_fit_dir)
message("Saved misspecification diagnostic tables to: ",
        misspec_diagnostic_table_dir)
message("Saved misspecification diagnostic figures to: ",
        misspec_diagnostic_figure_dir)
message("Finished misspecification diagnostics.")