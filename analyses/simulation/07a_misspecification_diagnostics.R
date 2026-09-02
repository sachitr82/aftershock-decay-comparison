#===============================================================================
# Fit convergence diagnostics
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

#-------------------------------------------------------------------------------
# Directories and fit files
#-------------------------------------------------------------------------------

main_fit_dir <- file.path(fit_dir, "main_simulation")
diag_dir <- file.path(main_fit_dir, "initialisation_diagnostic")
prior_diag_dir <- file.path(main_fit_dir, "prior_diagnostic")

dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prior_diag_dir, recursive = TRUE, showWarnings = FALSE)

fm <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"), stringsAsFactors = FALSE)
simulation_manifest <- read.csv(file.path(catalogue_dir, "simulation_manifest.csv"),
                                stringsAsFactors = FALSE)
fit_index <- read.csv(file.path(main_fit_dir, "fit_index.csv"), stringsAsFactors = FALSE)

fm$fit_file <- file.path(main_fit_dir, paste0("truth_", fm$truth_kernel), fm$file)

rs_link <- make_links_P0("rate_state")

#-------------------------------------------------------------------------------
# Diagnostic 1: terminal rate-state parameter geometry
# RS fits to OU-, MSE-, and RS-generated catalogues
#-------------------------------------------------------------------------------

rs_rows <- fm[fm$fitted_kernel == "rate_state", ]

rs_diag <- do.call(rbind, lapply(seq_len(nrow(rs_rows)), function(i) {
  r <- rs_rows[i, ]
  o <- readRDS(r$fit_file)
  B <- rs_link$B(o$fit$summary.fixed["th.B", "mean"])
  ta <- rs_link$ta(o$fit$summary.fixed["th.ta", "mean"])
  tau <- (1 - B) * ta / B
  data.frame(truth_kernel = r$truth_kernel, rep = r$rep, n_fit = r$n_fit, 
             converged = r$converged, hit_max = r$hit_max, 
             inla_failure = r$inla_failure, degenerate = r$degenerate, 
             min_sd = r$min_sd, B = B, one_minus_B = 1 - B, ta = ta, tau = tau)
}))

rs_diag$truth_kernel <- factor(rs_diag$truth_kernel, 
                               levels = c("ou", "mse", "rate_state"), 
                               labels = c("OU truth", "MSE truth", "RS truth"))
rs_diag$fit_result <- ifelse(rs_diag$converged & !rs_diag$degenerate, "Converged",
                             ifelse(rs_diag$degenerate, "Degenerate", "Non-converged"))
rs_diag$fit_result <- factor(rs_diag$fit_result, levels = c("Converged", 
                                                "Non-converged", "Degenerate"))

write.csv(rs_diag, file.path(main_fit_dir, "rs_parameter_geometry.csv"), 
          row.names = FALSE)

cat("\nCounts\n")
print(with(rs_diag, table(truth_kernel, fit_result)))

cat("\nParameter summary - RS fits\n")
print(rs_diag %>% 
        group_by(truth_kernel, fit_result) %>% 
        summarise(n = n(), median_B = median(B, na.rm = TRUE), 
                  median_1mB = median(one_minus_B, na.rm = TRUE), 
                  median_ta = median(ta, na.rm = TRUE), 
                  median_tau = median(tau, na.rm = TRUE), 
                  q25_ta = quantile(ta, 0.25, na.rm = TRUE), 
                  q75_ta = quantile(ta, 0.75, na.rm = TRUE), 
                  q25_tau = quantile(tau, 0.25, na.rm = TRUE), 
                  q75_tau = quantile(tau, 0.75, na.rm = TRUE), .groups = "drop"))

cat("\nNon-degenerate RS fits only\n")
print(rs_diag %>% 
        filter(!degenerate) %>% 
        group_by(truth_kernel, fit_result) %>%
        summarise(n = n(), median_B = median(B, na.rm = TRUE), 
                  median_1mB = median(one_minus_B, na.rm = TRUE), 
                  median_ta = median(ta, na.rm = TRUE), 
                  median_tau = median(tau, na.rm = TRUE), 
                  q25_ta = quantile(ta, 0.25, na.rm = TRUE), 
                  q75_ta = quantile(ta, 0.75, na.rm = TRUE), 
                  q25_tau = quantile(tau, 0.25, na.rm = TRUE), 
                  q75_tau = quantile(tau, 0.75, na.rm = TRUE), .groups = "drop"))

p_ta <- ggplot(rs_diag %>% filter(is.finite(ta), ta > 0), 
               aes(truth_kernel, ta, shape = fit_result)) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 2) +
  scale_y_log10() +
  labs(x = NULL, y = expression(t[a]~"(days)"), title = expression("Terminal "*t[a])) +
  theme_bw() +
  theme(legend.position = "bottom")

p_B <- ggplot(rs_diag %>% filter(is.finite(one_minus_B), one_minus_B > 0), 
              aes(truth_kernel, one_minus_B, shape = fit_result)) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 2) +
  scale_y_log10() +
  labs(x = NULL, y = expression(1-B), title = expression("Terminal "*(1-B))) +
  theme_bw() +
  theme(legend.position = "bottom")

p_tau <- ggplot(rs_diag %>% filter(is.finite(tau), tau > 0), 
                aes(truth_kernel, tau, shape = fit_result)) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 2) +
  scale_y_log10() +
  labs(x = NULL, y = expression(tau == (1-B)*t[a]/B), 
       title = expression("Implied short-time scale "*tau)) +
  theme_bw() +
  theme(legend.position = "bottom")

p_rs_geometry <- (p_ta | p_B | p_tau) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(p_rs_geometry)
ggsave(file.path(main_fit_dir, "rs_parameter_geometry.pdf"), 
       p_rs_geometry, width = 11, height = 4)
#-------------------------------------------------------------------------------
# Diagnostic 2: Can RS resemble MSE?
#-------------------------------------------------------------------------------

T <- T_fit_end
tg <- c(0, 10^seq(-5, log10(T), length.out = 2000))

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

norm_shape <- function(theta, kern) {
  g <- temporal_kernel(tg, theta, kern)
  G <- temporal_kernel_integral(0, T, theta, kern)
  g / G
}

h_mse_true <- norm_shape(list(d = d_true, rho = rho_true, gamma = gamma_true), "mse")
h_rs_true <- norm_shape(list(B = B_true, ta = ta_true), "rate_state")

tv <- function(h1, h2) 0.5 * trapz(tg, abs(h1 - h2))

obj_rs <- function(z) {
  B <- plogis(z[1]); ta <- exp(z[2])
  h <- norm_shape(list(B = B, ta = ta), "rate_state")
  tv(h_mse_true, h)
}

rs_starts <- rbind(c(qlogis(B_true), log(ta_true)), c(qlogis(0.999), log(50)), 
                   c(qlogis(0.9999), log(100)), c(qlogis(0.99999),log(500)),
                   c(qlogis(0.999999), log(2000)))
rs_opts <- lapply(seq_len(nrow(rs_starts)), function(i) 
  optim(rs_starts[i, ], obj_rs, method = "Nelder-Mead", control = list(maxit = 5000)))
rs_best <- rs_opts[[which.min(vapply(rs_opts, function(x) x$value, numeric(1)))]]

B_best <- plogis(rs_best$par[1])
ta_best <- exp(rs_best$par[2])
tau_best <- (1 - B_best) * ta_best / B_best
tv_rs_to_mse <- rs_best$value

obj_mse <- function(z) {
  d <- exp(z[1]); rho <- exp(z[2]); gamma <- plogis(z[3])
  h <- norm_shape(list(d = d, rho = rho, gamma = gamma), "mse")
  tv(h_rs_true, h)
}

mse_starts <- rbind(c(log(d_true), log(rho_true), qlogis(gamma_true)), 
                    c(log(0.03), log(1), qlogis(0.2)), 
                    c(log(0.1), log(2), qlogis(0.4)), 
                    c(log(0.01), log(0.5), qlogis(0.1)), 
                    c(log(0.3), log(3), qlogis(0.7)))
mse_opts <- lapply(seq_len(nrow(mse_starts)), function(i) 
  optim(mse_starts[i, ], obj_mse, method = "Nelder-Mead", control = list(maxit = 5000)))
mse_best <- mse_opts[[which.min(vapply(mse_opts, function(x) x$value, numeric(1)))]]

d_best <- exp(mse_best$par[1])
rho_best <- exp(mse_best$par[2])
gamma_best <- plogis(mse_best$par[3])
tv_mse_to_rs <- mse_best$value

cat("\nBest RS approximation to MSE truth\n")
cat(sprintf("TV = %.6f | B = %.10f | 1-B = %.3g | ta = %.3f | tau = %.6f\n",
            tv_rs_to_mse, B_best, 1 - B_best, ta_best, tau_best))
cat(sprintf("RS prior positions: logit(B) z = %.2f | log(ta) z = %.2f\n",
            (qlogis(B_best) - 7.5) / 1, (log(ta_best) - log(200)) / 0.5))

cat("\nBest MSE approximation to RS truth\n")
cat(sprintf("TV = %.6f | d = %.6f | rho = %.6f | gamma = %.6f\n",
            tv_mse_to_rs, d_best, rho_best, gamma_best))

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

ggplot(plot_dat, aes(t, h, linetype = form)) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~comparison) +
  labs(x = "Time since triggering event (days)", 
       y = "Normalised temporal triggering", linetype = NULL) +
  theme_bw()

#-------------------------------------------------------------------------------
# Diagnostic 2b: OU-to-RS control
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

B_ou_best <- plogis(rs_ou_best$par[1])
ta_ou_best <- exp(rs_ou_best$par[2])
tau_ou_best <- (1 - B_ou_best) * ta_ou_best / B_ou_best
tv_rs_to_ou <- rs_ou_best$value

cat("\nBest approximation to OU truth\n")
cat(sprintf("TV = %.6f | B = %.10f | 1-B = %.3g | ta = %.3f | tau = %.6f\n",
            tv_rs_to_ou, B_ou_best, 1 - B_ou_best, ta_ou_best, tau_ou_best))
cat(sprintf("RS prior positions: logit(B) z = %.2f | log(ta) z = %.2f\n",
            (qlogis(B_ou_best) - 7.5), (log(ta_ou_best) - log(200)) / 0.5))

cat("\nComparison with MSE truth\n")
cat(sprintf("RS -> OU  | TV = %.6f | ta = %.3f | tau = %.6f | ta prior z = %.2f\n",
            tv_rs_to_ou, ta_ou_best, tau_ou_best, (log(ta_ou_best) - log(200)) / 0.5))
cat(sprintf("RS -> MSE | TV = %.6f | ta = %.3f | tau = %.6f | ta prior z = %.2f\n",
            tv_rs_to_mse, ta_best, tau_best, (log(ta_best) - log(200)) / 0.5))

h_rs_ou_best <- norm_shape(list(B = B_ou_best, ta = ta_ou_best), "rate_state")

plot_ou_rs <- rbind(data.frame(t = tg[-1], h = h_ou_true[-1], form = "OU truth"),
                    data.frame(t = tg[-1], h = h_rs_ou_best[-1], form = "Best RS approximation"))

ggplot(plot_ou_rs, aes(t, h, linetype = form)) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "OU truth -> RS approximation",
       x = "Time since triggering event (days)",
       y = "Normalised temporal triggering", linetype = NULL) +
  theme_bw()

#-------------------------------------------------------------------------------
# Diagnostic 3: Initialisation sensitivity
#-------------------------------------------------------------------------------

test_rows <- fm %>%
  filter(truth_kernel == "mse", fitted_kernel == "rate_state", 
         fit_status == "hit_max", !degenerate) %>%
  arrange(runtime_minutes) %>%
  slice_head(n = 5)
print(test_rows[, c("rep", "n_fit", "n_iter", "min_sd")])

fit_diag <- function(fit) {
  log_i <- as.character(inlabru::bru_log(fit))
  
  inla_failure <- any(grepl(paste(
    "Problem in inla",
    "Giving up and returning last successfully obtained result",
    "inla-program exited with an error",
    "maximum number of tries has been reached",
    "Newton-Raphson optimizer did not converge",
    sep = "|"
  ), log_i))
  
  nan_inf_logl <- any(grepl("NAN/INF values in logl", log_i, fixed = TRUE))
  vb_aborted <- any(grepl("max_correction|vb.correction.*aborted", log_i))
  
  converged <- !inla_failure &&
    any(grepl("Convergence criterion met", log_i, fixed = TRUE))
  
  hit_max <- !inla_failure &&
    any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
  
  n_iter <- max(fit$bru_iinla$track$iteration, na.rm = TRUE)
  min_sd <- min(fit$summary.fixed$sd, na.rm = TRUE)
  degenerate <- min_sd < 1e-6
  
  status <- case_when(
    inla_failure ~ "inla_failure",
    converged ~ "converged",
    hit_max ~ "hit_max",
    TRUE ~ "unknown"
  )
  
  usable <- converged && !hit_max && !degenerate && !vb_aborted &&
    !inla_failure && !nan_inf_logl
  
  list(
    status = status,
    usable = usable,
    converged = converged,
    hit_max = hit_max,
    inla_failure = inla_failure,
    nan_inf_logl = nan_inf_logl,
    vb_aborted = vb_aborted,
    n_iter = n_iter,
    min_sd = min_sd,
    degenerate = degenerate
  )
}

warm_results <- lapply(seq_len(nrow(test_rows)), function(i) {
  r <- test_rows[i, ]
  
  sim_i <- simulation_manifest %>% filter(kernel == "mse", rep == r$rep)
  job_i <- fit_index %>% filter(truth_kernel == "mse", rep == r$rep, fitted_kernel == "rate_state")
  
  catalogue_i <- readRDS(file.path(catalogue_dir, sim_i$file))$catalogue
  catalogue_i <- catalogue_i %>% filter(ts >= T_fit_start, ts <= T_fit_end) %>% arrange(ts)
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  set.seed(job_i$fit_seed)
  
  link_i <- make_links_P0("rate_state")
  inv_i <- make_inverse_links_P0("rate_state")
  bru_i <- make_bru_options_P0("rate_state", rel_tol = fit_control$rel_tol, 
                               max_iter = fit_control$max_iter)
  
  bru_i$bru_initial$th.B <- inv_i$B(B_best)
  bru_i$bru_initial$th.ta <- inv_i$ta(ta_best)
  bru_i$control.compute <- list(config = TRUE, dic = TRUE, waic = TRUE, mlik = TRUE)
  
  start_time <- Sys.time()
  
  fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue_i, M0 = M0, 
                                       T1 = T_fit_start, T2 = T_fit_end, 
                                       link.functions = link_i, 
                                       coef.t. = temporal_binning$coef.t, 
                                       delta.t. = temporal_binning$delta.t, 
                                       N.max. = temporal_binning$N.max, 
                                       bru.opt = bru_i, kernel = "rate_state")
  
  runtime_i <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  diag_i <- fit_diag(fit_i)
  saveRDS(fit_i, file.path(diag_dir, sprintf("mse_rep_%04d_rs_warm_start.rds", 
                                             r$rep)))
  data.frame(
    rep = r$rep,
    n_fit = nrow(catalogue_i),
    baseline_min_sd = r$min_sd,
    warm_status = diag_i$status,
    warm_usable = diag_i$usable,
    warm_iter = diag_i$n_iter,
    warm_min_sd = diag_i$min_sd,
    warm_degenerate = diag_i$degenerate,
    runtime_minutes = runtime_i)
})

warm_summary <- bind_rows(warm_results)

print(warm_summary)
write.csv(warm_summary, file.path(diag_dir, "warm_start_comparison.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Diagnostic 4: RS prior sensitivity. Unif(1, 200), wider lognormal, 
# centred lognormal
#-------------------------------------------------------------------------------

ta_priors <- list(
  wide_lognormal = list(type = "lognormal", meanlog = log(200), sdlog = 2),
  mse_centered = list(type = "lognormal", meanlog = log(ta_best), sdlog = 0.5),
  uniform_1_200 = list(type = "uniform", min = 1, max = 200))

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

for (prior_name in names(ta_priors)) {
  pr <- ta_priors[[prior_name]]
  transforms <- make_rs_ta_links(pr)
  
  for (i in seq_len(nrow(test_rows))) {
    r <- test_rows[i, ]
    
    sim_i <- simulation_manifest %>% filter(kernel == "mse", rep == r$rep)
    job_i <- fit_index %>% filter(truth_kernel == "mse", rep == r$rep, 
                                  fitted_kernel == "rate_state")
    
    catalogue_i <- readRDS(file.path(catalogue_dir, sim_i$file))$catalogue
    catalogue_i <- catalogue_i %>% filter(ts >= T_fit_start, ts <= T_fit_end) %>% arrange(ts)
    catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
    
    set.seed(job_i$fit_seed)
    
    link_i <- transforms$link
    bru_i <- make_bru_options_P0("rate_state", rel_tol = fit_control$rel_tol, 
                                 max_iter = fit_control$max_iter)
    bru_i$bru_initial$th.ta <- transforms$inv_ta(initials[["rate_state"]]["ta"])
    bru_i$control.compute <- list(config = TRUE, dic = TRUE, waic = TRUE, mlik = TRUE)
    start_time <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue_i, M0 = M0, 
                                         T1 = T_fit_start, T2 = T_fit_end, 
                                         link.functions = link_i, 
                                         coef.t. = temporal_binning$coef.t, 
                                         delta.t. = temporal_binning$delta.t, 
                                         N.max. = temporal_binning$N.max, 
                                         bru.opt = bru_i, kernel = "rate_state")
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    d <- fit_diag(fit_i)
    
    prior_results[[paste(prior_name, r$rep, sep = "_")]] <- data.frame(
      prior = prior_name,
      rep = r$rep,
      status = d$status,
      usable = d$usable,
      converged = d$converged,
      n_iter = d$n_iter,
      min_sd = d$min_sd,
      degenerate = d$degenerate,
      runtime_minutes = runtime_i)
    saveRDS(fit_i, file.path(prior_diag_dir, sprintf("%s_rep_%04d.rds", prior_name, r$rep)))
  }
}

prior_summary <- bind_rows(prior_results)

print(prior_summary)

print(prior_summary %>%
        group_by(prior) %>%
        summarise(
          n = n(),
          usable = sum(usable),
          converged = sum(converged),
          degenerate = sum(degenerate),
          median_iter = median(n_iter),
          .groups = "drop"))
write.csv(prior_summary, file.path(prior_diag_dir, "prior_sensitivity.csv"), row.names = FALSE)