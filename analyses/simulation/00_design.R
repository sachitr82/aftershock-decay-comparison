#===============================================================================
# Synthetic study design
#===============================================================================

library(here)
library(ETAS.inlabru)

#-------------------------------------------------------------------------------
# Number of replicates
#-------------------------------------------------------------------------------

n_rep <- 100

#-------------------------------------------------------------------------------
# Magnitudes: cutoff and upper bound for randomly simulated magnitudes
#-------------------------------------------------------------------------------

M0 <- 2.5
Mmax <- 7.1

# Ridgecrest b-value (obtained from EDA)
b_true <- 0.793967
beta_true <- b_true * log(10)

#-------------------------------------------------------------------------------
# Observation times
#-------------------------------------------------------------------------------

mainshock_date <- as.Date("2019-07-06")

fit_start_date <- as.Date("2016-01-01")

# Convert dates to days relative to mainshock

T_fit_start <- as.numeric(fit_start_date - mainshock_date)

# Full Ridgecrest post-mainshock observation window
T_fit_end <- 2371

#-------------------------------------------------------------------------------
# Imposed Ridgecrest-like mainshock
#-------------------------------------------------------------------------------

mainshock_event <- data.frame(ts = 0, magnitudes = 7.1)

#-------------------------------------------------------------------------------
# Common parameters across models
#-------------------------------------------------------------------------------

# Empirically Ridgecrest pre-foreshock total event rate used as a simulation
# background-rate anchor. Not a direct estimate of mu.
mu_true <- 0.074

# Ridgecrest-informed anchor
alpha_true <- 1.89

#-------------------------------------------------------------------------------
# Temporal decay shape parameters
#-------------------------------------------------------------------------------

# Omori--Utsu (Ridgecrest-informed)
c_true <- 0.03
p_true <- 1.17

# Modified stretched exponential (d matching Omori c, gamma backed by literature)
d_true <- 0.03
gamma_true <- 0.22

# Rate-state (backed by literature)
ta_true <- 188

#-------------------------------------------------------------------------------
# Common mainshock triggering mass before the fitting boundary [0,T_fit_end)
# Enables controlled comparison: given the same imposed mainshock, choose
# parameters to return the same expected direct triggering before the fitting 
# boundary, while changing shape in time of that triggering
#-------------------------------------------------------------------------------

G_mainshock_obs_target <- ETAS.inlabru::temporal_decay_integral(
  a = 0, b = T_fit_end, theta = list(c = c_true, p = p_true), form = "ou")

rho_true <- uniroot(
  function(rho)  {
    ETAS.inlabru::temporal_decay_integral(
      a = 0, b = T_fit_end,
      theta = list(d = d_true, rho = rho,gamma = gamma_true), form = "mse") - 
      G_mainshock_obs_target
  },
  interval = c(1e-8, 100), tol = 1e-12)$root

B_true <- uniroot(
  function(B)  {
    ETAS.inlabru::temporal_decay_integral(
      a = 0, b = T_fit_end,
      theta = list(B = B,ta = ta_true), form = "rate_state") - 
      G_mainshock_obs_target
  },
  interval = c(1e-8, 1 - 1e-10), tol = 1e-12)$root

G_mainshock_obs <- c(
  ou = ETAS.inlabru::temporal_decay_integral(
    a = 0, b = T_fit_end, theta = list(c = c_true, p = p_true), form = "ou"),
  
  mse = ETAS.inlabru::temporal_decay_integral(
    a = 0, b = T_fit_end,
    theta = list(d = d_true, rho = rho_true ,gamma = gamma_true), form = "mse"),
  
  rate_state = ETAS.inlabru::temporal_decay_integral(
    a = 0, b = T_fit_end, theta = list(B = B_true ,ta = ta_true),
    form = "rate_state"))

stopifnot(max(abs(G_mainshock_obs - G_mainshock_obs_target)) < 1e-6)

#-------------------------------------------------------------------------------
# Mean magnitude-productivity multiplier under truncated GR
#-------------------------------------------------------------------------------

L_mag <- Mmax - M0

if (abs(alpha_true - beta_true) < 1e-10) {
  mean_productivity_multiplier <- beta_true * L_mag /(1 - exp(-beta_true * L_mag))
} else{
mean_productivity_multiplier <- 
  beta_true / (beta_true - alpha_true) * 
  (1 - exp(-(beta_true - alpha_true) * L_mag)) / (1 - exp(-beta_true * L_mag))
}

#-------------------------------------------------------------------------------
# Lifetime temporal decay masses
#-------------------------------------------------------------------------------

G_inf <- c(
  
  ou = ETAS.inlabru::temporal_decay_total_mass(
    theta = list(c = c_true, p = p_true),
    form = "ou"),
  
  mse = ETAS.inlabru::temporal_decay_total_mass(
    theta = list(d = d_true, rho = rho_true,gamma = gamma_true), form = "mse"),
  
  rate_state = ETAS.inlabru::temporal_decay_total_mass(
    theta = list(B = B_true, ta = ta_true), form = "rate_state")
)

#-------------------------------------------------------------------------------
# Common productivity scale 
#-------------------------------------------------------------------------------

# Calibration-informed synthetic design choice, fixed after generative validation
n_ou_anchor <- 0.90
K_true <- n_ou_anchor /(mean_productivity_multiplier * G_inf[["ou"]])

#-------------------------------------------------------------------------------
# Actual lifetime branching ratios
#-------------------------------------------------------------------------------

branching_ratio <- K_true * mean_productivity_multiplier * G_inf

stopifnot(all(branching_ratio < 1))

mainshock_direct_obs <- K_true * exp(alpha_true * (7.1 - M0)) * G_mainshock_obs

stopifnot(max(abs(mainshock_direct_obs - mainshock_direct_obs[1])) < 1e-6)

#------------------------------------------------------------------------------
# True parameter vectors
#------------------------------------------------------------------------------

truths <- list(
  
  ou = list(mu = mu_true, K = K_true, alpha = alpha_true, c = c_true, p = p_true),
  
  mse = list(mu = mu_true, K = K_true, alpha = alpha_true, d = d_true, 
             rho = rho_true,gamma = gamma_true),
  
  rate_state = list(mu = mu_true, K = K_true, alpha = alpha_true, B = B_true,
                    ta = ta_true)
)

#------------------------------------------------------------------------------
# Fixed fitting initial values
#------------------------------------------------------------------------------

initials <- list(
  
  ou = c(mu = 0.3, K = 0.1, alpha = 1, c = 0.2, p = 1.1),
  
  mse = c(mu = 0.3, K = 0.1, alpha = 1, d = 0.2, rho = 1, gamma = 0.5),
  
  rate_state = c(mu = 0.3, K = 0.1, alpha = 1, B = 100 / 100.2, ta = 100))

#------------------------------------------------------------------------------
# Baseline fitting priors
#------------------------------------------------------------------------------

bound_eps <- 1e-6
prior_baseline <- list(
  mu = list(dist = "gamma", shape = 0.5, rate = 0.5),
  K = list(dist = "lognormal", meanlog = -1, sdlog = 0.5),
  alpha = list(dist = "gamma", shape = 1, rate = 0.5),
  
  ou = list(
    c = list(dist = "uniform", min = bound_eps, max = 1),
    p = list(dist = "uniform", min = 1 + bound_eps, max = 2)
  ),
  
  mse = list(
    d = list(dist = "uniform", min = bound_eps, max = 1),
    rho = list(dist = "lognormal", meanlog = log(1.5), sdlog = 1),
    gamma = list(dist = "uniform", min = bound_eps, max = 1 - bound_eps)
  ),
  
  rate_state = list(
    B = list(dist = "logit_normal", mean = 7.5, sd = 1),
    ta = list(dist = "lognormal", meanlog = log(200), sdlog = 0.5)
  )
)

prior_calibration <- list(n_draws = 10000, T = T_fit_end,
                          seed = 800001)

#-------------------------------------------------------------------------------
# Fixed baseline priors
#-------------------------------------------------------------------------------

# Forward copula transformations:
# internal N(0,1) scale -> physical ETAS parameter scale

make_links_P0 <- function(form) {
  
  stopifnot(form %in% c("ou", "mse", "rate_state"))
  
  common <- list(
    mu = \(x) gamma_t(x, prior_baseline$mu$shape,prior_baseline$mu$rate),
    K = \(x) loggaus_t(x, prior_baseline$K$meanlog, prior_baseline$K$sdlog),
    alpha = \(x) gamma_t(x, prior_baseline$alpha$shape, prior_baseline$alpha$rate))
  
  if (form == "ou") {
    return(c(common,list(
      c_ = \(x) unif_t(x, prior_baseline$ou$c$min, prior_baseline$ou$c$max),
      p = \(x) unif_t(x, prior_baseline$ou$p$min, prior_baseline$ou$p$max))))
  }
  
  if (form == "mse") {
    return(c(common, list(
      d = \(x) unif_t(x, prior_baseline$mse$d$min, prior_baseline$mse$d$max),
      rho = \(x) loggaus_t(x, prior_baseline$mse$rho$meanlog, 
                           prior_baseline$mse$rho$sdlog),
      gamma = \(x) unif_t(x, prior_baseline$mse$gamma$min, 
                          prior_baseline$mse$gamma$max))))
  }
  
  if (form == "rate_state") {
    return(c(common, list(
      B = \(x) logitgaus_t(x, prior_baseline$rate_state$B$mean, 
                           prior_baseline$rate_state$B$sd),
      ta = \(x) loggaus_t(x, prior_baseline$rate_state$ta$meanlog, 
                          prior_baseline$rate_state$ta$sdlog))))
  }
}

#-------------------------------------------------------------------------------
# Inverse copula transformations:
# physical ETAS parameter scale -> internal N(0,1) scale
#-------------------------------------------------------------------------------

make_inverse_links_P0 <- function(form) {
  
  stopifnot(form %in% c("ou", "mse", "rate_state"))
  
  common <- list(
    mu = \(x) inv_gamma_t(x, prior_baseline$mu$shape, prior_baseline$mu$rate),
    
    K = \(x) inv_loggaus_t(x, prior_baseline$K$meanlog, prior_baseline$K$sdlog),
    
    alpha = \(x) inv_gamma_t(x, prior_baseline$alpha$shape, prior_baseline$alpha$rate))
  
  
  if (form == "ou") {
    return(c(common, list(
      c_ = \(x) inv_unif_t(x,  prior_baseline$ou$c$min, prior_baseline$ou$c$max),
      p = \(x) inv_unif_t(x, prior_baseline$ou$p$min, prior_baseline$ou$p$max))))
  }
  
  
  if (form == "mse") {
    return(c(common, list( 
      d = \(x) inv_unif_t(x,  prior_baseline$mse$d$min, prior_baseline$mse$d$max),
      rho = \(x) inv_loggaus_t(x, prior_baseline$mse$rho$meanlog,
                               prior_baseline$mse$rho$sdlog),
      gamma = \(x) inv_unif_t(x, prior_baseline$mse$gamma$min, 
                              prior_baseline$mse$gamma$max))))
  }
  
  if (form == "rate_state") { 
    return(c(common, list(
      B = \(x) inv_logitgaus_t(x, prior_baseline$rate_state$B$mean,
                               prior_baseline$rate_state$B$sd),
      ta = \(x) inv_loggaus_t(x, prior_baseline$rate_state$ta$meanlog,
                              prior_baseline$rate_state$ta$sdlog))))
  }
}

#-------------------------------------------------------------------------------
# inlabru fitting options under baseline priors
#-------------------------------------------------------------------------------

make_bru_options_P0 <- function(form, rel_tol = 0.1, max_iter = 100) {
  
  stopifnot(form %in% c("ou", "mse", "rate_state"))
  
  inv <- make_inverse_links_P0(form)
  init <- initials[[form]]
  
  if (form == "ou") {
    th_init <- list(
      th.mu = inv$mu(init["mu"]),
      th.K = inv$K(init["K"]),
      th.alpha = inv$alpha(init["alpha"]),
      th.c = inv$c_(init["c"]),
      th.p = inv$p(init["p"]))
  }
  
  if (form == "mse") {
    th_init <- list(
      th.mu = inv$mu(init["mu"]),
      th.K = inv$K(init["K"]),
      th.alpha = inv$alpha(init["alpha"]),
      th.d = inv$d(init["d"]),
      th.rho = inv$rho(init["rho"]),
      th.gamma = inv$gamma(init["gamma"]))
  }
  
  if (form == "rate_state") {
    th_init <- list(
      th.mu = inv$mu(init["mu"]),
      th.K = inv$K(init["K"]),
      th.alpha = inv$alpha(init["alpha"]),
      th.B = inv$B(init["B"]),
      th.ta = inv$ta(init["ta"]))
  }
  
  list(bru_verbose = 0, bru_max_iter = max_iter, 
       bru_method = list(rel_tol = rel_tol),
       bru_initial = th_init)
}

#-------------------------------------------------------------------------------
# Fitting numerical controls (validated in 04 and 05)
#-------------------------------------------------------------------------------

temporal_binning <- list(coef.t = 1, delta.t = 0.1, N.max = 14)

fit_control <- list(rel_tol = 0.1, max_iter = 100)

#-------------------------------------------------------------------------------
# Pilot seeds
#-------------------------------------------------------------------------------

pilot_seeds <- c(ou = 900001, mse = 900002, rate_state = 900003)

#-------------------------------------------------------------------------------
# Deterministic seeds
#-------------------------------------------------------------------------------

sim_index <- rbind(
  
  data.frame(form = "ou", rep = 1:n_rep, seed = 1001:(1000 + n_rep)),
  
  data.frame(form = "mse", rep = 1:n_rep,seed = 2001:(2000 + n_rep)),
  
  data.frame(form = "rate_state", rep = 1:n_rep,seed = 3001:(3000 + n_rep))
)

#-------------------------------------------------------------------------------
# Output directories
#-------------------------------------------------------------------------------

simulation_outputs_dir <- here("outputs", "simulation")

pilot_dir <- file.path(simulation_outputs_dir, "pilots")
catalogue_dir <- file.path(simulation_outputs_dir, "catalogues")
figure_dir <- file.path(simulation_outputs_dir, "figures")
table_dir <- file.path(simulation_outputs_dir, "tables")
fit_dir <- file.path(simulation_outputs_dir, "fits")

dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalogue_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)