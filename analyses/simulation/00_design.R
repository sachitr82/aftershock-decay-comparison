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
b_true <- 0.783
beta_true <- b_true * log(10)

#-------------------------------------------------------------------------------
# Observation times
#-------------------------------------------------------------------------------

mainshock_date <- as.Date("2019-07-06")

fit_start_date <- as.Date("2016-01-01")

# Boundary: events before this date are used for fitting
fit_end_date <- as.Date("2022-01-06")

# Boundary: forecast period ends immediately before this date
forecast_end_date <- as.Date("2024-01-06")

# Convert dates to days relative to mainshock

T_fit_start <- as.numeric(fit_start_date - mainshock_date)

T_fit_end <- as.numeric(fit_end_date - mainshock_date)

T_forecast_end <- as.numeric(forecast_end_date - mainshock_date)

#-------------------------------------------------------------------------------
# Imposed Ridgecrest-like mainshock
#-------------------------------------------------------------------------------

mainshock_event <- data.frame(ts = 0, magnitudes = 7.1)

#-------------------------------------------------------------------------------
# Common parameters across models
#-------------------------------------------------------------------------------

# Empirically Ridgecrest pre-foreshock total event rate used as a simulation
# background-rate anchor. Not a direct estimate of mu.
mu_true <- 0.075

# Ridgecrest-informed anchor
alpha_true <- 1.89

#-------------------------------------------------------------------------------
# Temporal kernel shape parameters
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
# Common mainshock triggering mass before the fitting boundary [0,915)
# Enables controlled comparison: given the same imposed mainshock, choose
# parameters to return the same expected direct triggering before the fitting 
# boundary, while changing shape in time of that triggering
#-------------------------------------------------------------------------------

G_mainshock_obs_target <- ETAS.inlabru::temporal_kernel_integral(
  a = 0, b = T_fit_end, theta = list(c = c_true, p = p_true), kernel = "ou")

rho_true <- uniroot(
  function(rho)  {
    ETAS.inlabru::temporal_kernel_integral(
      a = 0, b = T_fit_end,
      theta = list(d = d_true, rho = rho,gamma = gamma_true),kernel = "mse") - 
      G_mainshock_obs_target
  },
  interval = c(1e-8, 100), tol = 1e-12)$root

B_true <- uniroot(
  function(B)  {
    ETAS.inlabru::temporal_kernel_integral(
      a = 0, b = T_fit_end,
      theta = list(B = B,ta = ta_true),kernel = "rate_state") - 
      G_mainshock_obs_target
  },
  interval = c(1e-8, 1 - 1e-10), tol = 1e-12)$root

G_mainshock_obs <- c(
  ou = ETAS.inlabru::temporal_kernel_integral(
    a = 0, b = T_fit_end, theta = list(c = c_true, p = p_true),kernel = "ou"),
  
  mse = ETAS.inlabru::temporal_kernel_integral(
    a = 0, b = T_fit_end,
    theta = list(d = d_true, rho = rho_true ,gamma = gamma_true),kernel = "mse"),
  
  rate_state = ETAS.inlabru::temporal_kernel_integral(
    a = 0, b = T_fit_end, theta = list(B = B_true ,ta = ta_true),
    kernel = "rate_state"))

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
# Lifetime temporal kernel masses
#-------------------------------------------------------------------------------

G_inf <- c(
  
  ou = ETAS.inlabru::temporal_kernel_total_mass(
    theta = list(c = c_true, p = p_true),
    kernel = "ou"),
  
  mse = ETAS.inlabru::temporal_kernel_total_mass(
    theta = list(d = d_true, rho = rho_true,gamma = gamma_true), kernel = "mse"),
  
  rate_state = ETAS.inlabru::temporal_kernel_total_mass(
    theta = list(B = B_true, ta = ta_true), kernel = "rate_state")
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

#-------------------------------------------------------------------------------
# Pilot seeds
#-------------------------------------------------------------------------------

pilot_seeds <- c(ou = 900001, mse = 900002, rate_state = 900003)

#-------------------------------------------------------------------------------
# Deterministic seeds
#-------------------------------------------------------------------------------

sim_index <- rbind(
  
  data.frame(kernel = "ou", rep = 1:n_rep, seed = 1001:(1000 + n_rep)),
  
  data.frame(kernel = "mse", rep = 1:n_rep,seed = 2001:(2000 + n_rep)),
  
  data.frame(kernel = "rate_state", rep = 1:n_rep,seed = 3001:(3000 + n_rep))
)

#-------------------------------------------------------------------------------
# Output directories
#-------------------------------------------------------------------------------

pilot_dir <- here::here("results", "simulation", "pilots")
catalogue_dir <- here::here("results", "simulation", "catalogues")
figure_dir <- here::here("results", "simulation", "figures")
table_dir <- here::here("results", "simulation", "tables")
fit_dir <- here::here("results", "simulation", "fits")
forecast_fit_dir <- here::here("results", "simulation", "forecast_fits")

dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalogue_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(forecast_fit_dir, recursive = TRUE, showWarnings = FALSE)