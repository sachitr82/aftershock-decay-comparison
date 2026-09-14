#===============================================================================
# Prior calibration / prior-predictive check
#===============================================================================

#-------------------------------------------------------------------------------
# Load package and fixed simulation design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(here)
library(ggplot2)
library(patchwork)

source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Fixed prior-calibration seed and observation window
#-------------------------------------------------------------------------------

set.seed(prior_calibration$seed)

n <- prior_calibration$n_draws
T_obs <- prior_calibration$T

#-------------------------------------------------------------------------------
# Draw temporal parameters from baseline priors
#-------------------------------------------------------------------------------

draws <- list(
  ou = data.frame(
    c = runif(n, prior_baseline$ou$c$min, prior_baseline$ou$c$max),
    p = runif(n, prior_baseline$ou$p$min, prior_baseline$ou$p$max)),
  
  mse = data.frame(
    d = runif(n, prior_baseline$mse$d$min, prior_baseline$mse$d$max),
    rho = rlnorm(n, prior_baseline$mse$rho$meanlog, prior_baseline$mse$rho$sdlog),
    gamma = runif(n, prior_baseline$mse$gamma$min, prior_baseline$mse$gamma$max)),
  
  rate_state = data.frame(
    B = plogis(rnorm(n, prior_baseline$rate_state$B$mean,
                     prior_baseline$rate_state$B$sd)),
    ta = rlnorm(n, prior_baseline$rate_state$ta$meanlog,
                prior_baseline$rate_state$ta$sdlog))
)

#-------------------------------------------------------------------------------
# Function to convert prior draw rows into temporal parameter lists
#-------------------------------------------------------------------------------

get_theta <- function(form, x, i) {
  
  if (form == "ou")
    return(list(c = x$c[i], p = x$p[i]))
  
  if (form == "mse")
    return(list(d = x$d[i], rho = x$rho[i], gamma = x$gamma[i]))
  
  if (form == "rate_state")
    return(list(B = x$B[i], ta = x$ta[i]))
  
  stop("Unknown temporal decay form.")
}

#-------------------------------------------------------------------------------
# Function to calculate finite-window temporal mass G(0, T_obs)
#-------------------------------------------------------------------------------

calc_mass <- function(form, x) {
  
  vapply(seq_len(nrow(x)), function(i) {
    ETAS.inlabru::temporal_decay_integral(
      a = 0, b = T_obs, theta = get_theta(form, x, i), form = form)
  }, numeric(1))
}

#-------------------------------------------------------------------------------
# Calculate finite-window temporal mass for all prior draws
#-------------------------------------------------------------------------------

draws$ou$G_0_T <- calc_mass("ou", draws$ou)
draws$mse$G_0_T <- calc_mass("mse", draws$mse)
draws$rate_state$G_0_T <- calc_mass("rate_state", draws$rate_state)

#-------------------------------------------------------------------------------
# Function to summarise prior distributions
#-------------------------------------------------------------------------------

q_summary <- function(x) {
  
  q <- quantile(x, c(.05, .5, .95), na.rm = TRUE)
  c(q05 = q[1], median = q[2], q95 = q[3])
}

#-------------------------------------------------------------------------------
# Summarise marginal priors and finite-window temporal mass
#-------------------------------------------------------------------------------

prior_summary <- rbind(
  ou_c = q_summary(draws$ou$c),
  ou_p = q_summary(draws$ou$p),
  ou_G_0_T = q_summary(draws$ou$G_0_T),
  
  mse_d = q_summary(draws$mse$d),
  mse_rho = q_summary(draws$mse$rho),
  mse_gamma = q_summary(draws$mse$gamma),
  mse_G_0_T = q_summary(draws$mse$G_0_T),
  
  rs_B = q_summary(draws$rate_state$B),
  rs_ta = q_summary(draws$rate_state$ta),
  rs_G_0_T = q_summary(draws$rate_state$G_0_T)
)

prior_summary <- data.frame(quantity = rownames(prior_summary), 
                            prior_summary, row.names = NULL)

print(prior_summary)

#-------------------------------------------------------------------------------
# Save prior summary table
#-------------------------------------------------------------------------------

tablefile <- file.path(table_dir, "baseline_prior_summary.csv")

write.csv(prior_summary, tablefile, row.names = FALSE)

message("Saved ", tablefile)

#-------------------------------------------------------------------------------
# Function to evaluate complete temporal decay functions
#-------------------------------------------------------------------------------

idx_curve <- seq_len(n)

# Log-spaced time grid captures short- and long-time behaviour
t_grid <- exp(seq(log(1e-4), log(T_obs), length.out = 500))

calc_curves <- function(form, x, idx) {
  
  vapply(idx, function(i) {
    ETAS.inlabru::temporal_decay(dt = t_grid,
                                  theta = get_theta(form, x, i),
                                  form = form)
  }, 
  numeric(length(t_grid)))
}

#-------------------------------------------------------------------------------
# Evaluate raw temporal decay functions g(t)
#-------------------------------------------------------------------------------

g_plot <- list(ou = calc_curves("ou", draws$ou, idx_curve),
               mse = calc_curves("mse", draws$mse, idx_curve),
               rate_state = calc_curves("rate_state", draws$rate_state, idx_curve))

#-------------------------------------------------------------------------------
# Normalise each decay function by its own fitting-window mass
# Integral_0^T h(t) dt = 1 for every prior draw
#-------------------------------------------------------------------------------

h_plot <- list(
  ou = sweep(g_plot$ou, 2, draws$ou$G_0_T, "/"),
  mse = sweep(g_plot$mse, 2, draws$mse$G_0_T, "/"),
  rate_state = sweep(g_plot$rate_state, 2, draws$rate_state$G_0_T, "/")
)

##-------------------------------------------------------------------------------
# Plot fitting-window-normalised temporal decay functions
#-------------------------------------------------------------------------------

plot_prior_draws <- function(H, title) {
  
  plot_data <- data.frame(time = rep(t_grid, ncol(H)), 
                          h = pmax(as.vector(H), 1e-300),
                          draw = rep(seq_len(ncol(H)), each = length(t_grid)))
  
  ggplot(plot_data, aes(x = time, y = h, group = draw)) +
  geom_line(alpha = 0.01, linewidth = 0.2) +
  scale_x_log10() + scale_y_log10() +
  coord_cartesian(ylim = c(1e-10, 1e2)) +
  labs (title = title, x = "Time since parent event (days)", y = "h(t)") +
  theme_bw()
}

plots <- list(
  ou = plot_prior_draws(h_plot$ou, "OU baseline prior"),
  mse = plot_prior_draws(h_plot$mse, "MSE baseline prior"),
  rate_state = plot_prior_draws(h_plot$rate_state, "Rate-state baseline prior"))

figfile <- file.path(figure_dir, "baseline_prior_normalised_kernel_draws.png")

ggsave(filename = figfile, plot = patchwork::wrap_plots(plots, ncol = 1),
       width = 9, height = 9)

message("Saved ", figfile)
message("Finished baseline prior calibration.")