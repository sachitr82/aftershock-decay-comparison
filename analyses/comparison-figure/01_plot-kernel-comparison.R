#===============================================================================
# Comparison of the three candidate triggering kernels
#
#   (a) normalised densities  f_k(t) = g_k(t) / G_k(0, Inf)      [log-log]
#   (b) survival functions    S_k(t) = 1 - F_k(t)                [log-x]
#
# Normalised forms are plotted to restrict kernel comparisons to shape of decay.  
# Survival functions display the fraction of an event's triggering still
# outstanding at time t.
#
# Parameters choices are informed following Hainzl & Christophersen (2017, Figure 1).
#===============================================================================

#-------------------------------------------------------------------------------
# Package dependencies
#-------------------------------------------------------------------------------

library(ggplot2)
library(here)
library(patchwork)

#-------------------------------------------------------------------------------
# Normalised densities and survival functions: log1p / expm1 for numerical stability
#-------------------------------------------------------------------------------

# Omori-Utsu 
f_ou <- function(t, c, p) (p - 1) / c * (1 + t / c)^(-p)
S_ou <- function(t, c, p) (1 + t / c)^(1 - p)

# Modified stretched exponential 
#  delta(t) = (d + t)^gamma - d^gamma, written to avoid cancellation at small t
delta_ms <- function(t, d, gamma) {
  d^gamma * expm1(gamma * log1p(t / d))
}

f_ms <- function(t, d, lambda, gamma) {
  lambda * gamma * (d + t)^(gamma - 1) * exp(-lambda * delta_ms(t, d, gamma))
}

S_ms <- function(t, d, lambda, gamma) {
  exp(-lambda * delta_ms(t, d, gamma))
}

# Rate-state 
f_rs <- function(t, B, ta) {
  z <- exp(-t / ta)
  B * z / (-ta * log1p(-B) * (1 - B * z))
}

S_rs <- function(t, B, ta) {
  log1p(-B * exp(-t / ta)) / log1p(-B)
}

#-------------------------------------------------------------------------------
# Parameters (Hainzl & Christophersen 2017, Figure 1)
#-------------------------------------------------------------------------------

par_ou <- list(c = 0.01,  p   = 1.1)
par_ms <- list(d = 0.01, lambda = 1, gamma = 0.2)
par_rs <- list(B = 0.9999, ta = 100)

kernel_levels <- c("Omori\u2013Utsu",
                   "Modified stretched exponential",
                   "Rate-state")

t_grid <- 10^seq(-4, 4, length.out = 3000)

#-------------------------------------------------------------------------------
# Evaluate kernels and order factor levels for legend
#-------------------------------------------------------------------------------

kern <- rbind(
  data.frame(t = t_grid, kernel = kernel_levels[1],
             density  = f_ou(t_grid, par_ou$c, par_ou$p),
             survival = S_ou(t_grid, par_ou$c, par_ou$p)),
  data.frame(t = t_grid, kernel = kernel_levels[2],
             density  = f_ms(t_grid, par_ms$d, par_ms$lambda, par_ms$gamma),
             survival = S_ms(t_grid, par_ms$d, par_ms$lambda, par_ms$gamma)),
  data.frame(t = t_grid, kernel = kernel_levels[3],
             density  = f_rs(t_grid, par_rs$B, par_rs$ta),
             survival = S_rs(t_grid, par_rs$B, par_rs$ta))
)

kern$kernel <- factor(kern$kernel, level = kernel_levels)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

# Kernel identifiers - colour and line type
kernel_cols <- setNames(c("steelblue", "orange", "springgreen4"), kernel_levels)
kernel_lty  <- setNames(c("solid", "dashed", "dotted"), kernel_levels)

# Format tick labels as, e.g., 10^2 instead 1e2
log_lab <- scales::trans_format("log10", scales::math_format(10^.x))

# Common theme choices for both plots
base_theme <- theme_bw(base_size = 10) +
  theme(legend.title     = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(0.9, "lines"),
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.spacing.x = unit(0.6, "cm"),
        legend.text      = element_text(size = 11),
        panel.grid.minor = element_line(linewidth = 0.15),
        plot.title       = element_text(size = 10, face = "plain"))

# Logarithmic x-axis used in both plots
x_scale <- scale_x_log10(breaks = 10^seq(-4, 4, by = 1), labels = log_lab,
                         expand = expansion(mult = 0))

#-------------------------------------------------------------------------------
# Panel (a): normalised densities (log-log)
#-------------------------------------------------------------------------------

p_density <- ggplot(kern, aes(t, density, colour = kernel, linetype = kernel)) +
  geom_line(linewidth = 1) +
  x_scale + 
  scale_y_log10(breaks = 10^seq(-12, 2, by = 2), labels = log_lab) +
  coord_cartesian(ylim = c(1e-12, 1e1)) +
  scale_colour_manual(values = kernel_cols) +
  scale_linetype_manual(values = kernel_lty) +
  labs(title = "(a) Normalised triggering density",
       x = "Time since triggering event (days)",
       y = expression(f[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Panel (b): survival functions (log x - linear y)
#-------------------------------------------------------------------------------

p_survival <- ggplot(kern, aes(t, survival, colour = kernel, linetype = kernel)) +
  geom_line(linewidth = 1) +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.01))) +
  scale_colour_manual(values = kernel_cols) +
  scale_linetype_manual(values = kernel_lty) +
  labs(title = "(b) Fraction of triggering still to occur",
       x = "Time since triggering event (days)",
       y = expression(S[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Combine 
#-------------------------------------------------------------------------------

fig <- (p_density | p_survival) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(fig)

#-------------------------------------------------------------------------------
# Percentage of triggering beyond $0.1$ days 
#-------------------------------------------------------------------------------

t_eval <- 1e-1
cat(sprintf("\nPercentage of triggering outstanding beyond %g days:\n", t_eval),
    sprintf("%-32s %10.3g%%\n",
            kernel_levels,
            100 * c(S_ou(t_eval, par_ou$c,  par_ou$p),
                    S_ms(t_eval, par_ms$d,  par_ms$lambda, par_ms$gamma),
                    S_rs(t_eval, par_rs$B,  par_rs$ta))), sep = "")

#-------------------------------------------------------------------------------
# Percentage of triggering beyond $1$ days 
#-------------------------------------------------------------------------------

t_eval <- 1
cat(sprintf("\nPercentage of triggering outstanding beyond %g day:\n", t_eval),
    sprintf("%-32s %10.3g%%\n",
            kernel_levels,
            100 * c(S_ou(t_eval, par_ou$c,  par_ou$p),
                    S_ms(t_eval, par_ms$d,  par_ms$lambda, par_ms$gamma),
                    S_rs(t_eval, par_rs$B,  par_rs$ta))), sep = "")

#-------------------------------------------------------------------------------
# Percentage of triggering beyond $1,000$ days 
#-------------------------------------------------------------------------------

t_eval <- 1e3
cat(sprintf("\nPercentage of triggering outstanding beyond %g days:\n", t_eval),
    sprintf("%-32s %10.3g%%\n",
            kernel_levels,
            100 * c(S_ou(t_eval, par_ou$c,  par_ou$p),
                    S_ms(t_eval, par_ms$d,  par_ms$lambda, par_ms$gamma),
                    S_rs(t_eval, par_rs$B,  par_rs$ta))), sep = "")

#-------------------------------------------------------------------------------
# Save plots
#-------------------------------------------------------------------------------

fig_src <- here("outputs", "comparison-figure")
dir.create(fig_src, recursive = TRUE, showWarnings = FALSE)
ggsave(here(fig_src, "kernel-comparison.pdf"), fig, width = 7, height = 3.5)