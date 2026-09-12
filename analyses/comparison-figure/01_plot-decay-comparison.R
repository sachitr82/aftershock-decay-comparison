#===============================================================================
# Comparison of the three candidate decay functions
#
# Output is outputs/comparison-figure/decay-comparison.pdf, containing 2 plots:
#   (a) normalised densities  f_k(t) = g_k(t) / G_k(0, Inf)      [log-log]
#   (b) survival functions    S_k(t) = 1 - G_k(t) / G_k(0, Inf)  [log-x]
#
# Normalised forms restrict comparisons to shape of decay.  
# Survival functions display the fraction of an event's triggering still
# outstanding at time t.
#
# Parameter choices based on Hainzl & Christophersen (2017, Figure 1).
#===============================================================================

#-------------------------------------------------------------------------------
# Package dependencies
#-------------------------------------------------------------------------------

library(ggplot2)
library(here)
library(patchwork)
library(scales)

#-------------------------------------------------------------------------------
# Normalised densities and survival functions: log1p/expm1 for numerical stability
#-------------------------------------------------------------------------------

# Omori-Utsu 
f_ou <- function(t, c, p) {
  (p - 1) / c * (1 + t / c)^(-p)
}

S_ou <- function(t, c, p) {
  (1 + t / c)^(1 - p)
}

# Modified stretched exponential 
#  delta(t) = (d + t)^gamma - d^gamma, written to avoid cancellation at small t
delta_ms <- function(t, d, gamma) {
  d^gamma * expm1(gamma * log1p(t / d))
}

f_ms <- function(t, d, rho, gamma) {
  rho * gamma * (d + t)^(gamma - 1) * exp(-rho * delta_ms(t, d, gamma))
}

S_ms <- function(t, d, rho, gamma) {
  exp(-rho * delta_ms(t, d, gamma))
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
# Parameters (based on Hainzl & Christophersen 2017, Figure 1)
#-------------------------------------------------------------------------------

par_ou <- list(c = 0.01,  p  = 1.1)
par_ms <- list(d = 0.01, rho = 1, gamma = 0.2)
par_rs <- list(B = 0.9999, ta = 100)

decay_levels <- c("Omori\u2013Utsu",
                   "Modified stretched exponential",
                   "Rate-state")

t_grid <- 10^seq(-4, 4, length.out = 3000)

#-------------------------------------------------------------------------------
# Evaluate decay functions and order factor levels for legend
#-------------------------------------------------------------------------------

dec <- rbind(
  data.frame(t = t_grid, decay = decay_levels[1],
             density  = f_ou(t_grid, par_ou$c, par_ou$p),
             survival = S_ou(t_grid, par_ou$c, par_ou$p)),
  data.frame(t = t_grid, decay = decay_levels[2],
             density  = f_ms(t_grid, par_ms$d, par_ms$rho, par_ms$gamma),
             survival = S_ms(t_grid, par_ms$d, par_ms$rho, par_ms$gamma)),
  data.frame(t = t_grid, decay = decay_levels[3],
             density  = f_rs(t_grid, par_rs$B, par_rs$ta),
             survival = S_rs(t_grid, par_rs$B, par_rs$ta))
)

dec$decay <- factor(dec$decay, levels = decay_levels)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

# Decay identifiers - colour and line type
decay_cols <- setNames(c("steelblue", "orange", "springgreen4"), decay_levels)
decay_lty  <- setNames(c("solid", "dashed", "dotted"), decay_levels)
decay_scales <- list(scale_colour_manual(values = decay_cols),
                     scale_linetype_manual(values = decay_lty))

# Format tick labels as, e.g., 10^2 instead 1e2
log_lab <- trans_format("log10", math_format(10^.x))

# Common theme choices for both plots
base_theme <- theme_bw(base_size = 10) +
              theme(legend.title     = element_blank(),
                    legend.key.width = unit(0.9, "lines"),
                    legend.key.spacing.x = unit(0.5, "cm"),
                    legend.text      = element_text(size = 11),
                    panel.grid.minor = element_line(linewidth = 0.15))

# Logarithmic x-axis used in both plots
x_scale <- scale_x_log10(breaks = 10^seq(-4, 4, by = 1), labels = log_lab,
                         expand = expansion(mult = 0))

#-------------------------------------------------------------------------------
# Panel (a): normalised densities (log-log)
#-------------------------------------------------------------------------------

p_density <- ggplot(dec, aes(t, density, colour = decay, linetype = decay)) +
  geom_line(linewidth = 1) +
  x_scale + 
  scale_y_log10(breaks = 10^seq(-8, 2, by = 2), labels = log_lab) +
  coord_cartesian(ylim = c(1e-8, NA)) +
  decay_scales + 
  labs(x = "Time since triggering event (days)", y = expression(f[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Panel (b): survival functions (log x - linear y)
#-------------------------------------------------------------------------------

p_survival <- ggplot(dec, aes(t, survival, colour = decay, linetype = decay)) +
  geom_line(linewidth = 1) +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
  decay_scales +
  labs(x = "Time since triggering event (days)", y = expression(S[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Combine 
#-------------------------------------------------------------------------------

fig <- (p_density | p_survival) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", 
        plot.tag = element_text(size = 10, face = "plain"))

#-------------------------------------------------------------------------------
# Save plots
#-------------------------------------------------------------------------------

fig_src <- here("outputs", "comparison-figure")
dir.create(fig_src, recursive = TRUE, showWarnings = FALSE)
ggsave(here(fig_src, "decay-comparison.pdf"), fig, width = 7, height = 3.5)
message("Saved decay-comparison.pdf to ", fig_src)