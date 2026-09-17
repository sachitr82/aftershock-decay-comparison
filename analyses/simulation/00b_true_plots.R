#===============================================================================
# Comparison of generating temporal decay functions
#
#   (a) unnormalised temporal decay functions g_k(t)              [log-log]
#   (b) cumulative triggering masses G_k(0,t)                     [log-x]
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and fixed simulation design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(ggplot2)
library(here)
library(patchwork)
library(scales)
source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
## Temporal parameter lists
#-------------------------------------------------------------------------------

theta_ou <- list(c = c_true, p = p_true)

theta_mse <- list(d = d_true, rho = rho_true,  gamma = gamma_true)

theta_rs <- list(B = B_true, ta = ta_true)

#-------------------------------------------------------------------------------
# Plotting grid
#-------------------------------------------------------------------------------

t_plot_max <- max(1e4, T_fit_end)

t_grid <- sort(unique(c(10^seq(-4, log10(t_plot_max), length.out = 3000), T_fit_end)))

#-------------------------------------------------------------------------------
## Evaluate g_k(t)
#-------------------------------------------------------------------------------

g_ou <- ETAS.inlabru::temporal_decay(dt = t_grid, theta = theta_ou, 
                                      form = "ou")

g_mse <- ETAS.inlabru::temporal_decay(dt = t_grid, theta = theta_mse, 
                                       form = "mse")

g_rs <- ETAS.inlabru::temporal_decay(dt = t_grid, theta = theta_rs, 
                                      form = "rate_state")

#-------------------------------------------------------------------------------
## Evaluate G_k(0,t)
#-------------------------------------------------------------------------------

cum_mass <- function(t, theta, form) {
  vapply(t, function(tt) {
      ETAS.inlabru::temporal_decay_integral(a = 0, b = tt, theta = theta,
                                             form = form)
    },numeric(1))
}

G_ou <- cum_mass(t_grid, theta_ou, "ou")
G_mse <- cum_mass(t_grid, theta_mse, "mse")
G_rs <- cum_mass(t_grid, theta_rs, "rate_state")

#-------------------------------------------------------------------------------
# Assemble plotting data
#-------------------------------------------------------------------------------

decay_levels <- c("Omori\u2013Utsu",
                  "Modified stretched exponential",
                  "Rate-state")

dec <- rbind(
  data.frame(t = t_grid, decay = decay_levels[1], g = g_ou, G = G_ou),
  data.frame( t = t_grid, decay = decay_levels[2], g = g_mse, G = G_mse),
  data.frame(t = t_grid, decay = decay_levels[3], g = g_rs, G = G_rs)
)

dec$decay <- factor(dec$decay, levels = decay_levels)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

decay_cols <- setNames(c("steelblue", "orange", "springgreen4"), decay_levels)

decay_lty <- setNames(c("solid", "dashed", "dotted"), decay_levels)

log_lab <- scales::trans_format("log10", scales::math_format(10^.x))

base_theme <- theme_bw(base_size = 10) +
  theme(legend.title     = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(0.9, "lines"),
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.spacing.x = unit(0.6, "cm"),
        legend.text      = element_text(size = 11),
        panel.grid.minor = element_line(linewidth = 0.15),
        plot.title       = element_text(size = 10, face = "plain"))

x_scale <- scale_x_log10(breaks = 10^seq(-4, 4, by = 1), labels = log_lab,
                         expand = expansion(mult = 0))

#-------------------------------------------------------------------------------
## Panel (a): unnormalised temporal decay functions
#-------------------------------------------------------------------------------

p_trigger <- ggplot(dec, aes(t, g, colour = decay, linetype = decay)) +
  geom_line(linewidth = 0.9) +
  x_scale +
  scale_y_log10(breaks = 10^seq(-12, 0, by = 2), labels = log_lab) +
  coord_cartesian(ylim = c(1e-12, 1)) + 
  scale_colour_manual(values = decay_cols) +
  scale_linetype_manual(values = decay_lty) +
  labs(x = "Time since triggering event (days)", y = expression(g[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Panel (b): cumulative integrated triggering mass
#-------------------------------------------------------------------------------

p_mass <- ggplot(dec, aes(t, G, colour = decay, linetype = decay)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = T_fit_end, linetype = "longdash", linewidth = 0.45) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_colour_manual(values = decay_cols) +
  scale_linetype_manual(values = decay_lty) +
  labs(x = "Time since triggering event (days)", y = expression(G[k](0,t))) +
  base_theme

#-------------------------------------------------------------------------------
# Combine
#-------------------------------------------------------------------------------

fig <- (p_trigger | p_mass) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(legend.position = "bottom")) &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 10, face = "plain"))

#-------------------------------------------------------------------------------
# Save figure
#-------------------------------------------------------------------------------

figfile <- file.path(design_figure_dir, "synthetic-decay-design.pdf")

ggsave(filename = figfile, plot = fig,
       width = 7, height = 3.5)

message("Saved ", figfile)