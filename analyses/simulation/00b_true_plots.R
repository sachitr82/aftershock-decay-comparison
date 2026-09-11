#===============================================================================
# Comparison of generating temporal triggering functions
#
#   (a) unnormalised temporal triggering functions g_k(t)         [log-log]
#   (b) cumulative triggering masses G_k(0,t)                     [log-x]
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and fixed simulation design
#-------------------------------------------------------------------------------

library(ggplot2)
library(here)
library(patchwork)
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

g_ou <- ETAS.inlabru::temporal_kernel(dt = t_grid, theta = theta_ou, 
                                      kernel = "ou")

g_mse <- ETAS.inlabru::temporal_kernel(dt = t_grid, theta = theta_mse, 
                                       kernel = "mse")

g_rs <- ETAS.inlabru::temporal_kernel(dt = t_grid, theta = theta_rs, 
                                      kernel = "rate_state")

#-------------------------------------------------------------------------------
## Evaluate G_k(0,t)
#-------------------------------------------------------------------------------

cum_mass <- function(t, theta, kernel) {
  vapply(t, function(tt) {
      ETAS.inlabru::temporal_kernel_integral(a = 0, b = tt, theta = theta,
                                             kernel = kernel)
    },numeric(1))
}

G_ou <- cum_mass(t_grid, theta_ou, "ou")
G_mse <- cum_mass(t_grid, theta_mse, "mse")
G_rs <- cum_mass(t_grid, theta_rs, "rate_state")

#-------------------------------------------------------------------------------
# Assemble plotting data
#-------------------------------------------------------------------------------

kernel_levels <- c(
  "Omori\u2013Utsu",
  "Modified stretched exponential",
  "Rate-state")

kern <- rbind(
  data.frame(t = t_grid, kernel = kernel_levels[1], g = g_ou, G = G_ou),
  data.frame( t = t_grid, kernel = kernel_levels[2], g = g_mse, G = G_mse),
  data.frame(t = t_grid, kernel = kernel_levels[3], g = g_rs, G = G_rs)
)

kern$kernel <- factor(kern$kernel, levels = kernel_levels)

#-------------------------------------------------------------------------------
# Plotting parameters
#-------------------------------------------------------------------------------

kernel_cols <- setNames(c("steelblue", "orange", "springgreen4"), kernel_levels)

kernel_lty <- setNames(c("solid", "dashed", "dotted"), kernel_levels)

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
## Panel (a): unnormalised temporal triggering functions
#-------------------------------------------------------------------------------

p_trigger <- ggplot(kern, aes(t, g, colour = kernel, linetype = kernel)) +
  geom_line(linewidth = 0.9) +
  x_scale +
  scale_y_log10(breaks = 10^seq(-12, 0, by = 2), labels = log_lab) +
  coord_cartesian(ylim = c(1e-12, 1)) + 
  scale_colour_manual(values = kernel_cols) +
  scale_linetype_manual(values = kernel_lty) +
  labs(x = "Time since triggering event (days)", y = expression(g[k](t))) +
  base_theme

#-------------------------------------------------------------------------------
# Panel (b): cumulative integrated triggering mass
#-------------------------------------------------------------------------------

p_mass <- ggplot(kern, aes(t, G, colour = kernel, linetype = kernel)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = T_fit_end, linetype = "longdash", linewidth = 0.45) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_colour_manual(values = kernel_cols) +
  scale_linetype_manual(values = kernel_lty) +
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

print(fig)

#-------------------------------------------------------------------------------
# Report matched finite-window mass and lifetime masses
#-------------------------------------------------------------------------------

G_inf_plot <- c(
  ou = ETAS.inlabru::temporal_kernel_total_mass(theta = theta_ou,kernel = "ou"),
  mse = ETAS.inlabru::temporal_kernel_total_mass(theta = theta_mse,kernel = "mse"),
  rate_state = ETAS.inlabru::temporal_kernel_total_mass(theta = theta_rs, 
                                                        kernel = "rate_state")
)

print(G_inf_plot)

#-------------------------------------------------------------------------------
# Save figure
#-------------------------------------------------------------------------------

fig_dir <- here("outputs", "simulation", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(filename = here(fig_dir, "synthetic-triggering-design.pdf"), plot = fig,
       width = 7, height = 3.5)