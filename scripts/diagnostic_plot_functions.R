# scripts/diagnostic_plot_functions.R
# -------------------------------------------------------------------
# Plotting helpers for the Rubalcaba-style internal O2 model diagnostics.
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

legend_bar_guide <- function() {
  guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    label.position = "bottom",
    frame.colour = "grey50",
    ticks = TRUE,
    ticks.colour = "grey50",
    frame.linewidth = 0.5,
    barwidth = grid::unit(4.0, "cm"),
    barheight = grid::unit(0.35, "cm")
  )
}

plot_available_energy <- function(df, x = c("B_used", "f_ref", "T"), y = "pO2") {
  x <- match.arg(x)

  dd <- df %>% mutate(.z = E_net_norm)
  dd_na <- dd %>% filter(!is.finite(.z) | is.na(.z))
  dd_neg <- dd %>% filter(is.finite(.z) & .z < 0)
  dd_pos <- dd %>% filter(is.finite(.z) & .z >= 0) %>% mutate(.zpos = .z)

  ggplot() +
    geom_raster(data = dd_na, aes(.data[[x]], .data[[y]]), fill = "black", interpolate = TRUE) +
    geom_raster(data = dd_neg, aes(.data[[x]], .data[[y]]), fill = "gray60", interpolate = TRUE) +
    geom_raster(data = dd_pos, aes(.data[[x]], .data[[y]], fill = pmin(.zpos, 1)), interpolate = TRUE) +
    scale_fill_gradient(
      low = "white", high = "#67A5CA", na.value = "black",
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = label_number(accuracy = 0.1),
      oob = squish,
      name = expression("Available energy"),
      guide = legend_bar_guide()
    ) +
    geom_contour(
      inherit.aes = FALSE,
      data = dd_pos, aes(x = .data[[x]], y = .data[[y]], z = .zpos),
      breaks = seq(0, 1, by = 0.1), colour = "black", linewidth = 0.3, alpha = 1, linetype = "dashed"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x, y = y, title = "Available energy diagnostic") +
    theme_minimal(base_size = 9) +
    theme(aspect.ratio = 1, legend.position = "bottom", legend.box = "horizontal")
}

plot_metabolic_allocation <- function(df, x = c("B_used", "f_ref", "T"), y = "pO2") {
  x <- match.arg(x)
  dd <- df %>%
    transmute(xvar = .data[[x]], yvar = .data[[y]],
              Maintenance = M_m,
              Activity = M_act,
              SDA = D_SDA,
              Assimilation = consump) %>%
    pivot_longer(cols = c(Maintenance, Activity, SDA, Assimilation), names_to = "process", values_to = "value")

  ggplot(dd, aes(x = xvar, y = yvar, fill = value)) +
    geom_raster(interpolate = TRUE) +
    facet_wrap(~process, scales = "free") +
    scale_fill_viridis_c(option = "C", na.value = "black", guide = legend_bar_guide()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x, y = y, title = "Metabolic allocation diagnostics")
}

plot_process_reference <- function(df, x = c("B_used", "f_ref"), y = "pO2") {
  x <- match.arg(x)
  dd <- df %>%
    transmute(xvar = .data[[x]], yvar = .data[[y]],
              Encounter = Enc,
              Cmax = Cmax,
              C_pot = C_pot,
              C_real = C_real,
              g = g,
              f = f,
              `C_real/Cmax` = proc_real_frac,
              `O2 supply` = O2_supply,
              `O2 demand` = O2_demand,
              `pO2_int` = pO2_int) %>%
    pivot_longer(cols = -c(xvar, yvar), names_to = "diagnostic", values_to = "value")

  ggplot(dd, aes(x = xvar, y = yvar, fill = value)) +
    geom_raster(interpolate = TRUE) +
    facet_wrap(~diagnostic, scales = "free") +
    scale_fill_viridis_c(option = "D", na.value = "black", guide = legend_bar_guide()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x, y = y, title = "Process-level reference diagnostics")
}
