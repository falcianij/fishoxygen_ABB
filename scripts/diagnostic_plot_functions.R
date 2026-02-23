# scripts/diagnostic_plot_functions.R
# -------------------------------------------------------------------
# Plotting helpers for the Rubalcaba-style internal O2 model diagnostics.
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

plot_available_energy <- function(df, x = c("B_used", "f_ref", "T"), y = "pO2") {
  x <- match.arg(x)
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = E_net_norm)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(low = "#7f0000", mid = "white", high = "#08519c", midpoint = 0,
                         name = "E_avail / Cmax", na.value = "grey20") +
    labs(x = x, y = y, title = "Available energy diagnostic")
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
    scale_fill_viridis_c(option = "C", na.value = "grey20") +
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
    scale_fill_viridis_c(option = "D", na.value = "grey20") +
    labs(x = x, y = y, title = "Process-level reference diagnostics")
}
