# scripts/diagnostic_plot_functions.R
# -------------------------------------------------------------------
# Plotting + surface utilities for the ingestion-decision O2-closure model.
#
# Goals:
# - Keep Rmd thin: interpolation + masking live here.
# - Per-plot contour control (none / bins / breaks / multiple overlays).
# - Consistent with updated model scripts:
#     * f = Enc/(Enc + Cmax)  (encounter saturation)
#     * g = Hill(pO2_int)
#     * proc_real_frac = C_real/C_pot  (equals g in current model)
#     * D_SDA = a_SDA*C_real, M_exc = a_exc*C_real, A_assim=(1-a_SDA-a_exc)*C_real
#     * B_from_f inverts the full model (optimized U + O2 closure)
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# --------------------------- Guides / themes ---------------------------

legend_bar_guide <- function(barwidth_cm = 4.0, barheight_cm = 0.35) {
  guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    label.position = "bottom",
    frame.colour = "grey50",
    ticks = TRUE,
    ticks.colour = "grey50",
    frame.linewidth = 0.5,
    barwidth = grid::unit(barwidth_cm, "cm"),
    barheight = grid::unit(barheight_cm, "cm")
  )
}

theme_surface_default <- function(base_size = 9, square = TRUE) {
  th <- theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.4),
      axis.ticks = element_line(colour = "grey50"),
      axis.text = element_text(size = 7),
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.spacing.x = unit(0.8, "lines"),
      panel.spacing.y = unit(0.8, "lines")
    )
  if (isTRUE(square)) th <- th + theme(aspect.ratio = 1)
  th
}

# --------------------------- Mask helpers ---------------------------

# Standardize exclusion flags into two masks:
# - excl_o2: oxygen exclusion
# - excl_e : energetic exclusion (only if not oxygen exclusion)
add_masks <- function(df,
                      o2_col = "oxygen_exclusion",
                      e_col = "energetic_exclusion") {
  df %>%
    mutate(
      oxygen_exclusion = as.logical(round(as.numeric(.data[[o2_col]]))),
      energetic_exclusion = as.logical(round(as.numeric(.data[[e_col]]))),
      energetic_exclusion = if_else(is.na(energetic_exclusion), FALSE, energetic_exclusion),
      excl_o2 = oxygen_exclusion,
      excl_e = !excl_o2 & energetic_exclusion
    )
}

# Optional smoothing of mask boundaries by interpolating the mask as numeric.
# Requires interp_to_grid() from scripts/analysis_functions_ingestion.R to be sourced.
build_mask_surface <- function(mask_df, x, y, nx = 300, ny = 300) {
  md <- mask_df %>%
    mutate(excl_e_num = as.numeric(excl_e),
           excl_o2_num = as.numeric(excl_o2))
  
  e_s <- interp_to_grid(md, xcol = x, ycol = y, zcol = "excl_e_num", nx = nx, ny = ny)
  o_s <- interp_to_grid(md, xcol = x, ycol = y, zcol = "excl_o2_num", nx = nx, ny = ny)
  
  left_join(e_s, o_s, by = c(x, y)) %>%
    mutate(
      excl_e = is.finite(excl_e_num) & excl_e_num >= 0.5,
      excl_o2 = is.finite(excl_o2_num) & excl_o2_num >= 0.5
    )
}

# --------------------------- Interpolation helpers ---------------------------

# Interpolate a single z surface for plotting, optionally by facets.
# - Drops excluded cells and non-finite values before interpolation.
# - Returns a grid suitable for geom_raster.
interp_surface <- function(df, x, y, z, facet_cols = character(0), nx = 300, ny = 300) {
  if (length(facet_cols) == 0) {
    dd_valid <- df %>%
      filter(!excl_o2, !excl_e,
             is.finite(.data[[z]]),
             is.finite(.data[[x]]),
             is.finite(.data[[y]]))
    if (nrow(dd_valid) < 5) return(df)
    ip <- interp_to_grid(dd_valid, xcol = x, ycol = y, zcol = z, nx = nx, ny = ny)
    ip$excl_o2 <- FALSE
    ip$excl_e <- FALSE
    return(ip)
  }
  
  split_df <- split(df, interaction(df[, facet_cols, drop = FALSE], drop = TRUE))
  purrr::imap_dfr(split_df, function(dd, nm) {
    dd_valid <- dd %>%
      filter(!excl_o2, !excl_e,
             is.finite(.data[[z]]),
             is.finite(.data[[x]]),
             is.finite(.data[[y]]))
    if (nrow(dd_valid) < 5) return(dd)
    ip <- interp_to_grid(dd_valid, xcol = x, ycol = y, zcol = z, nx = nx, ny = ny)
    for (fc in facet_cols) ip[[fc]] <- dd[[fc]][1]
    ip$excl_o2 <- FALSE
    ip$excl_e <- FALSE
    ip
  })
}

# --------------------------- Contour specs ---------------------------

# Helper to build a contour spec list.
# z can be:
# - NULL (no contours)
# - character (name of column to contour)
# - a list of multiple contour layers: list(list(z="f", breaks=...), list(z="g", bins=...))
contour_spec <- function(z,
                         breaks = NULL,
                         bins = NULL,
                         colour = "black",
                         linewidth = 0.25,
                         linetype = "dashed",
                         alpha = 0.8) {
  list(z = z, breaks = breaks, bins = bins,
       colour = colour, linewidth = linewidth,
       linetype = linetype, alpha = alpha)
}

# Defaults: you can call default_contours(zname) and override as needed.
default_contours <- function(zname) {
  if (identical(zname, "f") || identical(zname, "g")) {
    contour_spec(z = zname, breaks = seq(0, 1, by = 0.2), linewidth = 0.25, alpha = 0.9)
  } else if (identical(zname, "E_net_norm")) {
    contour_spec(z = zname, breaks = seq(0, 1, by = 0.1), linewidth = 0.25, alpha = 0.9)
  } else {
    contour_spec(z = zname, bins = 5, linewidth = 0.25, alpha = 0.7)
  }
}

# --------------------------- Palette / scales ---------------------------

scale_fill_palette <- function(palette = c("blue", "viridis", "green", "purple", "red"),
                               fill_title = "",
                               limits = NULL,
                               oob = scales::squish,
                               na.value = "black",
                               guide = legend_bar_guide()) {
  palette <- match.arg(palette)
  if (palette == "blue") {
    return(scale_fill_gradient(low = "white", high = "#67A5CA",
                               name = fill_title, limits = limits,
                               oob = oob, na.value = na.value, guide = guide))
  }
  if (palette == "viridis") {
    return(scale_fill_viridis_c(option = "D",
                                name = fill_title, limits = limits,
                                oob = oob, na.value = na.value, guide = guide))
  }
  if (palette == "green") {
    return(scale_fill_gradient(high = "white", low = "#238b45",
                               name = fill_title, limits = limits,
                               oob = oob, na.value = na.value, guide = guide))
  }
  if (palette == "purple") {
    return(scale_fill_gradient(high = "white", low = "#6a51a3",
                               name = fill_title, limits = limits,
                               oob = oob, na.value = na.value, guide = guide))
  }
  scale_fill_gradient(low = "white", high = "#cb181d",
                      name = fill_title, limits = limits,
                      oob = oob, na.value = na.value, guide = guide)
}

# --------------------------- Main surface plotter ---------------------------

# plot_surface_masked() is your workhorse.
# It:
# - expects df has been add_masks()'d (has excl_o2/excl_e)
# - renders excluded regions (black/grey) and values elsewhere
# - optionally overlays contours (one or multiple) with per-plot control
#
# Key design choice: we *do not* interpolate by default; you choose whether to call interp_surface()
# and pass that in. That’s what makes it “pliable” plot-by-plot.
plot_surface_masked <- function(df, x, y, z,
                                palette = c("blue", "viridis", "green", "purple", "red"),
                                fill_title = "",
                                limits = NULL,
                                x_limits = NULL,
                                y_limits = NULL,
                                x_trans = c("identity", "log10"),
                                mask_df = NULL,
                                smooth_mask = FALSE,
                                mask_nx = 300,
                                mask_ny = 300,
                                # contours:
                                contours = NULL,       # NULL=no contours; "auto"=default_contours(z); or list(specs)
                                contour_data = NULL,   # if NULL use df filtered to valid; else provide a df
                                # theme:
                                theme = theme_surface_default()) {
  
  palette <- match.arg(palette)
  x_trans <- match.arg(x_trans)
  if (is.null(mask_df)) mask_df <- df
  
  # Mask plotting dataset (optionally smoothed)
  mask_plot <- if (isTRUE(smooth_mask)) build_mask_surface(mask_df, x = x, y = y, nx = mask_nx, ny = mask_ny) else mask_df
  
  # Value subsets (only non-excluded)
  val_df <- df %>%
    filter(!excl_o2, !excl_e) %>%
    mutate(.zval = .data[[z]])
  
  val_na  <- val_df %>% filter(!is.finite(.zval) | is.na(.zval))
  val_neg <- val_df %>% filter(is.finite(.zval) & .zval < 0)
  val_pos <- val_df %>% filter(is.finite(.zval) & .zval >= 0)
  
  p <- ggplot() +
    geom_raster(data = val_na,  aes(.data[[x]], .data[[y]]), fill = "black", interpolate = TRUE) +
    geom_raster(data = val_neg, aes(.data[[x]], .data[[y]]), fill = "grey70", interpolate = TRUE) +
    geom_raster(data = val_pos, aes(.data[[x]], .data[[y]], fill = .zval), interpolate = TRUE)
  
  # Contours
  if (!is.null(contours)) {
    if (identical(contours, "auto")) contours <- list(default_contours(z))
    
    # Allow single spec or list of specs
    if (!is.list(contours) || (is.list(contours) && !is.null(contours$z))) contours <- list(contours)
    
    cdat <- contour_data %||% val_pos
    
    for (sp in contours) {
      if (is.null(sp$z) || !nzchar(sp$z)) next
      if (!(sp$z %in% names(cdat))) next
      
      if (!is.null(sp$breaks)) {
        p <- p + geom_contour(
          data = cdat,
          aes(.data[[x]], .data[[y]], z = .data[[sp$z]]),
          breaks = sp$breaks,
          linetype = sp$linetype, linewidth = sp$linewidth,
          colour = sp$colour, alpha = sp$alpha
        )
      } else {
        p <- p + geom_contour(
          data = cdat,
          aes(.data[[x]], .data[[y]], z = .data[[sp$z]]),
          bins = sp$bins %||% 5,
          linetype = sp$linetype, linewidth = sp$linewidth,
          colour = sp$colour, alpha = sp$alpha
        )
      }
    }
  }
  
  # Exclusions on top (so they always mask)
  p <- p +
    geom_raster(data = mask_plot %>% filter(excl_e),  aes(.data[[x]], .data[[y]]), fill = "grey70", interpolate = TRUE) +
    geom_raster(data = mask_plot %>% filter(excl_o2), aes(.data[[x]], .data[[y]]), fill = "black",  interpolate = TRUE) +
    scale_fill_palette(palette = palette, fill_title = fill_title, limits = limits)
  
  # Axes
  if (x_trans == "log10") {
    p <- p + scale_x_log10(expand = c(0, 0), limits = x_limits, labels = label_scientific())
  } else {
    p <- p + scale_x_continuous(expand = c(0, 0), limits = x_limits)
  }
  
  p + scale_y_continuous(expand = c(0, 0), limits = y_limits) + theme
}

# --------------------------- Convenience diagnostics ---------------------------

# "available energy" uses E_net_norm with NA/neg/pos handling; now just a thin wrapper.
plot_available_energy <- function(df, x = c("B_used", "f_ref", "T"), y = "pO2",
                                  contours = "auto",
                                  ...) {
  x <- match.arg(x)
  plot_surface_masked(df, x = x, y = y, z = "E_net_norm",
                      palette = "blue", fill_title = "E/Cmax",
                      contours = contours,
                      ...) +
    labs(x = x, y = y, title = "Available energy (masked)")
}

# Metabolic allocation facets (uses updated names: M_m, M_act, D_SDA, A_assim, M_exc)
plot_metabolic_allocation <- function(df, x = c("B_used", "f_ref", "T"), y = "pO2",
                                      interpolate = TRUE, nx = 300, ny = 300,
                                      facet_scales = "free",
                                      theme = theme_surface_default()) {
  x <- match.arg(x)
  
  dd <- df %>%
    transmute(xvar = .data[[x]], yvar = .data[[y]],
              Maintenance = M_m,
              Activity = M_act,
              SDA = D_SDA,
              Excretion = M_exc,
              Assimilation = A_assim) %>%
    pivot_longer(cols = c(Maintenance, Activity, SDA, Excretion, Assimilation),
                 names_to = "process", values_to = "value")
  
  if (isTRUE(interpolate)) {
    # Interpolate per facet for smoother rasters
    dd <- dd %>% group_by(process) %>% group_modify(\(d, key) {
      d <- add_masks(dplyr::mutate(d, excl_o2 = FALSE, excl_e = FALSE,
                                   oxygen_exclusion = FALSE, energetic_exclusion = FALSE),
                     o2_col = "oxygen_exclusion", e_col = "energetic_exclusion")
      interp_surface(d, x = "xvar", y = "yvar", z = "value", facet_cols = character(0), nx = nx, ny = ny)
    }) %>% ungroup()
  }
  
  ggplot(dd, aes(x = xvar, y = yvar, fill = value)) +
    geom_raster(interpolate = TRUE) +
    facet_wrap(~process, scales = facet_scales) +
    scale_fill_viridis_c(option = "C", na.value = "black", guide = legend_bar_guide()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x, y = y, title = "Metabolic allocation") +
    theme
}

plot_process_reference <- function(df, x = c("B_used", "f_ref"), y = "pO2",
                                   facet_scales = "free",
                                   theme = theme_surface_default()) {
  x <- match.arg(x)
  
  dd <- df %>%
    transmute(xvar = .data[[x]], yvar = .data[[y]],
              Encounter = Enc,
              Cmax = Cmax,
              C_pot = C_pot,
              C_real = C_real,
              g = g,
              f = f,
              `C_real/C_pot` = proc_real_frac,
              `O2 supply` = O2_supply,
              `O2 demand` = O2_demand,
              `pO2_int` = pO2_int) %>%
    pivot_longer(cols = -c(xvar, yvar), names_to = "diagnostic", values_to = "value")
  
  ggplot(dd, aes(x = xvar, y = yvar, fill = value)) +
    geom_raster(interpolate = TRUE) +
    facet_wrap(~diagnostic, scales = facet_scales) +
    scale_fill_viridis_c(option = "D", na.value = "black", guide = legend_bar_guide()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = x, y = y, title = "Process-level diagnostics") +
    theme
}

# Bivariate f–g limitation mix (kept, but now accepts contour control)
plot_fg_bivariate <- function(df, x, y,
                              x_limits = NULL, y_limits = NULL,
                              x_trans = c("identity", "log10"),
                              mask_df = NULL,
                              # contour control:
                              contour_f = TRUE,
                              contour_g = TRUE,
                              contour_fg_eq = TRUE,
                              contour_fg_mult = TRUE,
                              breaks_fg = seq(0, 1, by = 0.2),
                              theme = theme_surface_default()) {
  
  x_trans <- match.arg(x_trans)
  if (is.null(mask_df)) mask_df <- df
  
  if (!all(c("f", "g", "excl_o2", "excl_e") %in% names(df))) {
    stop("plot_fg_bivariate() requires columns: f, g, excl_o2, excl_e (run add_masks() first).")
  }
  
  dd <- df %>%
    mutate(
      f = if_else(is.finite(f), pmin(pmax(f, 0), 1), NA_real_),
      g = if_else(is.finite(g), pmin(pmax(g, 0), 1), NA_real_),
      fg_eq = f - g,
      fg_mult = f * g,
      prey_lim = 1 - f,
      oxy_lim = 1 - g,
      r_ch = pmin(pmax(1 - (0.86 * prey_lim + 0.58 * oxy_lim), 0), 1),
      g_ch = pmin(pmax(1 - (0.26 * prey_lim + 0.81 * oxy_lim), 0), 1),
      b_ch = pmin(pmax(1 - (0.86 * prey_lim + 0.30 * oxy_lim), 0), 1),
      fg_col = rgb(red = replace_na(r_ch, 1), green = replace_na(g_ch, 1), blue = replace_na(b_ch, 1))
    )
  
  dd_val <- dd %>% filter(!excl_o2, !excl_e)
  dd_na  <- dd_val %>% filter(!is.finite(f) | !is.finite(g))
  dd_ok  <- dd_val %>% filter(is.finite(f) & is.finite(g))
  
  p <- ggplot() +
    geom_raster(data = dd_na, aes(.data[[x]], .data[[y]]), fill = "grey70", interpolate = TRUE) +
    geom_raster(data = dd_ok, aes(.data[[x]], .data[[y]], fill = fg_col), interpolate = TRUE) +
    scale_fill_identity(name = "Limitation mix", guide = "none")
  
  if (isTRUE(contour_f)) {
    p <- p + geom_contour(data = dd_ok, aes(.data[[x]], .data[[y]], z = f),
                          breaks = breaks_fg, colour = "#238b45",
                          linewidth = 0.35, linetype = "dashed", alpha = 1)
  }
  if (isTRUE(contour_g)) {
    p <- p + geom_contour(data = dd_ok, aes(.data[[x]], .data[[y]], z = g),
                          breaks = breaks_fg, colour = "#6a51a3",
                          linewidth = 0.35, linetype = "dashed", alpha = 1)
  }
  if (isTRUE(contour_fg_eq)) {
    p <- p + geom_contour(data = dd_ok, aes(.data[[x]], .data[[y]], z = fg_eq),
                          breaks = 0, colour = "black", linewidth = 0.75)
  }
  if (isTRUE(contour_fg_mult)) {
    p <- p + geom_contour(data = dd_ok, aes(.data[[x]], .data[[y]], z = fg_mult),
                          breaks = breaks_fg, colour = "black", linewidth = 0.35, linetype = "dashed")
  }
  
  # Masks on top
  p <- p +
    geom_raster(data = mask_df %>% filter(excl_e),  aes(.data[[x]], .data[[y]]), fill = "grey70", interpolate = TRUE) +
    geom_raster(data = mask_df %>% filter(excl_o2), aes(.data[[x]], .data[[y]]), fill = "black",  interpolate = TRUE)
  
  if (x_trans == "log10") {
    p <- p + scale_x_log10(expand = c(0, 0), limits = x_limits, labels = label_scientific())
  } else {
    p <- p + scale_x_continuous(expand = c(0, 0), limits = x_limits)
  }
  
  p + scale_y_continuous(expand = c(0, 0), limits = y_limits) + theme
}