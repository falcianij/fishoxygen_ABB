# scripts/analysis_functions_ingestion.R
# -------------------------------------------------------------------
# Analysis helpers for Rubalcaba-style internal O2 closure model.
# -------------------------------------------------------------------

library(interp)
library(dplyr)

compute_state_surface <- function(tr,
                                  masses,
                                  Tseq,
                                  pO2seq,
                                  prey_units = c("biomass", "prey_index", "feeding_level", "powerlaw"),
                                  preys = NULL,
                                  intercepts = 1.6,
                                  slopes = -0.08,
                                  m0 = 1,
                                  pair_mode = c("auto", "cross", "zip"),
                                  normalize_E = TRUE,
                                  progress = TRUE,
                                  T_ref = tr$T_ref,
                                  pO2_ref = 25,
                                  U_lo = 1e-5, U_mech = 5,
                                  enforce_O2_feasible = TRUE,
                                  u_prey = 0, D = 3) {

  prey_units <- match.arg(prey_units)
  pair_mode <- match.arg(pair_mode)

  make_spec_cross <- function() {
    if (prey_units %in% c("biomass", "prey_index", "feeding_level")) {
      tidyr::expand_grid(mass = masses, prey_arg = preys)
    } else {
      tidyr::expand_grid(mass = masses, intercept = intercepts, slope = slopes)
    }
  }
  make_spec_zip <- function() {
    if (prey_units %in% c("biomass", "prey_index", "feeding_level")) {
      tibble::tibble(mass = rep_len(masses, max(length(masses), length(preys))),
                     prey_arg = rep_len(preys, max(length(masses), length(preys))))
    } else {
      tibble::tibble(mass = rep_len(masses, max(length(masses), length(intercepts), length(slopes))),
                     intercept = rep_len(intercepts, max(length(masses), length(intercepts), length(slopes))),
                     slope = rep_len(slopes, max(length(masses), length(intercepts), length(slopes))))
    }
  }

  if (pair_mode == "auto") pair_mode <- if (length(masses) > 1 && !is.null(preys) && length(preys) == length(masses)) "zip" else "cross"
  spec <- if (pair_mode == "zip") make_spec_zip() else make_spec_cross()

  resolve_B <- function(m, row, Tval = NULL) {
    if (prey_units == "biomass") {
      list(B = row$prey_arg, mode = "biomass", f_ref = NA_real_)
    } else if (prey_units == "prey_index" || prey_units == "feeding_level") {
      f_ref <- pmin(pmax(row$prey_arg, .Machine$double.eps), 1 - .Machine$double.eps)
      T_use <- ifelse(is.null(Tval), T_ref, Tval)
      B <- B_from_f(f_ref = f_ref, m = m, tr = tr,
                    T_ref = T_use, pO2_ref = pO2_ref,
                    U_lo = U_lo, U_mech = U_mech,
                    enforce_O2_feasible = enforce_O2_feasible,
                    u_prey = u_prey, D = D)
      list(B = B, mode = prey_units, f_ref = row$prey_arg)
    } else {
      list(B = B_from_powerlaw_spectrum(m, row$intercept, row$slope, m0 = m0), mode = "powerlaw", f_ref = NA_real_)
    }
  }

  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(spec), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  blocks <- vector("list", nrow(spec))
  for (i in seq_len(nrow(spec))) {
    m <- spec$mass[i]
    grid <- tidyr::expand_grid(T = Tseq, pO2 = pO2seq)

    rr <- apply(as.matrix(grid), 1, function(row) {
      Tval <- as.numeric(row[[1]])
      pO2v <- as.numeric(row[[2]])
      B_row <- resolve_B(m, spec[i, , drop = FALSE], Tval = if (prey_units == "feeding_level") Tval else NULL)

      st <- find_U_opt(pO2_env = pO2v, T = Tval, m = m, prey = B_row$B, tr = tr,
                       u_prey = u_prey, D = D, U_lo = U_lo, U_mech = U_mech)

      c(B_used = B_row$B, f_ref = B_row$f_ref,
        U_opt = st$U_opt, E_net = st$E_net, Cmax = st$Cmax,
        Enc = st$Enc, C_pot = st$C_pot, C_real = st$C_real, I = st$I,
        f = st$f, g = st$g, pO2_int = st$pO2_int,
        consump = st$consump, A_assim = st$A_assim, M_m = st$M_m, M_act = st$M_act, D_SDA = st$D_SDA, M_exc = st$M_exc,
        O2_supply = st$O2_supply, O2_demand = st$O2_demand, O2_margin = st$O2_margin,
        oxygen_exclusion = as.numeric(st$oxygen_exclusion), energetic_exclusion = as.numeric(st$energetic_exclusion))
    })

    blk <- dplyr::bind_cols(tibble::tibble(mass = m), grid, as.data.frame(t(rr)))
    blocks[[i]] <- blk
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, i)
  }

  res <- dplyr::bind_rows(blocks)
  if (isTRUE(normalize_E)) {
    res <- dplyr::mutate(res,
      E_net_norm = dplyr::if_else(is.finite(E_net) & is.finite(Cmax) & Cmax != 0, E_net / Cmax, NA_real_),
      proc_real_frac = dplyr::if_else(is.finite(C_real) & is.finite(Cmax) & Cmax != 0, C_real / Cmax, NA_real_)
    )
  }
  tibble::as_tibble(res)
}

add_regime_labels_ingestion <- function(df,
                                        g_thresh = 0.95,
                                        f_low = 0.3,
                                        keep_helpers = FALSE) {
  out <- df %>%
    dplyr::mutate(
      oxygen_exclusion = as.logical(round(as.numeric(.data$oxygen_exclusion))),
      energetic_exclusion = as.logical(round(as.numeric(.data$energetic_exclusion))),
      exclusion = oxygen_exclusion | energetic_exclusion,
      prey_limited = !exclusion & is.finite(.data$f) & is.finite(.data$g) & (.data$f < f_low) & (.data$g >= g_thresh),
      oxygen_limited = !exclusion & is.finite(.data$g) & (.data$g < g_thresh) & is.finite(.data$f) & (.data$f >= f_low),
      co_limited = !exclusion & is.finite(.data$g) & is.finite(.data$f) & !prey_limited & !oxygen_limited,
      regime = dplyr::case_when(
        oxygen_exclusion ~ "Oxygen exclusion",
        energetic_exclusion ~ "Energetic exclusion",
        prey_limited ~ "Prey limited",
        oxygen_limited ~ "Oxygen-limited processing",
        co_limited ~ "Co-limited",
        TRUE ~ "Unclassified"
      )
    )

  if (!keep_helpers) out <- dplyr::select(out, -exclusion, -prey_limited, -oxygen_limited, -co_limited)
  out
}

interp_to_grid <- function(df, xcol, ycol, zcol, nx = 500, ny = 500) {
  x <- df[[xcol]]; y <- df[[ycol]]; z <- df[[zcol]]
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  xo <- seq(min(x[ok]), max(x[ok]), length.out = nx)
  yo <- seq(min(y[ok]), max(y[ok]), length.out = ny)
  ip <- interp::interp(x[ok], y[ok], z[ok], xo = xo, yo = yo, linear = TRUE, extrap = FALSE)
  out <- expand.grid(stats::setNames(list(ip$x, ip$y), c(xcol, ycol)))
  out[[zcol]] <- as.vector(ip$z)
  out
}
