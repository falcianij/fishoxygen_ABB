# scripts/analysis_functions.R
# -------------------------------------------------------------------
# Bioenergetics fish model functions for analysis
#
# Patch (2026-01-26):
#   - Add scenario controls: u_prey (m/s), D (2 or 3)
#   - Thread u_prey, D through B_from_f() and find_U_opt()
#   - Remove add_B_gradients()
#   - Simplify regimes: Oxygen limited when S≈D, prey limited otherwise; keep exclusions
# -------------------------------------------------------------------

library(interp)
library(dplyr)

# -------------------------------------------------------------------
# Compute state surface across (T, pO2) for masses and prey settings
# -------------------------------------------------------------------

compute_state_surface <- function(tr,
                                  masses,
                                  Tseq,
                                  pO2seq,
                                  # how to interpret "preys" / compute B
                                  prey_units = c("biomass", "prey_index", "feeding_level", "powerlaw"),
                                  preys      = NULL,         # biomass or f_ref; ignored for powerlaw
                                  intercepts = 1.6,          # powerlaw
                                  slopes     = -0.08,        # powerlaw
                                  m0         = 1,
                                  # pairing across masses×preys (or spectra)
                                  pair_mode  = c("auto","cross","zip"),
                                  normalize_E = TRUE,
                                  progress    = TRUE,
                                  # options passed to B_from_f
                                  T_ref = tr$T_ref, pO2_ref = 25,
                                  U_lo = 0, U_mech = 5, enforce_O2_feasible = TRUE,
                                  # scenario controls
                                  u_prey = 0, D = 3) {

  stopifnot(length(masses) > 0, length(Tseq) > 0, length(pO2seq) > 0)
  prey_units <- match.arg(prey_units)
  pair_mode  <- match.arg(pair_mode)
  if (!D %in% c(2,3)) stop("compute_state_surface(): D must be 2 or 3")

  # ---------- pairing helpers ----------
  make_spec_cross <- function() {
    if (prey_units %in% c("biomass","prey_index","feeding_level")) {
      if (is.null(preys)) stop("Provide 'preys' for biomass/prey_index/feeding_level modes.")
      tidyr::expand_grid(mass = masses, prey_arg = preys)
    } else {
      tidyr::expand_grid(mass = masses, intercept = intercepts, slope = slopes)
    }
  }
  make_spec_zip <- function() {
    if (prey_units %in% c("biomass","prey_index","feeding_level")) {
      if (is.null(preys)) stop("Provide 'preys' for biomass/prey_index/feeding_level modes.")
      mm <- masses; pp <- preys
      if (length(mm) == 1 && length(pp) > 1) mm <- rep(mm, length(pp))
      if (length(pp) == 1 && length(mm) > 1) pp <- rep(pp, length(mm))
      stopifnot(length(mm) == length(pp))
      tibble::tibble(mass = mm, prey_arg = pp)
    } else {
      mm <- masses; aa <- intercepts; bb <- slopes
      if (length(mm) == 1 && length(aa) > 1) mm <- rep(mm, length(aa))
      if (length(mm) == 1 && length(bb) > 1) mm <- rep(mm, length(bb))
      if (length(aa) == 1 && length(mm) > 1) aa <- rep(aa, length(mm))
      if (length(bb) == 1 && length(mm) > 1) bb <- rep(bb, length(mm))
      stopifnot(length(mm) == length(aa), length(mm) == length(bb))
      tibble::tibble(mass = mm, intercept = aa, slope = bb)
    }
  }
  if (pair_mode == "auto") {
    if (prey_units %in% c("biomass","prey_index","feeding_level") &&
        !is.null(preys) &&
        length(preys) == length(masses) && length(masses) > 1) {
      pair_mode <- "zip"
    } else if (prey_units == "powerlaw" &&
               (length(intercepts) == length(masses) ||
                length(slopes) == length(masses)) &&
               max(length(intercepts), length(slopes), length(masses)) > 1) {
      pair_mode <- "zip"
    } else {
      pair_mode <- "cross"
    }
  }
  spec <- if (pair_mode == "zip") make_spec_zip() else make_spec_cross()

  # ---------- resolve fixed-vs-dynamic B ----------
  resolve_B <- function(m, row) {
    if (prey_units == "biomass") {
      B  <- row$prey_arg
      meta <- list(mode = "biomass", dynamic = FALSE,
                   prey_B = B, f_ref = NA_real_,
                   pl_intercept = NA_real_, pl_slope = NA_real_, pl_m0 = NA_real_)
    } else if (prey_units == "prey_index") {
      f_ref <- row$prey_arg
      f_use <- pmin(pmax(f_ref, .Machine$double.eps), 1 - .Machine$double.eps)
      B  <- B_from_f(
        f_ref = f_use, m = m, tr = tr,
        T_ref = T_ref, pO2_ref = pO2_ref,
        U_lo = U_lo, U_mech = U_mech,
        enforce_O2_feasible = enforce_O2_feasible,
        u_prey = u_prey, D = D
      )
      meta <- list(mode = "prey_index", dynamic = FALSE,
                   prey_B = B, f_ref = f_ref,
                   pl_intercept = NA_real_, pl_slope = NA_real_, pl_m0 = NA_real_)
    } else if (prey_units == "feeding_level") {
      f_ref <- row$prey_arg
      meta <- list(mode = "feeding_level", dynamic = TRUE,
                   prey_B = NA_real_, f_ref = f_ref,
                   pl_intercept = NA_real_, pl_slope = NA_real_, pl_m0 = NA_real_)
      B <- NA_real_
    } else {
      a <- row$intercept; b <- row$slope
      B <- B_from_powerlaw_spectrum(m, intercept = a, slope = b, m0 = m0)
      meta <- list(mode = "powerlaw", dynamic = FALSE,
                   prey_B = B, f_ref = NA_real_,
                   pl_intercept = a, pl_slope = b, pl_m0 = m0)
    }
    list(B = B, meta = meta)
  }

  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(spec), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  # ---------- evaluate across T × pO2 ----------
  blocks <- vector("list", nrow(spec))
  for (i in seq_len(nrow(spec))) {
    m <- spec$mass[i]
    resB <- resolve_B(m, spec[i, , drop = FALSE])
    B_fixed <- resB$B
    meta    <- resB$meta

    grid <- tidyr::expand_grid(T = Tseq, pO2 = pO2seq)

    rr <- apply(as.matrix(grid), 1, function(row) {
      Tval <- as.numeric(row[[1]])
      pO2v <- as.numeric(row[[2]])

      # Choose B to use at this (Tval, pO2v)
      if (isTRUE(meta$dynamic) && identical(meta$mode, "feeding_level")) {
        f_ref <- meta$f_ref
        f_use <- pmin(pmax(f_ref, .Machine$double.eps), 1 - .Machine$double.eps)
        # B(T) satisfying f(T, pO2_ref) == f_ref (keep pO2_ref fixed)
        B_here <- B_from_f(
          f_ref = f_use, m = m, tr = tr,
          T_ref = Tval, pO2_ref = pO2_ref,
          U_lo = U_lo, U_mech = U_mech,
          enforce_O2_feasible = enforce_O2_feasible,
          u_prey = u_prey, D = D
        )
      } else {
        B_here <- B_fixed
      }

      st <- find_U_opt(pO2_env = pO2v, T = Tval, m = m, prey = B_here, tr = tr, u_prey = u_prey, D = D)

      Ustar <- as.numeric(st$U_opt)
      Enet  <- as.numeric(st$E_net)
      Cmax  <- as.numeric(st$Cmax)
      Cass  <- as.numeric(st$consump)
      Mm    <- as.numeric(st$M_m)
      Ma    <- as.numeric(st$M_act)
      f_eff <- as.numeric(st$f)

      Osup  <- suppressWarnings(as.numeric(st$O2_supply))
      Odem  <- suppressWarnings(as.numeric(st$O2_demand))
      Omarg <- suppressWarnings(as.numeric(st$O2_margin))

      if (!is.finite(Osup)) Osup <- Smax_whole(pO2v, Tval, m, tr)
      if (!is.finite(Odem) && is.finite(Ustar)) Odem <- O2_demand(Ustar, pO2v, Tval, m, B_here, tr, u_prey = u_prey, D = D)
      if (!is.finite(Omarg)) Omarg <- if (is.finite(Osup) && is.finite(Odem)) Osup - Odem else NA_real_

      c(B_used     = B_here,
        U_opt      = Ustar,
        E_net      = Enet,
        Cmax       = Cmax,
        consump    = Cass,
        M_m        = Mm,
        M_act      = Ma,
        f          = f_eff,
        O2_supply  = Osup,
        O2_demand  = Odem,
        O2_margin  = Omarg)
    })
    rr <- as.data.frame(t(rr))

    blk <- dplyr::bind_cols(
      tibble::tibble(
        mass = m,
        prey_unit     = meta$mode,
        prey_B        = meta$prey_B,
        f_ref         = meta$f_ref,
        pl_intercept  = meta$pl_intercept,
        pl_slope      = meta$pl_slope,
        pl_m0         = meta$pl_m0
      ),
      grid, rr
    )

    blocks[[i]] <- blk
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, i)
  }

  res <- dplyr::bind_rows(blocks)

  if (isTRUE(normalize_E) && all(c("E_net","Cmax") %in% names(res))) {
    res <- dplyr::mutate(
      res,
      E_net_norm = dplyr::if_else(is.finite(E_net) & is.finite(Cmax) & Cmax != 0, E_net / Cmax, NA_real_)
    )
  }

  tibble::as_tibble(res)
}


# -------------------------------------------------------------------
# Interpolate scattered (x,y,z) to a regular grid (for plotting)
# -------------------------------------------------------------------

interp_to_grid <- function(df, xcol, ycol, zcol, nx = 500, ny = 500) {
  stopifnot(all(c(xcol, ycol, zcol) %in% names(df)))

  x <- df[[xcol]]
  y <- df[[ycol]]
  z <- df[[zcol]]
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (!any(ok)) stop("No finite (x,y,z) to interpolate.")

  xo <- seq(min(x[ok]), max(x[ok]), length.out = nx)
  yo <- seq(min(y[ok]), max(y[ok]), length.out = ny)

  ip <- interp::interp(
    x[ok], y[ok], z[ok],
    xo = xo, yo = yo,
    linear = TRUE, extrap = FALSE
  )

  out <- expand.grid(stats::setNames(list(ip$x, ip$y), c(xcol, ycol)))
  out[[zcol]] <- as.vector(ip$z)

  out
}


# -------------------------------------------------------------------
# Regime labels (simple): oxygen limited when S ≈ D, prey limited otherwise; keep exclusions
# -------------------------------------------------------------------

add_regime_labels <- function(df, eps_O2 = 0.05, keep_helpers = FALSE) {
  base_needed <- c("O2_supply","O2_demand","U_opt","E_net")
  miss <- setdiff(base_needed, names(df))
  if (length(miss)) {
    stop("add_regime_labels(): missing columns: ", paste(miss, collapse = ", "))
  }

  out <- df |>
    dplyr::mutate(
      bad_S   = !is.finite(.data$O2_supply) | .data$O2_supply <= 0,
      bad_U   = !is.finite(.data$U_opt),
      bad_E   = !is.finite(.data$E_net) | .data$E_net <= 0,

      uO2     = .data$O2_demand / .data$O2_supply,
      near_O2 = is.finite(.data$uO2) & ((1 - .data$uO2) <= eps_O2),

      regime = dplyr::case_when(
        bad_U | bad_S | !is.finite(.data$uO2) ~ "Oxygen exclusion",
        bad_E                                ~ "Energetic exclusion",
        near_O2                              ~ "Oxygen limited",
        TRUE                                 ~ "Prey limited"
      ),
      regime   = factor(regime, levels = c("Oxygen exclusion","Energetic exclusion","Oxygen limited","Prey limited")),
      class_id = as.integer(regime)
    )

  if (!keep_helpers) {
    out <- dplyr::select(out, -bad_S, -bad_U, -bad_E, -uO2, -near_O2)
  }
  out
}
