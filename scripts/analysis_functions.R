# scripts/analysis_functions.R
# -------------------------------------------------------------------
# Bioenergetics fish model functions for analysis
# Dependencies:
# -------------------------------------------------------------------


library(interp)
library(dplyr)

#
# 
#

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
                                  U_lo = 1e-3, U_mech = 5, enforce_O2_feasible = TRUE) {
  
  stopifnot(length(masses) > 0, length(Tseq) > 0, length(pO2seq) > 0)
  prey_units <- match.arg(prey_units)
  pair_mode  <- match.arg(pair_mode)
  
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
      # OLD behavior: fixed B from f_ref at (T_ref, pO2_ref)
      f_ref <- row$prey_arg
      # allow edge cases f_ref in [0,1] by gentle clamp
      f_use <- pmin(pmax(f_ref, .Machine$double.eps), 1 - .Machine$double.eps)
      B  <- B_from_f(f_ref = f_use, m = m, tr = tr, T_ref = T_ref, pO2_ref = pO2_ref,
                     U_lo = U_lo, U_mech = U_mech, enforce_O2_feasible = enforce_O2_feasible)
      meta <- list(mode = "prey_index", dynamic = FALSE,
                   prey_B = B, f_ref = f_ref,
                   pl_intercept = NA_real_, pl_slope = NA_real_, pl_m0 = NA_real_)
    } else if (prey_units == "feeding_level") {
      # NEW behavior: dynamic B(T) so that f(T, pO2_ref) == f_ref at each T
      f_ref <- row$prey_arg
      meta <- list(mode = "feeding_level", dynamic = TRUE,
                   prey_B = NA_real_, f_ref = f_ref,
                   pl_intercept = NA_real_, pl_slope = NA_real_, pl_m0 = NA_real_)
      B <- NA_real_  # computed per row later
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
        B_here <- B_from_f(f_ref = f_use, m = m, tr = tr,
                           T_ref = Tval, pO2_ref = pO2_ref,
                           U_lo = U_lo, U_mech = U_mech, enforce_O2_feasible = enforce_O2_feasible)
      } else {
        B_here <- B_fixed
      }
      
      st <- find_U_opt(pO2_env = pO2v, T = Tval, m = m, prey = B_here, tr = tr)
      
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
      if (!is.finite(Odem) && is.finite(Ustar)) Odem <- O2_demand(Ustar, pO2v, Tval, m, B_here, tr)
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
        prey_B        = meta$prey_B,     # fixed-B if applicable; NA for dynamic feeding_level
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


#
#
#

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
  
  # Build output with desired column names
  out <- expand.grid(stats::setNames(list(ip$x, ip$y), c(xcol, ycol)))
  out[[zcol]] <- as.vector(ip$z)  # matches expand.grid order
  
  out
}


#
#
#

add_B_gradients <- function(surface,
                            tr,
                            method    = c("auto", "surface", "perturb"),
                            prey_axis = c("B", "f_ref"),
                            use_norm  = FALSE,
                            d_rel     = 0.05,
                            progress  = TRUE,
                            # args used if prey_axis == "f_ref"
                            T_ref     = tr$T_ref,
                            pO2_ref   = 25,
                            U_lo      = 1e-3,
                            U_mech    = 5,
                            enforce_O2_feasible = TRUE) {
  
  method    <- match.arg(method)
  prey_axis <- match.arg(prey_axis)
  
  # which columns represent "prey" and "energy"?
  prey_col <- if (prey_axis == "B") "B_used" else "f_ref"
  E_col    <- if (use_norm) "E_net_norm" else "E_net"
  
  needed <- c("mass", "T", "pO2", "U_opt", prey_col, E_col,
              "O2_supply", "O2_demand")
  miss   <- setdiff(needed, names(surface))
  if (length(miss)) {
    stop("add_B_gradients(): missing columns: ",
         paste(miss, collapse = ", "))
  }
  
  df <- surface
  
  # ----------------------------------------------------------
  # 1) Surface-based central differences along prey axis
  # ----------------------------------------------------------
  if (method %in% c("auto", "surface")) {
    df <- df %>%
      dplyr::group_by(.data$mass, .data$T, .data$pO2) %>%
      dplyr::arrange(.data[[prey_col]], .by_group = TRUE) %>%
      dplyr::mutate(
        nP   = dplyr::n(),
        prey = .data[[prey_col]],
        E    = .data[[E_col]],
        O    = .data$O2_demand / .data$O2_supply,
        
        dE_dP = dplyr::case_when(
          !is.finite(prey) | nP < 2 ~ NA_real_,
          dplyr::row_number() == 1 ~
            (dplyr::lead(E) - E) /
            (dplyr::lead(prey) - prey),
          dplyr::row_number() == nP ~
            (E - dplyr::lag(E)) /
            (prey - dplyr::lag(prey)),
          TRUE ~
            (dplyr::lead(E) - dplyr::lag(E)) /
            (dplyr::lead(prey) - dplyr::lag(prey))
        ),
        
        dU_dP = dplyr::case_when(
          !is.finite(prey) | nP < 2 ~ NA_real_,
          dplyr::row_number() == 1 ~
            (dplyr::lead(.data$U_opt) - .data$U_opt) /
            (dplyr::lead(prey) - prey),
          dplyr::row_number() == nP ~
            (.data$U_opt - dplyr::lag(.data$U_opt)) /
            (prey - dplyr::lag(prey)),
          TRUE ~
            (dplyr::lead(.data$U_opt) - dplyr::lag(.data$U_opt)) /
            (dplyr::lead(prey) - dplyr::lag(prey))
        ),
        
        dO_dP = dplyr::case_when(
          !is.finite(prey) | nP < 2 ~ NA_real_,
          dplyr::row_number() == 1 ~
            (dplyr::lead(O) - O) /
            (dplyr::lead(prey) - prey),
          dplyr::row_number() == nP ~
            (O - dplyr::lag(O)) /
            (prey - dplyr::lag(prey)),
          TRUE ~
            (dplyr::lead(O) - dplyr::lag(O)) /
            (dplyr::lead(prey) - dplyr::lag(prey))
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-nP, -prey, -E, -O)
  }
  
  # If user asked for surface-only, we’re done here
  if (method == "surface") {
    df <- df %>%
      dplyr::rename(
        dE_dB = dE_dP,
        dU_dB = dU_dP,
        dO_dB = dO_dP
      )
    return(df)
  }
  
  # ----------------------------------------------------------
  # 2) Perturb-based fallback (or full) using find_U_opt()
  # ----------------------------------------------------------
  # rows needing perturbation: either 'perturb' mode, or NA gradients from surface
  need_idx <- if (method == "perturb") {
    seq_len(nrow(df))
  } else {
    which(!is.finite(df$dE_dP) | !is.finite(df$dU_dP) | !is.finite(df$dO_dP))
  }
  
  if (length(need_idx) == 0) {
    df <- df %>%
      dplyr::rename(
        dE_dB = dE_dP,
        dU_dB = dU_dP,
        dO_dB = dO_dP
      )
    return(df)
  }
  
  # helper: get "E" from state list, respecting use_norm
  get_E_from_state <- function(st) {
    if (!use_norm) return(st$E_net)
    if (!is.finite(st$E_net) || !is.finite(st$Cmax) || st$Cmax == 0) return(NA_real_)
    st$E_net / st$Cmax
  }
  
  # ------ perturb in B-space (prey_axis = "B") -----------------------
  grad_B_perturb <- function(mass, T, pO2, B, tr, d_rel) {
    if (!is.finite(B) || B <= 0) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    B_hi <- B * (1 + d_rel)
    B_lo <- B * (1 - d_rel)
    if (B_lo <= 0) B_lo <- B * 0.5
    
    st_hi <- find_U_opt(pO2_env = pO2, T = T, m = mass, prey = B_hi, tr = tr)
    st_lo <- find_U_opt(pO2_env = pO2, T = T, m = mass, prey = B_lo, tr = tr)
    
    if (!is.finite(st_hi$U_opt) || !is.finite(st_lo$U_opt)) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    E_hi <- get_E_from_state(st_hi)
    E_lo <- get_E_from_state(st_lo)
    
    O_hi <- if (is.finite(st_hi$O2_supply) && st_hi$O2_supply > 0) {
      st_hi$O2_demand / st_hi$O2_supply
    } else NA_real_
    
    O_lo <- if (is.finite(st_lo$O2_supply) && st_lo$O2_supply > 0) {
      st_lo$O2_demand / st_lo$O2_supply
    } else NA_real_
    
    if (!is.finite(E_hi) || !is.finite(E_lo) ||
        !is.finite(O_hi) || !is.finite(O_lo)) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    dB <- B_hi - B_lo
    c(
      dE_dP = (E_hi - E_lo) / dB,
      dU_dP = (st_hi$U_opt - st_lo$U_opt) / dB,
      dO_dP = (O_hi - O_lo) / dB
    )
  }
  
  # ------ perturb in f_ref-space (prey_axis = "f_ref") ---------------
  grad_f_perturb <- function(mass, T, pO2, f_ref, tr, d_rel) {
    if (is.na(f_ref) || !is.finite(f_ref) || f_ref <= 0) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    # multiplicative perturb inside (0,1)
    f_hi <- min(f_ref * (1 + d_rel), 1 - 1e-8)
    f_lo <- max(f_ref * (1 - d_rel), 1e-8)
    
    B_hi <- B_from_f(f_ref = f_hi, m = mass, tr = tr,
                     T_ref = T_ref, pO2_ref = pO2_ref,
                     U_lo = U_lo, U_mech = U_mech,
                     enforce_O2_feasible = enforce_O2_feasible)
    B_lo <- B_from_f(f_ref = f_lo, m = mass, tr = tr,
                     T_ref = T_ref, pO2_ref = pO2_ref,
                     U_lo = U_lo, U_mech = U_mech,
                     enforce_O2_feasible = enforce_O2_feasible)
    
    if (!is.finite(B_hi) || !is.finite(B_lo) || B_hi <= 0 || B_lo <= 0) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    st_hi <- find_U_opt(pO2_env = pO2, T = T, m = mass, prey = B_hi, tr = tr)
    st_lo <- find_U_opt(pO2_env = pO2, T = T, m = mass, prey = B_lo, tr = tr)
    
    if (!is.finite(st_hi$U_opt) || !is.finite(st_lo$U_opt)) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    E_hi <- get_E_from_state(st_hi)
    E_lo <- get_E_from_state(st_lo)
    
    O_hi <- if (is.finite(st_hi$O2_supply) && st_hi$O2_supply > 0) {
      st_hi$O2_demand / st_hi$O2_supply
    } else NA_real_
    
    O_lo <- if (is.finite(st_lo$O2_supply) && st_lo$O2_supply > 0) {
      st_lo$O2_demand / st_lo$O2_supply
    } else NA_real_
    
    if (!is.finite(E_hi) || !is.finite(E_lo) ||
        !is.finite(O_hi) || !is.finite(O_lo)) {
      return(c(dE_dP = NA_real_, dU_dP = NA_real_, dO_dP = NA_real_))
    }
    
    df_ <- f_hi - f_lo
    c(
      dE_dP = (E_hi - E_lo) / df_,
      dU_dP = (st_hi$U_opt - st_lo$U_opt) / df_,
      dO_dP = (O_hi - O_lo) / df_
    )
  }
  
  # choose which perturbation kernel to use
  perturb_fun <- if (prey_axis == "B") grad_B_perturb else grad_f_perturb
  
  # progress bar for perturbation step
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = length(need_idx), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  grads <- vector("list", length(need_idx))
  for (k in seq_along(need_idx)) {
    i   <- need_idx[k]
    row <- df[i, ]
    
    if (prey_axis == "B") {
      g <- grad_B_perturb(
        mass  = row$mass,
        T     = row$T,
        pO2   = row$pO2,
        B     = row$B_used,
        tr    = tr,
        d_rel = d_rel
      )
    } else {
      g <- grad_f_perturb(
        mass  = row$mass,
        T     = row$T,
        pO2   = row$pO2,
        f_ref = row$f_ref,
        tr    = tr,
        d_rel = d_rel
      )
    }
    
    grads[[k]] <- c(row_id = i, g)
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, k)
  }
  
  grads <- do.call(rbind, grads)
  grads <- as.data.frame(grads)
  grads$row_id <- as.integer(grads$row_id)
  
  # plug perturb-based gradients back in
  df$dE_dP[grads$row_id] <- grads$dE_dP
  df$dU_dP[grads$row_id] <- grads$dU_dP
  df$dO_dP[grads$row_id] <- grads$dO_dP
  
  # final rename: keep old names for downstream code
  df <- df %>%
    dplyr::rename(
      dE_dB = dE_dP,
      dU_dB = dU_dP,
      dO_dB = dO_dP
    )
  
  df
}



#
#
#

add_regime_labels <- function(
    df,
    eps_O2     = 0.05,   # closeness of demand/supply to call it O2-limited
    dE_thr     = 0.1,    # threshold for O2–prey co-limited: dE_dB > dE_thr
    use_derivs = TRUE,   # if FALSE, ignore derivative-based regimes
    keep_helpers = FALSE
) {
  # --- base columns always required ---------------------------------
  base_needed <- c("O2_supply","O2_demand","U_opt","E_net","M_act","M_m","f")
  miss_base   <- setdiff(base_needed, names(df))
  if (length(miss_base)) {
    stop("add_regime_labels(): missing columns: ",
         paste(miss_base, collapse = ", "))
  }
  
  # --- derivative columns only needed if we’re using them -----------
  if (isTRUE(use_derivs)) {
    deriv_needed <- c("dE_dB","dU_dB")
    miss_deriv   <- setdiff(deriv_needed, names(df))
    if (length(miss_deriv)) {
      stop("add_regime_labels(): use_derivs = TRUE but missing columns: ",
           paste(miss_deriv, collapse = ", "))
    }
  }
  
  # --- helper flags used in both modes ------------------------------
  out <- df |>
    dplyr::mutate(
      bad_S   = !is.finite(.data$O2_supply) | .data$O2_supply <= 0,
      bad_U   = !is.finite(.data$U_opt),
      bad_E   = !is.finite(.data$E_net) | .data$E_net <= 0,
      
      uO2     = .data$O2_demand / .data$O2_supply,
      near_O2 = is.finite(.data$uO2) & ((1 - .data$uO2) <= eps_O2)
    )
  
  # --- SIMPLE MODE: no derivative-based regimes ---------------------
  if (!isTRUE(use_derivs)) {
    
    reg_levels_simple <- c(
      "Oxygen exclusion",
      "Energetic exclusion",
      "Oxygen limited",
      "Prey limited"
    )
    
    out <- out |>
      dplyr::mutate(
        regime = dplyr::case_when(
          bad_U                               ~ "Oxygen exclusion",
          bad_S | !is.finite(.data$uO2)       ~ "Oxygen exclusion",
          bad_E                               ~ "Energetic exclusion",
          near_O2                             ~ "Oxygen limited",
          TRUE                                ~ "Prey limited"
        ),
        regime   = factor(regime, levels = reg_levels_simple),
        class_id = as.integer(regime)
      )
    
  } else {
    
    # --- DETAILED MODE: use dE_dB and dU_dB --------------------------
    reg_levels_detailed <- c(
      "Oxygen exclusion",
      "Energetic exclusion",
      "Oxygen limited (prey sensitive)",
      "Oxygen limited (strict)",
      "Prey limited",
      "Prey saturated"
    )
    
    out <- out |>
      dplyr::mutate(
        # co-limited interior band: O2-limited *and* strong prey leverage on E
        O2_prey_co = near_O2 & is.finite(.data$dE_dB) & (.data$dE_dB > dE_thr),
        
        regime = dplyr::case_when(
          # exclusion regimes
          bad_U                               ~ "Oxygen exclusion",
          bad_S | !is.finite(.data$uO2)       ~ "Oxygen exclusion",
          bad_E                               ~ "Energetic exclusion",
          
          # oxygen-limited band
          O2_prey_co                          ~ "Oxygen limited (prey sensitive)",
          near_O2                             ~ "Oxygen limited (strict)",
          
          # outside O2-limited: prey vs metabolic limitation
          !near_O2 & is.finite(.data$dU_dB) & (.data$dU_dB >  0) ~ "Prey limited",
          !near_O2 & is.finite(.data$dU_dB) & (.data$dU_dB <= 0) ~ "Prey saturated",
          
          # fallback (if derivatives are NA for some reason)
          TRUE                                ~ "Prey limited"
        ),
        regime   = factor(regime, levels = reg_levels_detailed),
        class_id = as.integer(regime)
      )
  }
  
  if (!keep_helpers) {
    out <- dplyr::select(out, -bad_S, -bad_U, -bad_E, -near_O2, -uO2,
                         dplyr::any_of("O2_prey_co"))
  }
  
  out
}


#
#
#

add_regime_labels2 <- function(
    df,
    eps_O2     = 0.05,   # closeness of demand/supply to call it O2-limited
    dE_thr     = 0.1,    # threshold for O2–prey co-limited: dE_dB > dE_thr
    dO_thr     = 0.01,   # small threshold on |dO_dB| for "tightening" vs "relaxing"
    use_derivs = TRUE,
    keep_helpers = FALSE
) {
  # --- base columns always required ---------------------------------
  base_needed <- c("O2_supply","O2_demand","U_opt","E_net","M_act","M_m","f")
  miss_base   <- setdiff(base_needed, names(df))
  if (length(miss_base)) {
    stop("add_regime_labels2(): missing columns: ",
         paste(miss_base, collapse = ", "))
  }
  
  # --- derivative columns only needed if we’re using them -----------
  if (isTRUE(use_derivs)) {
    deriv_needed <- c("dE_dB","dU_dB","dO_dB")
    miss_deriv   <- setdiff(deriv_needed, names(df))
    if (length(miss_deriv)) {
      stop("add_regime_labels2(): use_derivs = TRUE but missing columns: ",
           paste(miss_deriv, collapse = ", "))
    }
  }
  
  # --- helper flags used in both modes ------------------------------
  out <- df |>
    dplyr::mutate(
      bad_S   = !is.finite(.data$O2_supply) | .data$O2_supply <= 0,
      bad_U   = !is.finite(.data$U_opt),
      bad_E   = !is.finite(.data$E_net) | .data$E_net <= 0,
      
      uO2     = .data$O2_demand / .data$O2_supply,
      near_O2 = is.finite(.data$uO2) & ((1 - .data$uO2) <= eps_O2)
    )
  
  # --- SIMPLE MODE: no derivative-based regimes ---------------------
  if (!isTRUE(use_derivs)) {
    
    reg_levels_simple <- c(
      "Oxygen exclusion",
      "Energetic exclusion",
      "Oxygen limited",
      "Prey limited"
    )
    
    out <- out |>
      dplyr::mutate(
        regime = dplyr::case_when(
          bad_U                               ~ "Oxygen exclusion",
          bad_S | !is.finite(.data$uO2)       ~ "Oxygen exclusion",
          bad_E                               ~ "Energetic exclusion",
          near_O2                             ~ "Oxygen limited",
          TRUE                                ~ "Prey limited"
        ),
        regime   = factor(regime, levels = reg_levels_simple),
        class_id = as.integer(regime)
      )
    
  } else {
    
    # --- DETAILED MODE: use dE_dB (co-limited) + dO_dB (direction in non-O2 region) ---
    reg_levels_detailed <- c(
      "Oxygen exclusion",
      "Energetic exclusion",
      "Oxygen limited (prey sensitive)",
      "Oxygen limited (strict)",
      "Prey limited (O2-tightening)",
      "Prey limited (O2-relaxing)",
      "Prey limited (O2-neutral)"
    )
    
    out <- out |>
      dplyr::mutate(
        # co-limited interior band: O2-limited *and* strong prey leverage on E
        O2_prey_co = near_O2 & is.finite(.data$dE_dB) & (.data$dE_dB > dE_thr),
        
        # outside O2-limited: split by sign of dO_dB
        O2_tighten = !near_O2 & is.finite(.data$dO_dB) & (.data$dO_dB >  dO_thr),
        O2_relax   = !near_O2 & is.finite(.data$dO_dB) & (.data$dO_dB < -dO_thr),
        
        regime = dplyr::case_when(
          # exclusion regimes
          bad_U                               ~ "Oxygen exclusion",
          bad_S | !is.finite(.data$uO2)       ~ "Oxygen exclusion",
          bad_E                               ~ "Energetic exclusion",
          
          # oxygen-limited band
          O2_prey_co                          ~ "Oxygen limited (prey sensitive)",
          near_O2                             ~ "Oxygen limited (strict)",
          
          # non–O2-limited, dO splits direction
          O2_tighten                          ~ "Prey limited (O2-tightening)",
          O2_relax                            ~ "Prey limited (O2-relaxing)",
          
          # non–O2-limited but |dO| small or NA
          !near_O2                            ~ "Prey limited (O2-relaxing)",
          
          # fallback
          TRUE                                ~ "Prey limited (O2-relaxing)"
        ),
        regime   = factor(regime, levels = reg_levels_detailed),
        class_id = as.integer(regime)
      )
  }
  
  if (!keep_helpers) {
    out <- dplyr::select(
      out,
      -bad_S, -bad_U, -bad_E, -near_O2, -uO2,
      dplyr::any_of(c("O2_prey_co","O2_tighten","O2_relax"))
    )
  }
  
  out
}