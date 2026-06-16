# scripts/analysis_functions_ingestion.R
# -------------------------------------------------------------------
# Analysis helpers for internal O2 closure + behavioural optimization model.
# Uses B_from_f() that inverts the FULL model (optimized U, O2 closure)
# so that "prey_index" and "feeding_level" mean what you intended.
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
                                  w0 = 1,
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
      tibble::tibble(
        mass = rep_len(masses, max(length(masses), length(preys))),
        prey_arg = rep_len(preys,  max(length(masses), length(preys)))
      )
    } else {
      tibble::tibble(
        mass = rep_len(masses, max(length(masses), length(intercepts), length(slopes))),
        intercept = rep_len(intercepts, max(length(masses), length(intercepts), length(slopes))),
        slope = rep_len(slopes, max(length(masses), length(intercepts), length(slopes)))
      )
    }
  }
  
  if (pair_mode == "auto") {
    pair_mode <- if (length(masses) > 1 && !is.null(preys) && length(preys) == length(masses)) "zip" else "cross"
  }
  spec <- if (pair_mode == "zip") make_spec_zip() else make_spec_cross()
  
  # ----------------------- timing helpers -----------------------
  t0 <- Sys.time()
  fmt_elapsed <- function(t_start) {
    dt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    mins <- floor(dt / 60)
    secs <- round(dt - 60 * mins)
    sprintf("%dm %02ds", mins, secs)
  }
  
  # ----------------------- caching for B_from_f -----------------------
  cache <- new.env(parent = emptyenv())
  
  key_pi <- function(w, f_ref) {
    paste0("pi|m=", signif(w, 12),
           "|f=", signif(f_ref, 12),
           "|Tref=", signif(T_ref, 12),
           "|pO2ref=", signif(pO2_ref, 12),
           "|Ulo=", signif(U_lo, 12),
           "|Umech=", signif(U_mech, 12),
           "|enf=", enforce_O2_feasible,
           "|uprey=", signif(u_prey, 12),
           "|D=", D)
  }
  key_fl <- function(w, f_ref, Tval) {
    paste0("fl|m=", signif(w, 12),
           "|f=", signif(f_ref, 12),
           "|T=", signif(Tval, 12),
           "|pO2ref=", signif(pO2_ref, 12),
           "|Ulo=", signif(U_lo, 12),
           "|Umech=", signif(U_mech, 12),
           "|enf=", enforce_O2_feasible,
           "|uprey=", signif(u_prey, 12),
           "|D=", D)
  }
  
  get_B_prey_index <- function(w, f_ref) {
    k <- key_pi(w, f_ref)
    if (exists(k, envir = cache, inherits = FALSE)) return(get(k, envir = cache, inherits = FALSE))
    B <- B_from_f(f_ref = f_ref, w = w, tr = tr,
                  T_ref = T_ref, pO2_ref = pO2_ref,
                  U_lo = U_lo, U_mech = U_mech,
                  enforce_O2_feasible = enforce_O2_feasible,
                  u_prey = u_prey, D = D)
    assign(k, B, envir = cache)
    B
  }
  
  get_B_feeding_level <- function(w, f_ref, Tval) {
    k <- key_fl(w, f_ref, Tval)
    if (exists(k, envir = cache, inherits = FALSE)) return(get(k, envir = cache, inherits = FALSE))
    B <- B_from_f(f_ref = f_ref, w = w, tr = tr,
                  T_ref = Tval, pO2_ref = pO2_ref,
                  U_lo = U_lo, U_mech = U_mech,
                  enforce_O2_feasible = enforce_O2_feasible,
                  u_prey = u_prey, D = D)
    assign(k, B, envir = cache)
    B
  }
  
  # ----------------------- progress bar -----------------------
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(spec), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  blocks <- vector("list", nrow(spec))
  
  for (i in seq_len(nrow(spec))) {
    w <- spec$mass[i]
    grid <- tidyr::expand_grid(T = Tseq, pO2 = pO2seq)
    
    if (prey_units == "biomass") {
      B_mode <- "biomass"
      B_const <- spec$prey_arg[i]
      f_ref_out <- NA_real_
      B_lookup <- NULL
    } else if (prey_units %in% c("prey_index", "feeding_level")) {
      B_mode <- prey_units
      f_ref <- pmin(pmax(spec$prey_arg[i], .Machine$double.eps), 1 - .Machine$double.eps)
      f_ref_out <- spec$prey_arg[i]
      
      if (prey_units == "prey_index") {
        B_const <- get_B_prey_index(w, f_ref)
        B_lookup <- NULL
      } else {
        B_vec <- vapply(Tseq, function(Tval) get_B_feeding_level(w, f_ref, Tval), numeric(1))
        B_lookup <- B_vec
        B_const <- NA_real_
      }
    } else {
      B_mode <- "powerlaw"
      B_const <- B_from_powerlaw_spectrum(w, spec$intercept[i], spec$slope[i], w0 = w0)
      f_ref_out <- NA_real_
      B_lookup <- NULL
    }
    
    rr <- apply(as.matrix(grid), 1, function(row) {
      Tval <- as.numeric(row[[1]])
      pO2v <- as.numeric(row[[2]])
      
      B_used <- if (B_mode == "feeding_level") {
        idx <- match(Tval, Tseq)
        if (is.na(idx)) idx <- which.min(abs(Tseq - Tval))
        B_lookup[[idx]]
      } else {
        B_const
      }
      
      st <- find_U_opt(pO2_env = pO2v, T = Tval, w = w, prey = B_used, tr = tr,
                       u_prey = u_prey, D = D, U_lo = U_lo, U_mech = U_mech)
      
      c(B_used = B_used, f_ref = f_ref_out,
        U_opt = st$U_opt, E_net = st$E_net, Cmax = st$Cmax,
        Cmax_potential = st$Cmax_potential, Enc = st$Enc, C_pot = st$C_pot, C_real = st$C_real, I = st$I,
        f = st$f, g = st$g, pO2_int = st$pO2_int,
        A_assim = st$A_assim,
        M_m = st$M_m, M_act = st$M_act, D_SDA = st$D_SDA, M_exc = st$M_exc,
        O2_supply = st$O2_supply, O2_demand = st$O2_demand, O2_margin = st$O2_margin,
        oxygen_exclusion = as.numeric(st$oxygen_exclusion),
        energetic_exclusion = as.numeric(st$energetic_exclusion))
    })
    
    blk <- dplyr::bind_cols(tibble::tibble(mass = w), grid, as.data.frame(t(rr)))
    blocks[[i]] <- blk
    
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, i)
  }
  
  res <- dplyr::bind_rows(blocks)
  
  if (isTRUE(normalize_E)) {
    res <- dplyr::mutate(
      res,
      E_net_norm = dplyr::if_else(
        is.finite(.data$E_net) & is.finite(.data$Cmax) & .data$Cmax != 0,
        .data$E_net / .data$Cmax,
        NA_real_
      ),
      Cmax_frac = dplyr::if_else(
        is.finite(.data$Cmax) & is.finite(.data$Cmax_potential) & .data$Cmax_potential != 0,
        .data$Cmax / .data$Cmax_potential,
        NA_real_
      )
    )
  }
  
  if (isTRUE(progress)) {
    message(" compute_state_surface(): done in ", fmt_elapsed(t0))
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
  
  if (sum(ok) < 3) {
    out <- df[ok, c(xcol, ycol, zcol), drop = FALSE]
    return(out)
  }
  
  xo <- seq(min(x[ok]), max(x[ok]), length.out = nx)
  yo <- seq(min(y[ok]), max(y[ok]), length.out = ny)
  
  ip <- interp::interp(x[ok], y[ok], z[ok], xo = xo, yo = yo,
                       linear = TRUE, extrap = FALSE, duplicate = "mean")
  
  out <- expand.grid(stats::setNames(list(xo, yo), c(xcol, ycol)))
  zvec <- as.vector(ip$z)
  expected <- nrow(out)
  
  if (length(zvec) != expected) {
    return(df[ok, c(xcol, ycol, zcol), drop = FALSE])
  }
  
  out[[zcol]] <- zvec
  out
}