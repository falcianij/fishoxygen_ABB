# scripts/model_functions.R
# -------------------------------------------------------------------
# Bioenergetics fish model core functions
# Minimal dependencies: base R only
#
# Patch (2026-01-26):
#   - Add RMS relative velocity: v_rel = sqrt(U^2 + u_prey^2)
#   - Add habitat dimensionality D in encounter scaling (D = 2 or 3)
#   - Thread u_prey and D through O2 demand, optimization, and prey-index calibrators
#   - NOTE: tr$b_e is interpreted as the 3D encounter exponent; 2D uses b_e/2
# -------------------------------------------------------------------

# Numerical tolerance (exported constant)
.TOL <- 1e-9


# ----- Relative speed (RMS) ------------------------------------------
v_rel_rms <- function(U, u_prey = 0) {
  sqrt(U^2 + u_prey^2)
}


# ----- Core rate functions (per-mass; m in grams) -------------------

Cmax_whole <- function(T, m, tr) {
  tr$a_c * m^tr$b_c * tr$Q10_c^((T - tr$T_ref)/10) / m
}

# Encounter (per-mass). Uses RMS relative speed and dimensionality-adjusted mass exponent.
Enc_whole <- function(U, T, m, prey, tr, u_prey = 0, D = 3) {
  if (!D %in% c(2, 3)) stop("Enc_whole(): D must be 2 or 3")

  # Interpret b_e as 3D exponent. Convert to D-dim via (D-1)/2.
  # If tr$b_e = 2/3 in 3D, then in 2D this becomes 1/3.
  b_e_eff <- tr$b_e * (D - 1) / 2

  vrel <- v_rel_rms(U, u_prey = u_prey)

  tr$a_e * m^b_e_eff * vrel * prey / m
}

Mm_whole <- function(T, m, tr) {
  tr$f_m * tr$a_c * m^tr$b_m * tr$Q10_m^((T - tr$T_ref)/10) / m
}

Ma_whole <- function(U, T, m, tr) {
  tr$a_a * m^((2 - tr$b_R)/3) * U^(3 - tr$b_R) * tr$Q10_a^((T - tr$T_ref)/10) / m
}

Smax_whole <- function(pO2_env, T, m, tr){
  (tr$a_g * m^(2/3) / tr$hm) * max(pO2_env, 1e-9) * tr$Q10_o^((T - tr$T_ref)/10) / m
}


# ----- Oxygen demand / objective ------------------------------------

O2_demand <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3){
  Mm   <- Mm_whole(T, m, tr)
  Ma   <- Ma_whole(U, T, m, tr)
  Cmax <- Cmax_whole(T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey, tr, u_prey = u_prey, D = D)
  f    <- Enc / (Enc + Cmax)
  Csda <- tr$alpha * f * Cmax
  Mm + Ma + Csda
}

E_avail <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3){
  Mm   <- Mm_whole(T, m, tr)
  Ma   <- Ma_whole(U, T, m, tr)
  Cmax <- Cmax_whole(T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey, tr, u_prey = u_prey, D = D)
  f    <- Enc / (Enc + Cmax)
  Cass <- tr$epsAssim * f * Cmax
  Cass - (Mm + Ma)
}


# ----- O2 capping & optimization ------------------------------------

cap_U_by_O2 <- function(pO2_env, T, m, prey, tr, U_lo, U_mech, u_prey = 0, D = 3){
  Smax <- Smax_whole(pO2_env, T, m, tr)
  F <- function(U) O2_demand(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D) - Smax
  f_lo <- F(U_lo);  f_hi <- F(U_mech)

  if (f_lo > 0) return(NA_real_)          # idling infeasible
  if (f_hi <= 0) return(U_mech)           # cap slack
  if (abs(f_lo) <= .TOL) return(U_lo)     # boundary

  uniroot(F, c(U_lo, U_mech))$root
}

O2_balance <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3){
  supply <- Smax_whole(pO2_env, T, m, tr)
  demand <- O2_demand(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  list(supply = supply, demand = demand, margin = supply - demand)
}

find_U_opt <- function(pO2_env, T, m, prey, tr, u_prey = 0, D = 3){
  U_lo   <- 0
  U_mech <- 10

  U_cap <- cap_U_by_O2(pO2_env, T, m, prey, tr, U_lo, U_mech, u_prey = u_prey, D = D)

  if (!is.finite(U_cap)) {
    Mm <- Mm_whole(T, m, tr)
    bal_idle <- O2_balance(U_lo, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
    return(list(
      M_m = Mm, M_act = NA_real_, Cmax = NA_real_, consump = NA_real_, E_net = NA_real_, U_opt = NA_real_, f = NA_real_,
      O2_supply = bal_idle$supply, O2_demand = bal_idle$demand, O2_margin = bal_idle$margin
    ))
  }

  if (U_cap <= U_lo + .TOL) {
    Ustar <- U_lo
    Mm   <- Mm_whole(T, m, tr)
    Ma   <- Ma_whole(Ustar, T, m, tr)
    Cmax <- Cmax_whole(T, m, tr)
    Enc  <- Enc_whole(Ustar, T, m, prey, tr, u_prey = u_prey, D = D)
    f    <- Enc / (Enc + Cmax)
    Cass <- tr$epsAssim * f * Cmax
    Enet <- Cass - (Mm + Ma)
    bal  <- O2_balance(Ustar, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
    return(list(
      M_m = Mm, M_act = Ma, Cmax = Cmax, consump = Cass, E_net = Enet, U_opt = Ustar, f = f,
      O2_supply = bal$supply, O2_demand = bal$demand, O2_margin = bal$margin
    ))
  }

  opt   <- optimize(function(U) E_avail(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D),
                    interval = c(U_lo, U_cap), maximum = TRUE)
  Ustar <- opt$maximum

  Mm      <- Mm_whole(T, m, tr)
  Ma      <- Ma_whole(Ustar, T, m, tr)
  Cmax    <- Cmax_whole(T, m, tr)
  Enc     <- Enc_whole(Ustar, T, m, prey, tr, u_prey = u_prey, D = D)
  f       <- Enc / (Enc + Cmax)
  Cassim  <- tr$epsAssim * f * Cmax
  Cunassim<- (1 - tr$epsAssim) * f * Cmax
  Enet    <- Cassim - (Mm + Ma)
  bal     <- O2_balance(Ustar, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)

  list(
    M_m = Mm, M_act = Ma, Cmax = Cmax, Enc = Enc, f = f,
    consump = Cassim, unassim = Cunassim, E_net = Enet, U_opt = Ustar,
    O2_supply = bal$supply, O2_demand = bal$demand, O2_margin = bal$margin
  )
}


# ----- Reference cruising speed --------------------------------------

U_ref_from_mass <- function(m) 10^(-0.55) * (m^(1/6))


# ----- Calibrators & spectra -----------------------------------------

B_from_f <- function(f_ref, m, tr,
                     T_ref = tr$T_ref, pO2_ref = 25,
                     U_lo = 0, U_mech = 5,
                     enforce_O2_feasible = TRUE,
                     eps_f = 1e-8,
                     u_prey = 0, D = 3) {
  # Vectorized over f_ref; m is assumed scalar (as in original use)
  # f_ref = 0  -> B = 0              (no encounter needed)
  # f_ref = 1  -> use 1 - eps_f      (asymptotic odds, finite and huge)
  # Returns NA if O2-feasible speed is not available even at U_lo.
  #
  # This version is consistent with Enc_whole(): uses RMS relative velocity and D-adjusted b_e.

  if (!is.finite(m) || m <= 0) stop("m must be a positive finite mass.")
  if (!is.finite(eps_f) || eps_f <= 0 || eps_f >= 0.1)
    stop("eps_f should be small and positive, e.g. 1e-8.")
  if (!D %in% c(2,3)) stop("B_from_f(): D must be 2 or 3")

  Uref <- U_ref_from_mass(m)
  Cmax_ref <- Cmax_whole(T_ref, m, tr)

  b_e_eff <- tr$b_e * (D - 1) / 2
  vref    <- v_rel_rms(Uref, u_prey = u_prey)

  one_B <- function(f) {
    if (!is.finite(f)) return(NA_real_)
    if (f <= 0) return(0)

    f_eff <- if (f >= 1) (1 - eps_f) else f
    Enc_needed <- (f_eff / (1 - f_eff)) * Cmax_ref

    denom_ref  <- tr$a_e * m^b_e_eff * vref / m
    Bf         <- Enc_needed / denom_ref
    if (!is.finite(denom_ref) || denom_ref <= 0) {
      # With zero relative speed (U=0 and u_prey=0), positive feeding level is unattainable.
      return(Inf)
    }

    if (enforce_O2_feasible) {
      U_cap <- cap_U_by_O2(pO2_ref, T_ref, m, Bf, tr, U_lo, U_mech, u_prey = u_prey, D = D)
      if (!is.finite(U_cap)) return(NA_real_)
      if (U_cap < Uref) {
        vcap <- v_rel_rms(U_cap, u_prey = u_prey)
        denom_cap <- tr$a_e * m^b_e_eff * vcap / m
        Bf <- Enc_needed / denom_cap
        if (!is.finite(denom_cap) || denom_cap <= 0) {
          return(Inf)
        }
      }
    }

    if (!is.finite(Bf)) return(NA_real_)
    if (Bf < 0) Bf <- 0
    Bf
  }

  vapply(f_ref, one_B, numeric(1))
}


# B_f(m) = intercept * (m / m0)^slope
B_from_powerlaw_spectrum <- function(m, intercept = 1.6, slope = -0.08, m0 = 1) {
  if (any(!is.finite(m)) || any(m <= 0)) stop("m must be positive and finite")
  if (!is.finite(intercept) || intercept <= 0) stop("intercept must be positive and finite")
  if (!is.finite(m0)  || m0 <= 0)  stop("m0 must be positive and finite")
  out <- intercept * (m / m0)^slope
  pmax(out, 1e-12)
}


allocation_from_forced_f <- function(T, pO2_env, f_target, m, tr,
                                     U_choice = NULL, U_lo = 0, U_mech = 5,
                                     enforce_O2_cap = TRUE,
                                     u_prey = 0, D = 3) {
  if (is.null(U_choice)) U_choice <- U_ref_from_mass(m)
  stopifnot(f_target > 0, f_target < 1)
  if (!D %in% c(2,3)) stop("allocation_from_forced_f(): D must be 2 or 3")

  Cmax   <- Cmax_whole(T, m, tr)
  EncReq <- (f_target / (1 - f_target)) * Cmax

  b_e_eff <- tr$b_e * (D - 1) / 2
  v_choice <- v_rel_rms(U_choice, u_prey = u_prey)
  denom  <- tr$a_e * m^b_e_eff * v_choice / m
  B_need <- EncReq / denom

  U_use <- U_choice
  if (enforce_O2_cap) {
    U_cap <- cap_U_by_O2(pO2_env, T, m, B_need, tr, U_lo, U_mech, u_prey = u_prey, D = D)
    if (is.finite(U_cap) && U_cap < U_choice) {
      U_use <- U_cap
      v_cap <- v_rel_rms(U_use, u_prey = u_prey)
      denom2 <- tr$a_e * m^b_e_eff * v_cap / m
      if (!is.finite(denom2) || denom2 <= 0) B_need <- Inf else B_need <- EncReq / denom2
    }
  }

  Mm   <- Mm_whole(T, m, tr)
  Ma   <- Ma_whole(U_use, T, m, tr)
  Cass <- tr$epsAssim * f_target * Cmax

  list(M_m = Mm, M_act = Ma, consump = Cass,
       E_net = Cass - (Mm + Ma), U = U_use, prey = B_need, f = f_target)
}
