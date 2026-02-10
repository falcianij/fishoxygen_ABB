# scripts/model_functions_ingestion.R
# -------------------------------------------------------------------
# Bioenergetics fish model core functions (ingestion-decision version)
# Minimal dependencies: base R only
#
# Key changes (2026-01-27):
#   - Encounter provides a *potential* ingestion ceiling I_pot(U) (not forced ingestion)
#   - Realized ingestion I is chosen as I(U) = min( I_pot(U), I_O2(U) )
#       where I_O2(U) = max( Smax - Mm - Ma(U), 0 ) / alpha
#   - SDA oxygen demand depends on realized ingestion: O2_SDA = alpha * I
#   - Optimization is 1D over U (I collapses analytically)
#   - Relative velocity uses RMS: v_rel = sqrt(U^2 + u_prey^2)
#   - Habitat dimensionality D enters encounter via effective exponent:
#         b_e_eff = b_e_3D * (D-1)/2   (so 3D -> b_e, 2D -> b_e/2)
# -------------------------------------------------------------------

# Numerical tolerance
.TOL <- 1e-9

# ----- Relative speed (RMS) ------------------------------------------
v_rel_rms <- function(U, u_prey = 0) {
  sqrt(U^2 + u_prey^2)
}

# ----- Core rate functions (per-mass; m in grams) ---------------------

Cmax_whole <- function(T, m, tr) {
  tr$a_c * m^tr$b_c * tr$Q10_c^((T - tr$T_ref)/10) / m
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

# Encounter "capacity" (per-mass). Uses RMS relative speed and D-adjusted exponent.
Enc_whole <- function(U, T, m, prey, tr, u_prey = 0, D = 3) {
  if (!D %in% c(2, 3)) stop("Enc_whole(): D must be 2 or 3")

  # Interpret tr$b_e as 3D exponent; convert to D-dim via (D-1)/2.
  b_e_eff <- tr$b_e * (D - 1) / 2

  vrel <- v_rel_rms(U, u_prey = u_prey)
  tr$a_e * m^b_e_eff * vrel * prey / m
}

# ----- Potential ingestion ceiling from encounter ----------------------

# Potential feeding level (encounter-limited saturation)
f_pot_from_Enc <- function(Enc, Cmax) {
  Enc / (Enc + Cmax)
}

# Potential ingestion (per-mass) if fully engaged with prey, but still subject to Cmax.
I_pot_whole <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  Cmax <- Cmax_whole(T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey, tr, u_prey = u_prey, D = D)
  fpot <- f_pot_from_Enc(Enc, Cmax)
  list(I_pot = fpot * Cmax, f_pot = fpot, Enc = Enc, Cmax = Cmax)
}

# ----- Oxygen-limited ingestion ceiling --------------------------------

I_O2_cap_whole <- function(U, pO2_env, T, m, tr) {
  Smax <- Smax_whole(pO2_env, T, m, tr)
  Mm   <- Mm_whole(T, m, tr)
  Ma   <- Ma_whole(U, T, m, tr)
  list(Smax = Smax, Mm = Mm, Ma = Ma, margin0 = Smax - (Mm + Ma))
}

# If alpha <= 0 (no SDA oxygen cost), the oxygen ceiling on ingestion is effectively infinite.
I_O2_from_margin <- function(margin0, alpha) {
  if (!is.finite(alpha) || alpha <= 0) return(Inf)
  pmax(margin0, 0) / alpha
}

# ----- Objective and oxygen bookkeeping --------------------------------

# Realized ingestion given U (collapses the inner problem analytically)
I_star_given_U <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3,
                           tol_bind = 1e-10) {
  
  ip  <- I_pot_whole(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  cap <- I_O2_cap_whole(U, pO2_env, T, m, tr)
  
  # Oxygen ceiling on ingestion
  I_O2 <- I_O2_from_margin(cap$margin0, tr$alpha)
  
  # ---- SAFE min(): avoid NA poisoning ---------------------------------
  I_pot <- ip$I_pot
  if (!is.finite(I_pot) && !is.finite(I_O2)) {
    I <- NA_real_
  } else if (!is.finite(I_pot)) {
    I <- I_O2
  } else if (!is.finite(I_O2)) {
    I <- I_pot
  } else {
    I <- min(I_pot, I_O2)
  }
  
  # realized feeding level
  f_real <- if (is.finite(ip$Cmax) && ip$Cmax > 0 && is.finite(I)) I / ip$Cmax else NA_real_
  
  # ---- Oxygen bookkeeping ---------------------------------------------
  # If alpha <= 0, no SDA oxygen cost
  alpha <- tr$alpha
  O2_sda <- if (is.finite(alpha) && alpha > 0 && is.finite(I)) alpha * I else 0
  
  O2_dem <- cap$Mm + cap$Ma + O2_sda
  O2_mar <- cap$Smax - O2_dem
  
  # ---- Snap to the oxygen-binding manifold when O2 cap is active -------
  # If I is set by I_O2 (within relative tolerance), margin should be exactly 0.
  if (is.finite(I) && is.finite(I_O2) && is.finite(alpha) && alpha > 0) {
    if (abs(I - I_O2) <= tol_bind * (1 + abs(I_O2))) {
      O2_mar <- 0
      O2_dem <- cap$Smax
    }
  }
  
  list(
    I         = I,
    I_pot     = I_pot,
    I_O2      = I_O2,
    f         = f_real,
    f_pot     = ip$f_pot,
    Enc       = ip$Enc,
    Cmax      = ip$Cmax,
    Mm        = cap$Mm,
    Ma        = cap$Ma,
    Smax      = cap$Smax,
    O2_demand = O2_dem,
    O2_margin = O2_mar
  )
}


# Net energy available (same units as metabolic terms; uses epsAssim on realized ingestion)
E_avail <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  st <- I_star_given_U(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  # If baseline margin0 < 0, I_star_given_U will set I=0 but O2_margin may be negative;
  # E_avail is still defined, but such U should be excluded by cap_U_by_O2_idle.
  tr$epsAssim * st$I - (st$Mm + st$Ma)
}

# ----- O2 feasibility for movement (idle-feeding) ----------------------

# Maximum U such that breathing is feasible even with I = 0 (no SDA).
cap_U_by_O2_idle <- function(pO2_env, T, m, tr, U_lo = 0, U_mech = 10) {
  Smax <- Smax_whole(pO2_env, T, m, tr)
  Mm   <- Mm_whole(T, m, tr)

  # If even at U_lo (usually 0) baseline is infeasible, no movement is feasible.
  F <- function(U) (Mm + Ma_whole(U, T, m, tr)) - Smax
  f_lo <- F(U_lo); f_hi <- F(U_mech)

  if (f_lo > 0) return(NA_real_)     # oxygen exclusion (even at idle movement)
  if (f_hi <= 0) return(U_mech)      # feasible all the way to mechanical cap
  if (abs(f_lo) <= .TOL) return(U_lo)

  uniroot(F, c(U_lo, U_mech))$root
}

# ----- Main optimizer --------------------------------------------------

find_U_opt <- function(pO2_env, T, m, prey, tr, u_prey = 0, D = 3,
                       U_lo = 0, U_mech = 10) {
  if (!D %in% c(2,3)) stop("find_U_opt(): D must be 2 or 3")

  U_cap <- cap_U_by_O2_idle(pO2_env, T, m, tr, U_lo = U_lo, U_mech = U_mech)

  # Oxygen exclusion: even idle movement infeasible
  if (!is.finite(U_cap)) {
    # For context, compute baseline at U_lo with I forced to 0
    Mm   <- Mm_whole(T, m, tr)
    Ma   <- Ma_whole(U_lo, T, m, tr)
    Smax <- Smax_whole(pO2_env, T, m, tr)
    return(list(
      M_m = Mm, M_act = NA_real_, Cmax = NA_real_,
      Enc = NA_real_,
      I = 0, I_pot = NA_real_, I_O2 = NA_real_,
      f = 0, f_pot = NA_real_,
      consump = 0, E_net = NA_real_,
      O2_supply = Smax,
      O2_demand = Mm + Ma,
      O2_margin = Smax - (Mm + Ma),
      U_opt = NA_real_
    ))
  }

  # Optimize E_avail over feasible U
  opt <- optimize(
    f = function(U) E_avail(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D),
    interval = c(U_lo, U_cap),
    maximum = TRUE
  )
  Ustar <- opt$maximum

  st <- I_star_given_U(Ustar, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)

  Cassim <- tr$epsAssim * st$I
  Enet   <- Cassim - (st$Mm + st$Ma)

  list(
    M_m = st$Mm,
    M_act = st$Ma,
    Cmax = st$Cmax,
    Enc  = st$Enc,
    I    = st$I,
    I_pot = st$I_pot,
    I_O2  = st$I_O2,
    f    = st$f,
    f_pot = st$f_pot,
    consump = Cassim,
    E_net = Enet,
    U_opt = Ustar,
    O2_supply = st$Smax,
    O2_demand = st$O2_demand,
    O2_margin = st$O2_margin
  )
}

# ----- Reference cruising speed --------------------------------------

U_ref_from_mass <- function(m) 10^(-0.55) * (m^(1/6))

# ----- Calibrators & spectra -----------------------------------------

# "Prey index" calibration: find B needed to achieve a target *potential* feeding level f_ref
# at a reference environment (T_ref, pO2_ref) and reference speed Uref (or O2-capped U if requested).
B_from_f <- function(f_ref, m, tr,
                     T_ref = tr$T_ref, pO2_ref = 25,
                     U_lo = 0, U_mech = 5,
                     enforce_O2_feasible = TRUE,
                     eps_f = 1e-8,
                     u_prey = 0, D = 3) {

  if (!is.finite(m) || m <= 0) stop("m must be a positive finite mass.")
  if (!D %in% c(2,3)) stop("B_from_f(): D must be 2 or 3")
  if (!is.finite(eps_f) || eps_f <= 0 || eps_f >= 0.1)
    stop("eps_f should be small and positive, e.g. 1e-8.")

  Uref <- U_ref_from_mass(m)
  Cmax_ref <- Cmax_whole(T_ref, m, tr)

  one_B <- function(f) {
    if (!is.finite(f)) return(NA_real_)
    if (f <= 0) return(0)

    f_eff <- if (f >= 1) (1 - eps_f) else f

    # Required encounter (per mass) to achieve f_eff at reference conditions
    Enc_needed <- (f_eff / (1 - f_eff)) * Cmax_ref

    # Encounter coefficient at reference speed (using RMS vrel)
    b_e_eff <- tr$b_e * (D - 1) / 2
    vref    <- v_rel_rms(Uref, u_prey = u_prey)
    denom_ref <- tr$a_e * m^b_e_eff * vref / m

    Bf <- Enc_needed / denom_ref

    if (enforce_O2_feasible) {
      # Ensure the reference speed itself is feasible for idle movement (I=0)
      U_cap <- cap_U_by_O2_idle(pO2_ref, T_ref, m, tr, U_lo = U_lo, U_mech = U_mech)
      if (!is.finite(U_cap)) return(NA_real_)
      if (U_cap < Uref) {
        # Recompute B at the O2-feasible cap speed
        vcap <- v_rel_rms(U_cap, u_prey = u_prey)
        denom_cap <- tr$a_e * m^b_e_eff * vcap / m
        Bf <- Enc_needed / denom_cap
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

# Force a target *realized* feeding level f_target by choosing B; respects idle-movement O2 cap on U.
allocation_from_forced_f <- function(T, pO2_env, f_target, m, tr,
                                     U_choice = NULL, U_lo = 0, U_mech = 5,
                                     enforce_O2_cap = TRUE,
                                     u_prey = 0, D = 3) {
  if (!D %in% c(2,3)) stop("allocation_from_forced_f(): D must be 2 or 3")
  if (is.null(U_choice)) U_choice <- U_ref_from_mass(m)
  stopifnot(f_target > 0, f_target < 1)

  Cmax   <- Cmax_whole(T, m, tr)
  I_req  <- f_target * Cmax

  # Convert required ingestion to required *potential* encounter (since I <= I_pot)
  # Here we assume engagement is possible; i.e., the required f is a potential feeding level.
  EncReq <- (f_target / (1 - f_target)) * Cmax

  b_e_eff <- tr$b_e * (D - 1) / 2
  vch     <- v_rel_rms(U_choice, u_prey = u_prey)
  denom   <- tr$a_e * m^b_e_eff * vch / m
  B_need  <- EncReq / denom

  # enforce idle-movement oxygen cap on the chosen speed
  U_use <- U_choice
  if (enforce_O2_cap) {
    U_cap <- cap_U_by_O2_idle(pO2_env, T, m, tr, U_lo = U_lo, U_mech = U_mech)
    if (is.finite(U_cap) && U_cap < U_choice) {
      U_use <- U_cap
      vcap  <- v_rel_rms(U_use, u_prey = u_prey)
      denom2 <- tr$a_e * m^b_e_eff * vcap / m
      B_need <- EncReq / denom2
    }
  }

  # Evaluate realized outcome at (U_use, B_need) under ingestion-decision dynamics
  st <- find_U_opt(pO2_env, T, m, prey = B_need, tr, u_prey = u_prey, D = D, U_lo = U_lo, U_mech = U_mech)

  list(M_m = st$M_m, M_act = st$M_act, consump = st$consump,
       E_net = st$E_net, U = U_use, prey = B_need, f = st$f,
       I = st$I, I_pot = st$I_pot, I_O2 = st$I_O2)
}
