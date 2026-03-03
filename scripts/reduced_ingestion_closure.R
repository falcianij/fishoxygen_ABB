# reduced_ingestion_closure.R
# -------------------------------------------------------------------
# Reduced fish bioenergetics + oxygen closure (no U optimization)
# Generated: 2026-02-06
#
# Goal
#   Replace adaptive behavior (optimize U) with:
#     f(B,pO2,T,m) = min( f_B(B,T,m; U_ref(m)),  f_O2(pO2,T,m) )
#
#   - f_B: Holling type-II encounter-limited feeding level evaluated at
#          reference cruising speed U_ref(m) = 10^(-0.55) * m^(1/6)
#   - f_O2: piecewise-linear ramp in pO2 between pO2_crit and pO2_limit:
#          f_O2 = 0 below pO2_crit
#          f_O2 = 1 above pO2_limit
#          linear in between
#     where:
#       pO2_crit  = 2.02 * m^(0.08) * Q10_p^((T-T_ref)/10)
#       pO2_limit = 3.7 * m^(0.08) * Q10_p^((T-T_ref)/10)
#     and Q10_p defaults to 1.33
#
# Notes
#   - Keeps your existing rate functions (Cmax, Mm, Ma, Smax, Enc) and uses
#     them at U = U_ref(m) for energetic accounting.
#   - Oxygen exclusion is diagnosed via idle-feasible movement at U_ref
#     ignoring SDA (movement feasibility):
#       Smax - (Mm + Ma(U_ref)) < 0  -> oxygen exclusion
#     (SDA oxygen cost is not used for exclusion because f_O2 already caps feeding.)
#
# Dependencies
#   - Base R only.
# -------------------------------------------------------------------

.TOL <- 1e-9
k_B <- 8.617e-5

get_par <- function(tr, names, default = NULL, required = FALSE) {
  for (nm in names) {
    if (!is.null(tr[[nm]]) && is.finite(tr[[nm]])) return(tr[[nm]])
  }
  if (required) stop("Missing required parameter (tr): ", paste(names, collapse = "/"))
  default
}

thermal_factor <- function(T, tr, E_name, Q10_name = NULL) {
  T_ref <- get_par(tr, c("T_ref", "Tref"), default = 10)
  E_x <- get_par(tr, c(E_name), default = NA_real_)
  if (is.finite(E_x)) {
    T_K <- T + 273.15
    T_refK <- T_ref + 273.15
    return(exp(-(E_x / k_B) * (1 / T_K - 1 / T_refK)))
  }
  if (!is.null(Q10_name)) {
    Q10_x <- get_par(tr, c(Q10_name), default = 1)
    return(Q10_x^((T - T_ref) / 10))
  }
  1
}


# ----- Relative speed (RMS) ------------------------------------------
v_rel_rms <- function(U, u_prey = 0) sqrt(U^2 + u_prey^2)

# ----- Reference cruising speed --------------------------------------
U_ref_from_mass <- function(w) 10^(-0.55) * (w^(1/6))

# ----- Core rate functions (per-mass; m in grams) ---------------------
Cmax_whole <- function(T, w, tr) {
  a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
  b_C <- get_par(tr, c("b_C", "b_c"), required = TRUE)
  a_C * w^b_C * thermal_factor(T, tr, "E_C", "Q10_c") / w
}

Mm_whole <- function(T, w, tr) {
  alpha_M <- get_par(tr, c("alpha_M", "f_m"), default = 0.1)
  a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
  b_M <- get_par(tr, c("b_M", "b_m"), default = 0.75)
  (alpha_M * a_C) * w^b_M * thermal_factor(T, tr, "E_M", "Q10_m") / w
}

Ma_whole <- function(U, T, w, tr) {
  a_A <- get_par(tr, c("a_A", "a_a"), required = TRUE)
  b_R <- get_par(tr, c("b_R"), default = 1)
  b_a <- get_par(tr, c("b_a"), default = (2 - b_R) / 3)
  b_U <- get_par(tr, c("b_U", "p_A"), default = 3 - b_R)
  a_A * w^b_a * U^b_U * thermal_factor(T, tr, "E_A", "Q10_a") / w
}

kO_whole <- function(T, w, tr) {
  a_O <- get_par(tr, c("a_O", "a_o"), required = TRUE)
  b_d <- get_par(tr, c("b_d"), default = 0.75)
  b_v <- get_par(tr, c("b_v"), default = 0.15)
  b_l <- get_par(tr, c("b_l"), default = 0.5 * b_d)
  b_O <- get_par(tr, c("b_O"), default = b_d + 0.5 * b_v - 0.5 * b_l)
  a_O * w^b_O * thermal_factor(T, tr, "E_O", "Q10_o") / w
}

Smax_whole <- function(pO2_env, T, w, tr){
  kO_whole(T, w, tr) * max(pO2_env, 1e-9)
}

# Encounter "capacity" (per-mass). Uses RMS relative speed and D-adjusted exponent.
Enc_whole <- function(U, T, w, prey, tr, u_prey = 0, D = 3) {
  if (!D %in% c(2, 3)) stop("Enc_whole(): D must be 2 or 3")
  b_E <- get_par(tr, c("b_E", "b_e"), required = TRUE)
  b_e_eff <- b_E * (D - 1) / 2
  vrel <- v_rel_rms(U, u_prey = u_prey)
  get_par(tr, c("a_E", "a_e"), required = TRUE) * w^b_e_eff * vrel * prey / w
}

# ----- Feeding-level closures ----------------------------------------

# Food-limited Holling-II feeding level evaluated at U_ref(m)
f_B_from_B <- function(B, T, m, tr, u_prey = 0, D = 3, U = NULL) {
  if (is.null(U)) U <- U_ref_from_mass(m)
  Cmax <- Cmax_whole(T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey = B, tr, u_prey = u_prey, D = D)
  if (!is.finite(Cmax) || Cmax <= 0) return(NA_real_)
  fB <- Enc / (Enc + Cmax)
  pmin(pmax(fB, 0), 1)
}

# Oxygen-limited feeding level as a ramp between (pO2_crit, pO2_limit)
pO2_crit_from_mT <- function(m, T, Q10_p = 1.33, T_ref = 10,
                            a_crit = 2.02, b = 0.08) {
  a_crit * m^b * Q10_p^((T - T_ref)/10)
}

pO2_limit_from_mT <- function(m, T, Q10_p = 1.33, T_ref = 10,
                             a_lim = 3.7, b = 0.08) {
  a_lim * m^b * Q10_p^((T - T_ref)/10)
}

f_O2_from_pO2 <- function(pO2, T, m, Q10_p = 1.33, T_ref = 10,
                          a_crit = 2.02, a_lim = 3.7, b = 0.08) {
  pcrit <- pO2_crit_from_mT(m, T, Q10_p = Q10_p, T_ref = T_ref, a_crit = a_crit, b = b)
  plim  <- pO2_limit_from_mT(m, T, Q10_p = Q10_p, T_ref = T_ref, a_lim  = a_lim,  b = b)

  if (!is.finite(pcrit) || !is.finite(plim)) return(NA_real_)
  if (plim <= pcrit + .TOL) return(ifelse(pO2 >= pcrit, 1, 0))

  ifelse(
    pO2 <= pcrit, 0,
    ifelse(pO2 >= plim, 1, (pO2 - pcrit) / (plim - pcrit))
  )
}

# Combined closure
f_reduced <- function(B, pO2, T, m, tr,
                      u_prey = 0, D = 3, U = NULL,
                      Q10_p = 1.33, T_ref_p = tr$T_ref,
                      a_crit = 2.02, a_lim = 3.7, b = 0.08) {

  fB  <- f_B_from_B(B, T, m, tr, u_prey = u_prey, D = D, U = U)
  fO2 <- f_O2_from_pO2(pO2, T, m, Q10_p = Q10_p, T_ref = T_ref_p,
                       a_crit = a_crit, a_lim = a_lim, b = b)

  if (!is.finite(fB) && !is.finite(fO2)) return(NA_real_)
  if (!is.finite(fB))  return(pmin(pmax(fO2, 0), 1))
  if (!is.finite(fO2)) return(pmin(pmax(fB, 0), 1))
  pmin(fB, fO2)
}

# ----- Single-point evaluation (B,pO2,T,m) ----------------------------

eval_state_reduced <- function(pO2_env, T, m, prey_B, tr,
                               u_prey = 0, D = 3,
                               Q10_p = 1.33, T_ref_p = tr$T_ref,
                               a_crit = 2.02, a_lim = 3.7, b = 0.08) {

  U <- U_ref_from_mass(m)

  Cmax <- Cmax_whole(T, m, tr)
  Mm   <- Mm_whole(T, m, tr)
  Ma   <- Ma_whole(U, T, m, tr)
  Smax <- Smax_whole(pO2_env, T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey = prey_B, tr, u_prey = u_prey, D = D)

  fB   <- if (is.finite(Cmax) && Cmax > 0) Enc/(Enc + Cmax) else NA_real_
  fO2  <- f_O2_from_pO2(pO2_env, T, m, Q10_p = Q10_p, T_ref = T_ref_p,
                        a_crit = a_crit, a_lim = a_lim, b = b)
  f    <- f_reduced(prey_B, pO2_env, T, m, tr,
                    u_prey = u_prey, D = D, U = U,
                    Q10_p = Q10_p, T_ref_p = T_ref_p,
                    a_crit = a_crit, a_lim = a_lim, b = b)

  I <- if (is.finite(Cmax) && Cmax > 0 && is.finite(f)) f * Cmax else NA_real_
  alpha_assim <- get_par(tr, c("alpha_assim", "epsAssim"), default = 0.7)
  nu_gain <- if (is.finite(I)) alpha_assim * I else NA_real_
  nu_net  <- if (is.finite(nu_gain)) nu_gain - (Mm + Ma) else NA_real_

  alpha <- get_par(tr, c("alpha_SDA", "alpha", "a_SDA", "a_sda"), default = 0.2)
  O2_sda <- if (is.finite(alpha) && alpha > 0 && is.finite(I)) alpha * I else 0
  O2_dem <- Mm + Ma + O2_sda
  O2_mar <- Smax - O2_dem

  # movement feasibility at U_ref with I=0 (no SDA)
  O2_margin_idle <- Smax - (Mm + Ma)

  pcrit <- pO2_crit_from_mT(m, T, Q10_p = Q10_p, T_ref = T_ref_p, a_crit = a_crit, b = b)
  plim  <- pO2_limit_from_mT(m, T, Q10_p = Q10_p, T_ref = T_ref_p, a_lim  = a_lim,  b = b)

  list(
    pO2 = pO2_env, T = T, mass = m, B_used = prey_B, D = D, u_prey = u_prey,
    U_ref = U,
    f = f, f_B = pmin(pmax(fB,0),1), f_O2 = pmin(pmax(fO2,0),1),
    pO2_crit = pcrit, pO2_limit = plim,
    Enc = Enc, Cmax = Cmax, I = I,
    M_m = Mm, M_act = Ma, nu_gain = nu_gain, nu_net = nu_net, consump = nu_gain, E_net = nu_net,
    O2_supply = Smax, O2_demand = O2_dem, O2_margin = O2_mar,
    O2_margin_idle = O2_margin_idle
  )
}

# ----- Surface computation (B vs pO2 at fixed T and mass) --------------
compute_state_surface_reduced <- function(tr,
                                         masses,
                                         Tseq,
                                         pO2seq,
                                         Bseq,
                                         u_prey = 0, D = 3,
                                         Q10_p = 1.33, T_ref_p = tr$T_ref,
                                         a_crit = 2.02, a_lim = 3.7, b = 0.08,
                                         normalize_nu = TRUE,
                                         progress = TRUE) {

  stopifnot(length(masses) > 0, length(Tseq) > 0, length(pO2seq) > 0, length(Bseq) > 0)
  if (!D %in% c(2,3)) stop("compute_state_surface_reduced(): D must be 2 or 3")

  grid_all <- expand.grid(
    mass = masses,
    T    = Tseq,
    pO2  = pO2seq,
    B    = Bseq,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  n <- nrow(grid_all)
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  out <- vector("list", n)

  for (i in seq_len(n)) {
    row <- grid_all[i, ]
    st <- eval_state_reduced(
      pO2_env = row$pO2,
      T       = row$T,
      m       = row$mass,
      prey_B  = row$B,
      tr      = tr,
      u_prey  = u_prey,
      D       = D,
      Q10_p   = Q10_p,
      T_ref_p = T_ref_p,
      a_crit  = a_crit,
      a_lim   = a_lim,
      b       = b
    )
    out[[i]] <- as.data.frame(st, stringsAsFactors = FALSE)

    if (isTRUE(progress) && (i %% 2000 == 0 || i == n)) utils::setTxtProgressBar(pb, i)
  }

  res <- do.call(rbind, out)

  if (isTRUE(normalize_nu)) {
    res$nu_net_norm <- ifelse(is.finite(res$nu_net) & is.finite(res$Cmax) & res$Cmax != 0,
                              res$nu_net / res$Cmax, NA_real_)
    res$E_net_norm <- ifelse(is.finite(res$E_net) & is.finite(res$Cmax) & res$Cmax != 0,
                             res$E_net / res$Cmax, NA_real_)
  }

  res
}

# ----- Regime labels (reduced closure) ---------------------------------
add_regime_labels_reduced <- function(df,
                                      eps_excl = 1e-6,
                                      eps_bind = 1e-6,
                                      eps_fsat = 0.99,
                                      keep_helpers = FALSE) {

  needed <- c("f","E_net","Cmax","O2_margin","O2_margin_idle")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))

  bad <- !is.finite(df$O2_margin_idle)

  oxy_excl  <- bad | (df$O2_margin_idle < -eps_excl)
  ener_excl <- !oxy_excl & is.finite(df$E_net) & (df$E_net <= 0)

  # approximate oxygen binding of this closed system via O2_margin with SDA included
  oxy_limited <- !oxy_excl & !ener_excl & is.finite(df$O2_margin) & (abs(df$O2_margin) <= eps_bind)

  sat <- !oxy_excl & !ener_excl &
    is.finite(df$Cmax) & df$Cmax > 0 &
    is.finite(df$f) & (df$f >= eps_fsat) &
    !oxy_limited

  enc_limited <- !oxy_excl & !ener_excl & !sat & !oxy_limited &
    is.finite(df$f) & (df$f > 0)

  regime <- ifelse(oxy_excl, "Oxygen exclusion",
            ifelse(ener_excl, "Energetic exclusion",
            ifelse(sat, "Saturated (Cmax-limited)",
            ifelse(oxy_limited, "Oxygen-limited feeding",
            ifelse(enc_limited, "Encounter-limited feeding",
                   "No feeding")))))

  if (isTRUE(keep_helpers)) {
    df$oxy_excl <- oxy_excl
    df$ener_excl <- ener_excl
    df$oxy_limited <- oxy_limited
    df$sat <- sat
    df$enc_limited <- enc_limited
  }

  df$regime <- regime
  df
}
