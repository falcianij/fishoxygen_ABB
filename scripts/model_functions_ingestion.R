# scripts/model_functions_ingestion.R
# -------------------------------------------------------------------
# Bioenergetics fish model core functions
# Rubalcaba-style internal O2 closure + Hill-limited processing,
# with Ware-style behavioural optimization over cruising speed U.
# -------------------------------------------------------------------

.TOL <- 1e-9
k_B <- 8.617e-5

`%||%` <- function(x, y) if (!is.null(x)) x else y

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

v_rel_rms <- function(U, u_prey = 0) sqrt(U^2 + u_prey^2)
U_ref_from_mass <- function(m) 10^(-0.55) * (m^(1 / 6))

Cmax_whole <- function(T, m, tr) {
  a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
  b_C <- get_par(tr, c("b_C", "b_c"), required = TRUE)
  F_C <- thermal_factor(T, tr, "E_C", "Q10_c")
  a_C * m^b_C * F_C / m
}

Mm_whole <- function(T, m, tr) {
  a_M <- get_par(tr, c("a_M"), default = NA_real_)
  if (!is.finite(a_M)) {
    f_m <- get_par(tr, c("f_m"), default = 0.1)
    a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
    a_M <- f_m * a_C
  }
  b_M <- get_par(tr, c("b_M", "b_m"), default = 0.75)
  F_M <- thermal_factor(T, tr, "E_M", "Q10_m")
  a_M * m^b_M * F_M / m
}

Ma_whole <- function(U, T, m, tr) {
  a_A <- get_par(tr, c("a_A", "a_a"), required = TRUE)
  b_A <- get_par(tr, c("b_A"), default = NA_real_)
  p_A <- get_par(tr, c("p_A"), default = NA_real_)
  if (!is.finite(b_A)) {
    b_R <- get_par(tr, c("b_R"), default = 1)
    b_A <- (2 - b_R) / 3
  }
  if (!is.finite(p_A)) {
    b_R <- get_par(tr, c("b_R"), default = 1)
    p_A <- 3 - b_R
  }
  F_A <- thermal_factor(T, tr, "E_A", "Q10_a")
  a_A * m^b_A * U^p_A * F_A / m
}

kO_whole <- function(T, m, tr) {
  a_O <- get_par(tr, c("a_O"), default = NA_real_)
  if (!is.finite(a_O)) {
    a_g <- get_par(tr, c("a_g"), default = 1)
    hm <- get_par(tr, c("hm"), default = 1)
    a_O <- a_g / hm
  }
  b_O <- get_par(tr, c("b_O"), default = 2 / 3)
  F_O <- thermal_factor(T, tr, "E_O", "Q10_o")
  a_O * m^b_O * F_O / m
}

Enc_whole <- function(U, T, m, prey, tr, u_prey = 0, D = 3) {
  if (!D %in% c(2, 3)) stop("Enc_whole(): D must be 2 or 3")
  a_e <- get_par(tr, c("a_e"), required = TRUE)
  b_e <- get_par(tr, c("b_e"), required = TRUE)
  b_e_eff <- b_e * (D - 1) / 2
  a_e * m^b_e_eff * v_rel_rms(U, u_prey) * prey / m
}

hill_g <- function(p, K, h) {
  p_use <- pmax(p, 0)
  p_use^h / (p_use^h + K^h)
}

solve_pO2_int <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  alpha <- get_par(tr, c("alpha"), default = 0)
  K <- get_par(tr, c("K", "K_O2"), default = 2)
  h <- get_par(tr, c("h", "h_O2"), default = 2)

  Cmax <- Cmax_whole(T, m, tr)
  Enc <- Enc_whole(U, T, m, prey, tr, u_prey = u_prey, D = D)
  f <- if (is.finite(Cmax) && Cmax > 0) Enc / (Enc + Cmax) else NA_real_
  C_pot <- Cmax * f

  Mm <- Mm_whole(T, m, tr)
  Ma <- Ma_whole(U, T, m, tr)
  kO <- kO_whole(T, m, tr)

  F <- function(p_int) {
    g <- hill_g(p_int, K = K, h = h)
    kO * (pO2_env - p_int) - (Mm + Ma + alpha * g * C_pot)
  }

  f0 <- F(0)
  f1 <- F(pO2_env)
  if (!is.finite(f0) || !is.finite(f1) || f0 < 0) {
    return(list(feasible = FALSE, pO2_int = NA_real_, g = NA_real_,
                C_pot = C_pot, C_real = NA_real_, f = f,
                Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
                O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_))
  }

  p_int <- if (abs(f0) <= .TOL) 0 else uniroot(F, c(0, pO2_env))$root
  g <- hill_g(p_int, K = K, h = h)
  C_real <- g * C_pot
  O2_supply <- kO * (pO2_env - p_int)
  O2_demand <- Mm + Ma + alpha * C_real

  list(feasible = TRUE, pO2_int = p_int, g = g,
       C_pot = C_pot, C_real = C_real, f = f,
       Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
       O2_supply = O2_supply, O2_demand = O2_demand,
       O2_margin = O2_supply - O2_demand)
}

E_avail <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  eps <- get_par(tr, c("eps", "epsAssim"), default = 0.7)
  st <- solve_pO2_int(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  if (!isTRUE(st$feasible)) return(NA_real_)
  eps * st$C_real - (st$Mm + st$Ma)
}

cap_U_by_O2_idle <- function(pO2_env, T, m, tr, U_lo = 1e-5, U_mech = 10) {
  kO <- kO_whole(T, m, tr)
  Mm <- Mm_whole(T, m, tr)
  F <- function(U) kO * pO2_env - (Mm + Ma_whole(U, T, m, tr))
  f_lo <- F(U_lo)
  f_hi <- F(U_mech)
  if (f_lo < 0) return(NA_real_)
  if (f_hi >= 0) return(U_mech)
  uniroot(F, c(U_lo, U_mech))$root
}

find_U_opt <- function(pO2_env, T, m, prey, tr, u_prey = 0, D = 3,
                       U_lo = 1e-5, U_mech = 10) {
  U_cap <- cap_U_by_O2_idle(pO2_env, T, m, tr, U_lo = U_lo, U_mech = U_mech)
  if (!is.finite(U_cap)) {
    Mm <- Mm_whole(T, m, tr)
    Ma <- Ma_whole(U_lo, T, m, tr)
    return(list(M_m = Mm, M_act = Ma, Cmax = NA_real_, Enc = NA_real_,
                C_pot = NA_real_, C_real = NA_real_, I = NA_real_,
                f = NA_real_, g = NA_real_, pO2_int = NA_real_,
                consump = NA_real_, E_net = NA_real_, U_opt = NA_real_,
                O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_,
                oxygen_exclusion = TRUE, energetic_exclusion = NA, D_SDA = NA_real_))
  }

  obj <- function(U) E_avail(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  opt <- optimize(obj, interval = c(U_lo, U_cap), maximum = TRUE)
  Ustar <- opt$maximum
  st <- solve_pO2_int(Ustar, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)

  eps <- get_par(tr, c("eps", "epsAssim"), default = 0.7)
  alpha <- get_par(tr, c("alpha"), default = 0)
  Cassim <- eps * st$C_real
  E_net <- Cassim - (st$Mm + st$Ma)

  list(M_m = st$Mm, M_act = st$Ma, Cmax = st$Cmax, Enc = st$Enc,
       C_pot = st$C_pot, C_real = st$C_real, I = st$C_real,
       f = st$f, g = st$g, pO2_int = st$pO2_int,
       consump = Cassim, E_net = E_net, U_opt = Ustar,
       O2_supply = st$O2_supply, O2_demand = st$O2_demand, O2_margin = st$O2_margin,
       oxygen_exclusion = FALSE,
       energetic_exclusion = is.finite(E_net) && E_net <= 0,
       D_SDA = alpha * st$C_real)
}

B_from_f <- function(f_ref, m, tr,
                     T_ref = get_par(tr, c("T_ref", "Tref"), default = 10), pO2_ref = 25,
                     U_lo = 1e-5, U_mech = 5,
                     enforce_O2_feasible = TRUE,
                     eps_f = 1e-8,
                     u_prey = 0, D = 3) {
  Uref <- U_ref_from_mass(m)
  Cmax_ref <- Cmax_whole(T_ref, m, tr)
  b_e <- get_par(tr, c("b_e"), required = TRUE)
  b_e_eff <- b_e * (D - 1) / 2
  a_e <- get_par(tr, c("a_e"), required = TRUE)

  vapply(f_ref, function(f) {
    if (!is.finite(f) || f <= 0) return(0)
    f_eff <- if (f >= 1) (1 - eps_f) else f
    Enc_needed <- (f_eff / (1 - f_eff)) * Cmax_ref

    U_use <- Uref
    if (enforce_O2_feasible) {
      U_cap <- cap_U_by_O2_idle(pO2_ref, T_ref, m, tr, U_lo = U_lo, U_mech = U_mech)
      if (!is.finite(U_cap)) return(NA_real_)
      U_use <- min(Uref, U_cap)
    }

    denom <- a_e * m^b_e_eff * v_rel_rms(U_use, u_prey = u_prey) / m
    pmax(Enc_needed / denom, 0)
  }, numeric(1))
}

B_from_powerlaw_spectrum <- function(m, intercept = 1.6, slope = -0.08, m0 = 1) {
  pmax(intercept * (m / m0)^slope, 1e-12)
}

allocation_from_forced_f <- function(T, pO2_env, f_target, m, tr,
                                     U_choice = NULL, U_lo = 1e-5, U_mech = 5,
                                     enforce_O2_cap = TRUE,
                                     u_prey = 0, D = 3) {
  if (is.null(U_choice)) U_choice <- U_ref_from_mass(m)
  Cmax <- Cmax_whole(T, m, tr)
  EncReq <- (f_target / (1 - f_target)) * Cmax
  b_e <- get_par(tr, c("b_e"), required = TRUE)
  b_e_eff <- b_e * (D - 1) / 2
  a_e <- get_par(tr, c("a_e"), required = TRUE)

  U_use <- U_choice
  if (enforce_O2_cap) {
    U_cap <- cap_U_by_O2_idle(pO2_env, T, m, tr, U_lo = U_lo, U_mech = U_mech)
    if (is.finite(U_cap)) U_use <- min(U_use, U_cap)
  }

  denom <- a_e * m^b_e_eff * v_rel_rms(U_use, u_prey = u_prey) / m
  B_need <- EncReq / denom
  st <- find_U_opt(pO2_env, T, m, prey = B_need, tr, u_prey = u_prey, D = D,
                   U_lo = U_lo, U_mech = U_mech)

  list(M_m = st$M_m, M_act = st$M_act, consump = st$consump,
       E_net = st$E_net, U = U_use, prey = B_need, f = st$f,
       C_pot = st$C_pot, C_real = st$C_real, g = st$g, pO2_int = st$pO2_int)
}
