# scripts/model_functions_ingestion.R
# -------------------------------------------------------------------
# Bioenergetics fish model core functions (ingestion-decision form)
# Variant A:
#   nu_avail = (1-alpha_SDA-alpha_exc) * C_real - (M_M + M_A)
# where
#   f      = Enc / (Enc + Cmax)
#   C_pot  = f * Cmax
#   g      = Hill(pO2_int; K_g, h_g)
#   C_real = g * C_pot
#
# Internal O2 closure (Rubalcaba-inspired reduced form):
#   O2_supply = kO*(pO2_env - pO2_int)
#   O2_demand = M_M + M_A + alpha_SDA * C_real
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

  # Arrhenius preferred
  if (is.finite(E_x)) {
    T_K <- T + 273.15
    T_refK <- T_ref + 273.15
    return(exp(-(E_x / k_B) * (1 / T_K - 1 / T_refK)))
  }

  # Optional Q10 fallback
  if (!is.null(Q10_name)) {
    Q10_x <- get_par(tr, c(Q10_name), default = 1)
    return(Q10_x^((T - T_ref) / 10))
  }

  1
}

get_partition_fracs <- function(tr) {
  alpha_SDA <- get_par(tr, c("alpha_SDA", "a_SDA", "a_sda", "alpha"), default = 0.2)
  alpha_exc <- get_par(tr, c("alpha_exc", "a_exc", "eps_exc", "epsExc"), default = 0.1)

  alpha_SDA <- pmax(alpha_SDA, 0)
  alpha_exc <- pmax(alpha_exc, 0)
  alpha_assim <- pmax(1 - alpha_SDA - alpha_exc, 0)

  list(alpha_SDA = alpha_SDA, alpha_exc = alpha_exc, alpha_assim = alpha_assim)
}

v_rel_rms <- function(U, u_prey = 0) sqrt(U^2 + u_prey^2)
U_ref_from_mass <- function(w) 10^(-0.55) * (w^(1 / 6))

# NOTE: "whole" here returns mass-specific rates (per biomass).
Cmax_whole <- function(T, w, tr) {
  a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
  b_C <- get_par(tr, c("b_C", "b_c"), required = TRUE)
  F_C <- thermal_factor(T, tr, "E_C", "Q10_c")
  a_C * w^b_C * F_C / w
}

Mm_whole <- function(T, w, tr) {
  a_M <- get_par(tr, c("a_M", "a_m"), default = NA_real_)

  if (!is.finite(a_M)) {
    alpha_M <- get_par(tr, c("alpha_M", "f_m"), default = 0.1)
    a_C <- get_par(tr, c("a_C", "a_c"), required = TRUE)
    a_M <- alpha_M * a_C
  }

  b_M <- get_par(tr, c("b_M", "b_m"), default = 0.75)
  F_M <- thermal_factor(T, tr, "E_M", "Q10_m")
  a_M * w^b_M * F_M / w
}

Ma_whole <- function(U, T, w, tr) {
  a_A <- get_par(tr, c("a_A", "a_a"), required = TRUE)

  b_a <- get_par(tr, c("b_a", "b_A_mass"), default = NA_real_)
  b_U <- get_par(tr, c("b_U", "p_A", "b_A_speed"), default = NA_real_)

  if (!is.finite(b_a) || !is.finite(b_U)) {
    b_R <- get_par(tr, c("b_R"), default = 1)
    if (!is.finite(b_a)) b_a <- (2 - b_R) / 3
    if (!is.finite(b_U)) b_U <- 3 - b_R
  }

  F_A <- thermal_factor(T, tr, "E_A", "Q10_a")
  a_A * w^b_a * U^b_U * F_A / w
}

Enc_whole <- function(U, T, w, prey, tr, u_prey = 0, D = 3) {
  if (!D %in% c(2, 3)) stop("Enc_whole(): D must be 2 or 3")
  a_E <- get_par(tr, c("a_E", "a_e"), required = TRUE)
  b_E <- get_par(tr, c("b_E", "b_e"), required = TRUE)
  b_E_eff <- b_E * (D - 1) / 2
  a_E * w^b_E_eff * v_rel_rms(U^0.2, u_prey) * prey / w
}

hill_g <- function(p, K_g, h_g) {
  p_use <- pmax(p, 0)
  p_use^h_g / (p_use^h_g + K_g^h_g)
}

kO_whole <- function(T, w, tr) {
  a_O <- get_par(tr, c("a_O", "a_o"), required = TRUE)

  b_d <- get_par(tr, c("b_d"), default = 0.75)
  b_v <- get_par(tr, c("b_v"), default = 0.15)
  b_l <- get_par(tr, c("b_l"), default = 0.5 * b_d)

  b_O_default <- b_d + 0.5 * b_v - 0.5 * b_l
  b_O <- get_par(tr, c("b_O"), default = b_O_default)

  F_O <- thermal_factor(T, tr, "E_O", "Q10_o")
  a_O * w^b_O * F_O / w
}

solve_pO2_int <- function(U, pO2_env, T, w, prey, tr, u_prey = 0, D = 3) {
  fr <- get_partition_fracs(tr)
  alpha_SDA <- fr$alpha_SDA

  K_g <- get_par(tr, c("K_g", "K", "K_O2"), default = 2)
  h_g <- get_par(tr, c("h_g", "h", "h_O2"), default = 2)

  Cmax <- Cmax_whole(T, w, tr)
  Enc <- Enc_whole(U, T, w, prey, tr, u_prey = u_prey, D = D)

  f <- if (is.finite(Cmax) && Cmax > 0) Enc / (Enc + Cmax) else NA_real_
  C_pot <- Cmax * f

  Mm <- Mm_whole(T, w, tr)
  Ma <- Ma_whole(U, T, w, tr)
  kO <- kO_whole(T, w, tr)

  F <- function(p_int) {
    g <- hill_g(p_int, K_g = K_g, h_g = h_g)
    C_real <- g * C_pot
    kO * (pO2_env - p_int) - (Mm + Ma + alpha_SDA * C_real)
  }

  f0 <- F(0)
  f1 <- F(pO2_env)

  if (!is.finite(f0) || !is.finite(f1) || f0 < 0 || f1 > 0) {
    return(list(
      feasible = FALSE,
      pO2_int = NA_real_, g = NA_real_,
      C_pot = C_pot, C_real = NA_real_, f = f,
      Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
      alpha_assim = NA_real_, D_SDA = NA_real_, M_exc = NA_real_,
      O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_
    ))
  }

  p_int <- if (abs(f0) <= .TOL) 0 else uniroot(F, c(0, pO2_env))$root
  g <- hill_g(p_int, K_g = K_g, h_g = h_g)
  C_real <- g * C_pot

  alpha_assim <- fr$alpha_assim
  alpha_exc <- fr$alpha_exc

  nu_gain <- alpha_assim * C_real
  D_SDA <- alpha_SDA * C_real
  M_exc <- alpha_exc * C_real

  O2_supply <- kO * (pO2_env - p_int)
  O2_demand <- Mm + Ma + D_SDA

  list(
    feasible = TRUE,
    pO2_int = p_int, g = g,
    C_pot = C_pot, C_real = C_real, f = f,
    Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
    nu_gain = nu_gain, D_SDA = D_SDA, M_exc = M_exc,
    O2_supply = O2_supply, O2_demand = O2_demand, O2_margin = O2_supply - O2_demand
  )
}

nu_avail <- function(U, pO2_env, T, w, prey, tr, u_prey = 0, D = 3) {
  st <- solve_pO2_int(U, pO2_env, T, w, prey, tr, u_prey = u_prey, D = D)
  if (!isTRUE(st$feasible)) return(NA_real_)
  st$nu_gain - (st$Mm + st$Ma)
}

cap_U_by_O2_idle <- function(pO2_env, T, w, tr, U_lo = 1e-5, U_mech = 10) {
  kO <- kO_whole(T, w, tr)
  Mm <- Mm_whole(T, w, tr)

  F <- function(U) kO * pO2_env - (Mm + Ma_whole(U, T, w, tr))

  f_lo <- F(U_lo)
  f_hi <- F(U_mech)

  if (!is.finite(f_lo) || f_lo < 0) return(NA_real_)
  if (is.finite(f_hi) && f_hi >= 0) return(U_mech)
  uniroot(F, c(U_lo, U_mech))$root
}

find_U_opt <- function(pO2_env, T, w, prey, tr, u_prey = 0, D = 3,
                       U_lo = 1e-5, U_mech = 10) {

  U_cap <- cap_U_by_O2_idle(pO2_env, T, w, tr, U_lo = U_lo, U_mech = U_mech)
  if (!is.finite(U_cap)) {
    Mm <- Mm_whole(T, w, tr)
    Ma <- Ma_whole(U_lo, T, w, tr)
    return(list(
      M_m = Mm, M_act = Ma, Cmax = NA_real_, Enc = NA_real_,
      C_pot = NA_real_, C_real = NA_real_, I = NA_real_,
      f = NA_real_, g = NA_real_, pO2_int = NA_real_,
      nu_net = NA_real_, U_opt = NA_real_,
      O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_,
      oxygen_exclusion = TRUE, energetic_exclusion = NA,
      D_SDA = NA_real_, M_exc = NA_real_, nu_gain = NA_real_,
      consump = NA_real_, E_net = NA_real_, A_assim = NA_real_
    ))
  }

  obj <- function(U) nu_avail(U, pO2_env, T, w, prey, tr, u_prey = u_prey, D = D)
  opt <- optimize(obj, interval = c(U_lo, U_cap), maximum = TRUE)
  Ustar <- opt$maximum

  st <- solve_pO2_int(Ustar, pO2_env, T, w, prey, tr, u_prey = u_prey, D = D)
  nu_net <- if (isTRUE(st$feasible)) st$nu_gain - (st$Mm + st$Ma) else NA_real_

  list(
    M_m = st$Mm, M_act = st$Ma, Cmax = st$Cmax, Enc = st$Enc,
    C_pot = st$C_pot, C_real = st$C_real, I = st$C_real,
    f = st$f, g = st$g, pO2_int = st$pO2_int,
    nu_gain = st$nu_gain, nu_net = nu_net, U_opt = Ustar,
    O2_supply = st$O2_supply, O2_demand = st$O2_demand, O2_margin = st$O2_margin,
    oxygen_exclusion = FALSE,
    energetic_exclusion = is.finite(nu_net) && nu_net <= 0,
    D_SDA = st$D_SDA,
    M_exc = st$M_exc,
    # backward-compat aliases
    consump = st$nu_gain,
    E_net = nu_net,
    A_assim = st$nu_gain
  )
}

B_from_f <- function(f_ref, w, tr,
                     T_ref = get_par(tr, c("T_ref", "Tref"), default = 10),
                     pO2_ref = 25,
                     U_lo = 1e-5, U_mech = 5,
                     enforce_O2_feasible = TRUE,
                     eps_f = 1e-8,
                     u_prey = 0, D = 3,
                     B_lo = 1e-12,
                     B_hi = 1e12,
                     max_expand = 40) {

  n <- max(length(f_ref), length(w), length(T_ref), length(pO2_ref))
  f_vec <- rep_len(f_ref, n)
  w_vec <- rep_len(w, n)
  T_vec <- rep_len(T_ref, n)
  p_vec <- rep_len(pO2_ref, n)

  f_at_B <- function(B, pO2_env, T, w) {
    if (!is.finite(B) || B <= 0) return(0)
    st <- find_U_opt(pO2_env = pO2_env, T = T, w = w, prey = B, tr = tr,
                     u_prey = u_prey, D = D, U_lo = U_lo, U_mech = U_mech)
    if (isTRUE(st$oxygen_exclusion) || !is.finite(st$f)) return(NA_real_)
    st$f
  }

  out <- vapply(seq_len(n), function(i) {
    f_t <- f_vec[[i]]
    wi <- w_vec[[i]]
    Ti <- T_vec[[i]]
    pi <- p_vec[[i]]

    if (!is.finite(f_t) || f_t <= 0 || !is.finite(wi) || wi <= 0) return(0)
    f_t <- min(max(f_t, eps_f), 1 - eps_f)

    if (isTRUE(enforce_O2_feasible)) {
      U_cap <- cap_U_by_O2_idle(pi, Ti, wi, tr, U_lo = U_lo, U_mech = U_mech)
      if (!is.finite(U_cap)) return(NA_real_)
    }

    G <- function(log10B) {
      B <- 10^log10B
      fB <- f_at_B(B, pi, Ti, wi)
      if (!is.finite(fB)) return(NA_real_)
      fB - f_t
    }

    lo <- log10(B_lo)
    hi <- log10(B_hi)
    g_lo <- G(lo)
    g_hi <- G(hi)

    if (!is.finite(g_lo) || !is.finite(g_hi)) {
      U_guess <- U_ref_from_mass(wi)
      Cmax_ref <- Cmax_whole(Ti, wi, tr)
      Enc_needed <- (f_t / (1 - f_t)) * Cmax_ref

      b_E <- get_par(tr, c("b_E", "b_e"), required = TRUE)
      b_E_eff <- b_E * (D - 1) / 2
      a_E <- get_par(tr, c("a_E", "a_e"), required = TRUE)
      denom <- a_E * wi^b_E_eff * v_rel_rms(U_guess, u_prey = u_prey) / wi
      B_guess <- pmax(Enc_needed / denom, B_lo)

      lo <- log10(B_guess) - 6
      hi <- log10(B_guess) + 6
      g_lo <- G(lo)
      g_hi <- G(hi)
    }

    expand_iter <- 0
    while (is.finite(g_lo) && is.finite(g_hi) && sign(g_lo) == sign(g_hi) && expand_iter < max_expand) {
      lo <- lo - 1
      hi <- hi + 1
      g_lo <- G(lo)
      g_hi <- G(hi)
      expand_iter <- expand_iter + 1
    }

    if (!is.finite(g_lo) || !is.finite(g_hi) || sign(g_lo) == sign(g_hi)) {
      return(NA_real_)
    }

    root <- uniroot(function(x) G(x), c(lo, hi))$root
    10^root
  }, numeric(1))

  out
}

B_from_powerlaw_spectrum <- function(w, intercept = 1.6, slope = -0.08, w0 = 1) {
  pmax(intercept * (w / w0)^slope, 1e-12)
}

allocation_from_forced_f <- function(T, pO2_env, f_target, w, tr,
                                     U_choice = NULL, U_lo = 1e-5, U_mech = 5,
                                     enforce_O2_cap = TRUE,
                                     u_prey = 0, D = 3) {
  if (is.null(U_choice)) U_choice <- U_ref_from_mass(w)

  Cmax <- Cmax_whole(T, w, tr)
  EncReq <- (f_target / (1 - f_target)) * Cmax

  b_E <- get_par(tr, c("b_E", "b_e"), required = TRUE)
  b_E_eff <- b_E * (D - 1) / 2
  a_E <- get_par(tr, c("a_E", "a_e"), required = TRUE)

  U_use <- U_choice
  if (enforce_O2_cap) {
    U_cap <- cap_U_by_O2_idle(pO2_env, T, w, tr, U_lo = U_lo, U_mech = U_mech)
    if (is.finite(U_cap)) U_use <- min(U_use, U_cap)
  }

  denom <- a_E * w^b_E_eff * v_rel_rms(U_use, u_prey = u_prey) / w
  B_need <- EncReq / denom

  st <- find_U_opt(pO2_env, T, w, prey = B_need, tr, u_prey = u_prey, D = D,
                   U_lo = U_lo, U_mech = U_mech)

  list(M_m = st$M_m, M_act = st$M_act, nu_gain = st$nu_gain,
       nu_net = st$nu_net, U = U_use, prey = B_need, f = st$f,
       C_pot = st$C_pot, C_real = st$C_real, g = st$g, pO2_int = st$pO2_int,
       consump = st$consump, E_net = st$E_net)
}
