# scripts/model_functions_ingestion.R
# -------------------------------------------------------------------
# Bioenergetics fish model core functions (ingestion-decision form)
# Variant A (as requested):
#   E_avail = (1-a_SDA-a_exc) * C_real - (M_M + M_A)
# where
#   f      = Enc / (Enc + Cmax)
#   C_pot  = f * Cmax
#   g      = Hill(pO2_int; K, h)
#   C_real = g * C_pot
#
# Internal O2 closure (Rubalcaba-style):
#   O2_supply = kO*(pO2_env - pO2_int)
#   O2_demand = M_M + M_A + a_SDA * C_real        # excretion not respired
#   Solve O2_supply - O2_demand = 0 for pO2_int in [0, pO2_env]
#
# Behavioural optimization (Ware-style):
#   Choose U to maximize E_avail(U) subject to oxygen feasibility.
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

# Partition fractions on intake (Variant A):
# - a_SDA: fraction of C_real that is respired as processing cost (counts in O2 demand)
# - a_exc: fraction of C_real lost to excretion/egestion (NOT respired here)
# - a_assim: remaining fraction that becomes net assimilated energy
get_partition_fracs <- function(tr) {
  a_SDA <- get_par(tr, c("a_SDA", "a_sda", "alpha", "eps_SDA", "epsSDA"), default = 0.2)
  a_exc <- get_par(tr, c("a_exc", "eps_exc", "epsExc"), default = 0.1)
  a_SDA <- pmax(a_SDA, 0)
  a_exc <- pmax(a_exc, 0)
  a_assim <- pmax(1 - a_SDA - a_exc, 0)
  list(a_SDA = a_SDA, a_exc = a_exc, a_assim = a_assim)
}

v_rel_rms <- function(U, u_prey = 0) sqrt(U^2 + u_prey^2)
U_ref_from_mass <- function(m) 10^(-0.55) * (m^(1 / 6))

# NOTE: "whole" here returns mass-specific rates (per biomass), consistent with your existing code.
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

# Solve for internal oxygen pO2_int such that supply = demand.
# Returns f = Enc/(Enc+Cmax) (encounter saturation), g = Hill(pO2_int).
solve_pO2_int <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  fr <- get_partition_fracs(tr)
  a_SDA <- fr$a_SDA
  
  K <- get_par(tr, c("K", "K_O2"), default = 2)
  h <- get_par(tr, c("h", "h_O2"), default = 2)
  
  Cmax <- Cmax_whole(T, m, tr)
  Enc  <- Enc_whole(U, T, m, prey, tr, u_prey = u_prey, D = D)
  
  f <- if (is.finite(Cmax) && Cmax > 0) Enc / (Enc + Cmax) else NA_real_
  C_pot <- Cmax * f
  
  Mm <- Mm_whole(T, m, tr)
  Ma <- Ma_whole(U, T, m, tr)
  kO <- kO_whole(T, m, tr)
  
  F <- function(p_int) {
    g <- hill_g(p_int, K = K, h = h)
    C_real <- g * C_pot
    kO * (pO2_env - p_int) - (Mm + Ma + a_SDA * C_real)
  }
  
  f0 <- F(0)
  f1 <- F(pO2_env)  # supply=0; should be negative unless demand<=0
  
  # Feasibility at this U requires that even with max gradient (p_int=0),
  # oxygen supply can cover baseline costs (and any SDA at g(0)=0, which is 0).
  # Also require a sign change to bracket a root.
  if (!is.finite(f0) || !is.finite(f1) || f0 < 0 || f1 > 0) {
    return(list(
      feasible = FALSE,
      pO2_int = NA_real_, g = NA_real_,
      C_pot = C_pot, C_real = NA_real_, f = f,
      Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
      A_assim = NA_real_, D_SDA = NA_real_, M_exc = NA_real_,
      O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_
    ))
  }
  
  p_int <- if (abs(f0) <= .TOL) 0 else uniroot(F, c(0, pO2_env))$root
  g <- hill_g(p_int, K = K, h = h)
  C_real <- g * C_pot
  
  # Energetics (Variant A):
  A_assim <- fr$a_assim * C_real
  D_SDA   <- fr$a_SDA   * C_real  # respired processing cost
  M_exc   <- fr$a_exc   * C_real  # non-respired loss
  
  O2_supply <- kO * (pO2_env - p_int)
  O2_demand <- Mm + Ma + D_SDA
  
  list(
    feasible = TRUE,
    pO2_int = p_int, g = g,
    C_pot = C_pot, C_real = C_real, f = f,
    Mm = Mm, Ma = Ma, kO = kO, Cmax = Cmax, Enc = Enc,
    A_assim = A_assim, D_SDA = D_SDA, M_exc = M_exc,
    O2_supply = O2_supply, O2_demand = O2_demand, O2_margin = O2_supply - O2_demand
  )
}

# Net available energy for optimization
E_avail <- function(U, pO2_env, T, m, prey, tr, u_prey = 0, D = 3) {
  st <- solve_pO2_int(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  if (!isTRUE(st$feasible)) return(NA_real_)
  st$A_assim - (st$Mm + st$Ma)
}

# Oxygen feasibility cap on U ignoring feeding (i.e., g=0 => no SDA).
# Solve kO*pO2_env - (Mm + Ma(U)) = 0.
cap_U_by_O2_idle <- function(pO2_env, T, m, tr, U_lo = 1e-5, U_mech = 10) {
  kO <- kO_whole(T, m, tr)
  Mm <- Mm_whole(T, m, tr)
  
  F <- function(U) kO * pO2_env - (Mm + Ma_whole(U, T, m, tr))
  
  f_lo <- F(U_lo)
  f_hi <- F(U_mech)
  
  if (!is.finite(f_lo) || f_lo < 0) return(NA_real_)  # can't even cover baseline at minimal U
  if (is.finite(f_hi) && f_hi >= 0) return(U_mech)    # feasible up to mechanical cap
  uniroot(F, c(U_lo, U_mech))$root
}

# Optimize U to maximize E_avail(U) subject to O2 feasibility.
find_U_opt <- function(pO2_env, T, m, prey, tr, u_prey = 0, D = 3,
                       U_lo = 1e-5, U_mech = 10) {
  
  U_cap <- cap_U_by_O2_idle(pO2_env, T, m, tr, U_lo = U_lo, U_mech = U_mech)
  if (!is.finite(U_cap)) {
    Mm <- Mm_whole(T, m, tr)
    Ma <- Ma_whole(U_lo, T, m, tr)
    return(list(
      M_m = Mm, M_act = Ma, Cmax = NA_real_, Enc = NA_real_,
      C_pot = NA_real_, C_real = NA_real_, I = NA_real_,
      f = NA_real_, g = NA_real_, pO2_int = NA_real_,
      consump = NA_real_, E_net = NA_real_, U_opt = NA_real_,
      O2_supply = NA_real_, O2_demand = NA_real_, O2_margin = NA_real_,
      oxygen_exclusion = TRUE, energetic_exclusion = NA,
      D_SDA = NA_real_, M_exc = NA_real_, A_assim = NA_real_
    ))
  }
  
  obj <- function(U) E_avail(U, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  opt <- optimize(obj, interval = c(U_lo, U_cap), maximum = TRUE)
  Ustar <- opt$maximum
  
  st <- solve_pO2_int(Ustar, pO2_env, T, m, prey, tr, u_prey = u_prey, D = D)
  E_net <- if (isTRUE(st$feasible)) st$A_assim - (st$Mm + st$Ma) else NA_real_
  
  list(
    M_m = st$Mm, M_act = st$Ma, Cmax = st$Cmax, Enc = st$Enc,
    C_pot = st$C_pot, C_real = st$C_real, I = st$C_real,
    f = st$f, g = st$g, pO2_int = st$pO2_int,
    consump = st$A_assim, E_net = E_net, U_opt = Ustar,
    O2_supply = st$O2_supply, O2_demand = st$O2_demand, O2_margin = st$O2_margin,
    oxygen_exclusion = FALSE,
    energetic_exclusion = is.finite(E_net) && E_net <= 0,
    D_SDA = st$D_SDA,
    M_exc = st$M_exc,
    A_assim = st$A_assim
  )
}

# -------------------------------------------------------------------
# Prey standardization: invert the FULL bioenergetic-behavioural model.
#
# B_from_f is defined as:
#   Find prey biomass B such that, at (T_ref, pO2_ref),
#   the OPTIMIZED fish has encounter feeding level f equal to f_ref:
#       f(U*(B), B) = f_ref
# where U*(B) = argmax_U E_avail(U; B).
#
# This is what you described as "prey_index" (fixed reference) and
# "feeding_level" (allow T_ref to vary with T in surface computations).
# -------------------------------------------------------------------
B_from_f <- function(f_ref, m, tr,
                     T_ref = get_par(tr, c("T_ref", "Tref"), default = 10),
                     pO2_ref = 25,
                     U_lo = 1e-5, U_mech = 5,
                     enforce_O2_feasible = TRUE,
                     eps_f = 1e-8,
                     u_prey = 0, D = 3,
                     # numerical controls
                     B_lo = 1e-12,
                     B_hi = 1e12,
                     max_expand = 40) {
  
  n <- max(length(f_ref), length(m), length(T_ref), length(pO2_ref))
  f_vec <- rep_len(f_ref, n)
  m_vec <- rep_len(m, n)
  T_vec <- rep_len(T_ref, n)
  p_vec <- rep_len(pO2_ref, n)
  
  # Evaluate f at the OPTIMAL U for a candidate prey biomass B
  f_at_B <- function(B, pO2_env, T, m) {
    if (!is.finite(B) || B <= 0) return(0)
    st <- find_U_opt(pO2_env = pO2_env, T = T, m = m, prey = B, tr = tr,
                     u_prey = u_prey, D = D, U_lo = U_lo, U_mech = U_mech)
    if (isTRUE(st$oxygen_exclusion) || !is.finite(st$f)) return(NA_real_)
    st$f
  }
  
  out <- vapply(seq_len(n), function(i) {
    f_t <- f_vec[[i]]
    mi  <- m_vec[[i]]
    Ti  <- T_vec[[i]]
    pi  <- p_vec[[i]]
    
    if (!is.finite(f_t) || f_t <= 0 || !is.finite(mi) || mi <= 0) return(0)
    f_t <- min(max(f_t, eps_f), 1 - eps_f)
    
    if (isTRUE(enforce_O2_feasible)) {
      U_cap <- cap_U_by_O2_idle(pi, Ti, mi, tr, U_lo = U_lo, U_mech = U_mech)
      if (!is.finite(U_cap)) return(NA_real_)
    }
    
    # Root-find in log10(B) for stability
    G <- function(log10B) {
      B <- 10^log10B
      fB <- f_at_B(B, pi, Ti, mi)
      if (!is.finite(fB)) return(NA_real_)
      fB - f_t
    }
    
    lo <- log10(B_lo); hi <- log10(B_hi)
    g_lo <- G(lo); g_hi <- G(hi)
    
    # If evaluation fails at extremes, seed around a crude analytic guess
    if (!is.finite(g_lo) || !is.finite(g_hi)) {
      U_guess <- U_ref_from_mass(mi)
      Cmax_ref <- Cmax_whole(Ti, mi, tr)
      Enc_needed <- (f_t / (1 - f_t)) * Cmax_ref
      
      b_e <- get_par(tr, c("b_e"), required = TRUE)
      b_e_eff <- b_e * (D - 1) / 2
      a_e <- get_par(tr, c("a_e"), required = TRUE)
      denom <- a_e * mi^b_e_eff * v_rel_rms(U_guess, u_prey = u_prey) / mi
      B_guess <- pmax(Enc_needed / denom, B_lo)
      
      lo <- log10(B_guess) - 6
      hi <- log10(B_guess) + 6
      g_lo <- G(lo); g_hi <- G(hi)
    }
    
    # Expand bracket if needed
    expand_iter <- 0
    while (is.finite(g_lo) && is.finite(g_hi) && sign(g_lo) == sign(g_hi) && expand_iter < max_expand) {
      lo <- lo - 1
      hi <- hi + 1
      g_lo <- G(lo); g_hi <- G(hi)
      expand_iter <- expand_iter + 1
    }
    
    if (!is.finite(g_lo) || !is.finite(g_hi) || sign(g_lo) == sign(g_hi)) {
      return(NA_real_)  # no bracket => no solution under this full model definition
    }
    
    root <- uniroot(function(x) G(x), c(lo, hi))$root
    10^root
  }, numeric(1))
  
  out
}

B_from_powerlaw_spectrum <- function(m, intercept = 1.6, slope = -0.08, m0 = 1) {
  pmax(intercept * (m / m0)^slope, 1e-12)
}

# Convenience: given a target f at (T,pO2) and a chosen U (optional),
# compute prey biomass needed and then evaluate the full optimized model at that prey.
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