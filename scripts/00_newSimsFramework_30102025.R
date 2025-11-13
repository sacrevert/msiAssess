# MSI + shrinkage + Freeman (R + JAGS)
# Sticky MNAR, multiple DGPs, empirical input (with/without SEs),
# convergence checks, and plotting.
#rm(list=ls())

suppressPackageStartupMessages({
  library(R2jags)
  library(splines)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# Utility ####
`%||%` <- function(x, y) if (is.null(x)) y else x
logit  <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))

# Measurement-scale helpers (loge | log10 | logit) for empirical inputs
.to_loge <- function(x, m.scale = c("loge","log10","logit")) {
  m.scale <- match.arg(m.scale)
  if (m.scale == "loge")  return(x)
  if (m.scale == "log10") return(x * log(10))
  if (m.scale == "logit") return(qlogis(x))
}
.from_loge <- function(x, m.scale = c("loge","log10","logit")) {
  m.scale <- match.arg(m.scale)
  if (m.scale == "loge")  return(x)
  if (m.scale == "log10") return(x / log(10))
  if (m.scale == "logit") return(plogis(x))
}

# 0) Common plotting theme ####
.msitheme <- function() theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank())

# 1) DGP generator for the common mean log-growth mu_t ####
# (DGP = "data generating process")
# Returns length n_years-1 vector (growth years).
# Non-congenial options added to avoid biasing in favour of splines.
generate_mu <- function(n_years,
                        mode = c("spline","rw","ar1","changepoint","seasonal","mixture"),
                        # spline
                        df_mu = 6,
                        # RW / AR1
                        sigma_eta = 0.05, phi = 0.7,
                        # change-points
                        n_cp = 2, jump_sd = 0.15, piecewise_linear = FALSE,
                        # seasonal
                        A = 0.15, period = 8, phi0 = runif(1, 0, 2*pi),
                        # mixture (K=2 by default)
                        K_guilds = 2) {
  mode <- match.arg(mode)
  T1 <- n_years - 1
  out_rw <- function() {
    mu <- numeric(T1); mu[1] <- rnorm(1, 0, sigma_eta)
    for (t in 2:T1) mu[t] <- mu[t-1] + rnorm(1, 0, sigma_eta)
    mu
  }
  out_ar1 <- function() {
    mu <- numeric(T1); mu[1] <- rnorm(1, 0, sigma_eta/sqrt(max(1e-8, 1-phi^2)))
    for (t in 2:T1) mu[t] <- phi*mu[t-1] + rnorm(1, 0, sigma_eta)
    mu
  }
  out_cp <- function() {
    mu <- rep(0, T1)
    cps <- if (n_cp>0) sort(sample(2:(T1-1), size = n_cp, replace = FALSE)) else integer(0)
    levels <- c(0, cumsum(rnorm(length(cps), 0, jump_sd)))
    seg <- cut(seq_len(T1), breaks = c(0, cps, T1), labels = FALSE, include.lowest = TRUE)
    mu_lvl <- levels[seg]
    if (!piecewise_linear) return(mu_lvl)
    mu_lin <- mu_lvl
    br <- c(1, cps, T1)
    for (j in seq_len(length(br)-1)) {
      a <- br[j]; b <- br[j+1]
      slope <- rnorm(1, 0, jump_sd/((b-a+1)/2))
      mu_lin[a:b] <- mu_lin[a:b] + slope * (seq(a, b) - a)
    }
    mu_lin
  }
  out_seas <- function() A * sin(2*pi*seq_len(T1)/period + phi0) + rnorm(T1, 0, sigma_eta)
  
  if (mode == "spline") {
    Z <- splines::ns(1:T1, df = df_mu)
    b <- rnorm(ncol(Z), 0, 0.15)
    beta0 <- rnorm(1, 0, 0.05)
    return(as.vector(beta0 + Z %*% b))
  }
  if (mode == "rw")          return(out_rw())
  if (mode == "ar1")         return(out_ar1())
  if (mode == "changepoint") return(out_cp())
  if (mode == "seasonal")    return(out_seas())
  if (mode == "mixture")     return(replicate(K_guilds, out_rw(), simplify = FALSE))
}

# 2) Simulation on the log scale with options + STICKINESS (state-dependent) ####
# innov_dist: "normal" or "student_t"
# sp_trend:   "none" or "random_slope" (per-species additive slope gamma_s)
# inclusion_bias: list(enabled=FALSE, a0, a1, rho1, rho0, p_init)
#   I_{s,t} ~ Bernoulli( logit^{-1}(a0 + a1 * r_{s,t-1} + rho0 + (rho1-rho0) I_{s,t-1}) )
#   r_{s,t-1} = mu_{g_s,t-1} + delta_{t-1} + gamma_s (guild-specific mean if mixture DGP)
# If dgp_mode = "mixture", species are split across K guilds with separate mu^{(k)}.

simulate_species_data <- function(n_species = 30, n_years = 30, seed = 232680,
                                  # DGP for the common mean growth
                                  dgp_mode = c("spline","rw","ar1","changepoint","seasonal","mixture"),
                                  df_mu = 6, sigma_eta = 0.05, phi = 0.7,
                                  n_cp = 2, jump_sd = 0.15, piecewise_linear = FALSE,
                                  A = 0.15, period = 8, phi0, K_guilds = 2,
                                  # species/state variation
                                  sigma_sp = 0.15, sigma_delta = 0.05, sd_alpha0 = 0.4,
                                  innov_dist = c("normal","student_t"), df_u = 3,
                                  sp_trend = c("none","random_slope"), sigma_gamma = 0.03, log_sd_se = 0.35,
                                  use_delta = TRUE,
                                  # observation model
                                  sigma_obs_mean = 0.25, prop_missing = 0.0,
                                  # inclusion MNAR
                                  inclusion_bias = list(enabled=FALSE, a0=-0.5, a1=3.0, rho1=1.0, rho0=0.0, p_init=1)) {
  set.seed(seed)
  dgp_mode  <- match.arg(dgp_mode)
  innov_dist<- match.arg(innov_dist)
  sp_trend  <- match.arg(sp_trend)
  years <- 1:n_years
  
  # Generate mean growth(s) at interval resolution t=1..T-1
  mu_obj <- generate_mu(n_years, mode = dgp_mode, df_mu = df_mu, sigma_eta = sigma_eta,
                        phi = phi, n_cp = n_cp, jump_sd = jump_sd, piecewise_linear = piecewise_linear,
                        A = A, period = period, phi0 = phi0, K_guilds = K_guilds)
  
  if (dgp_mode == "mixture") {
    mu_list <- mu_obj # list length K_guilds, each length T-1
    guild <- sample(seq_len(K_guilds), n_species, replace = TRUE)
    mu_for_species <- function(s) mu_list[[ guild[s] ]]
  } else {
    mu_true <- mu_obj # length T-1
    guild <- rep(1L, n_species)
    mu_for_species <- function(s) mu_true
  }
  
  # Common annual residuals (length T-1)
  delta_true <- if (use_delta) rnorm(n_years - 1, 0, sigma_delta) else rep(0, n_years - 1)
  
  # Species-specific additive slope (per-interval drift)
  gamma_s <- if (sp_trend == "random_slope") rnorm(n_species, 0, sigma_gamma) else rep(0, n_species)
  
  # Latent log-indices l_{s,t}
  l_true <- matrix(NA_real_, n_species, n_years)
  alpha <- rnorm(n_species, 0, sd_alpha0)
  l_true[,1] <- alpha
  
  # store the species-specific structural means mu_vec for later “truth” summaries
  MU_spec <- matrix(NA_real_, n_species, n_years - 1)
  
  for (t in 2:n_years) {
    # species-specific structural mean growth at interval t-1
    MU_spec[, t-1] <- vapply(seq_len(n_species), function(s) mu_for_species(s)[t-1], numeric(1))
    # process innovation
    u <- if (innov_dist == "normal") {
      rnorm(n_species, 0, sigma_sp)
    } else {
      u_raw <- rt(n_species, df = df_u)
      sd_t <- sqrt(df_u / max(2, df_u - 2))
      (u_raw / sd_t) * sigma_sp
    }
    # evolve
    l_true[,t] <- l_true[,t-1] + MU_spec[, t-1] + delta_true[t-1] + gamma_s + u
  }
  
  # Observations and SEs
  se <- matrix(sigma_obs_mean * exp(rnorm(n_species * n_years, 0, log_sd_se)), n_species, n_years)
  y  <- l_true + matrix(rnorm(n_species * n_years, 0, 1), n_species) * se
  
  # Inclusion (state-dependent stickiness)
  I <- matrix(1L, n_species, n_years)
  if (isTRUE(inclusion_bias$enabled)) {
    a0  <- inclusion_bias$a0 %||% -0.5
    a1  <- inclusion_bias$a1 %||%  3.0
    rho1<- inclusion_bias$rho1 %||% 1.0
    rho0<- inclusion_bias$rho0 %||% 0.0
    p1  <- inclusion_bias$p_init %||% 1
    I[,1] <- rbinom(n_species, 1, p1)
    for (t in 2:n_years) {
      r_det <- MU_spec[, t-1] + delta_true[t-1] + gamma_s
      eta <- a0 + a1 * r_det + rho0 + (rho1 - rho0) * I[,t-1]
      p   <- invlogit(eta)
      I[,t] <- rbinom(n_species, 1, p)
    }
    y[I == 0]  <- NA_real_
    se[I == 0] <- NA_real_
  }
  
  if (prop_missing > 0) {
    miss <- matrix(runif(n_species * n_years) < prop_missing, n_species, n_years)
    y[miss]  <- NA_real_;  se[miss] <- NA_real_
    if (isTRUE(inclusion_bias$enabled)) I[miss] <- 0L
  }
  
  ## Truth constructs ####
  # Realised species growths g_true (S x (T-1))
  g_true <- l_true[, 2:n_years, drop = FALSE] - l_true[, 1:(n_years-1), drop = FALSE]
  
  # “Common” (population-mean) generating growth: mean over species of MU_spec + delta
  mu_common_true <- colMeans(MU_spec) + delta_true # length T-1
  M_common_true  <- c(0, cumsum(mu_common_true)) # length T
  
  # Realised aggregate growth (mean of actual Delta_l) and its cumulative
  g_mean_true    <- colMeans(g_true) # length T-1
  M_realised_true<- c(0, cumsum(g_mean_true)) # length T
  
  # Selected truths under inclusion MNAR (pairwise presence for growth intervals)
  pair_mask <- I[, 2:n_years, drop = FALSE] & I[, 1:(n_years-1), drop = FALSE]
  denom <- colSums(pair_mask); denom[denom == 0] <- NA
  mu_selected_true <- colSums(g_true * pair_mask, na.rm = TRUE) / denom
  M_selected_true  <- c(0, cumsum(mu_selected_true))
  
  # For Bayesian geomean MNAR selected truth
  lmnar_selected_true <- ifelse(I == 1L, l_true, NA_real_)
  
  list(
    years = years, y = y, se = se, I = I,
    l_true = l_true,
    # generator truths
    mu_common_true = mu_common_true, # growth truth for plots (log scale)
    M_common_true  = M_common_true, # cumulative truth (log scale)
    # realised truths (optional overlays)
    g_mean_true    = g_mean_true,
    M_realised_true= M_realised_true,
    # selected-only truths (optional overlays for MNAR diagnostics)
    mu_selected_true = mu_selected_true, # for log-growth
    M_selected_true  = M_selected_true, # for cumulative M
    lmnar_selected_true = lmnar_selected_true, # for MSI
    # meta
    mu_true = if (dgp_mode=="mixture") NULL else mu_true %||% NA,
    mu_list = if (dgp_mode=="mixture") mu_list else NULL,
    guild = guild, delta_true = delta_true, gamma_s = gamma_s,
    innov_dist = innov_dist, df_u = df_u,
    inclusion_bias = inclusion_bias, use_delta = use_delta,
    meta = list(dgp_mode = dgp_mode, K_guilds = K_guilds)
  )
}


# 2b) Approximate calibration for a0 given target inclusion pi^* ####
# This could be used to set a0 in the MNAR sim setting, given a proportion of species
# that one wants included (i.e. there is an interaction between intercept and proportion missing)
# so, given an average species growth rate and other variation, there is a need to
# set a0 appropriately
approx_a0_for_target_pi <- function(mu_true, delta_true, gamma_s, a1, rho, pi_star) {
  r_bar <- mean(outer(mu_true + delta_true, rep(1, length(gamma_s))) + outer(rep(1, length(mu_true)), gamma_s))
  logit(pi_star) - a1 * r_bar - rho * pi_star
}

# 3) Basis constructor and JAGS data packing for partial/nopool models ####
make_basis <- function(years, df_mu = 6) ns(years[-1], df = df_mu)

# Enhanced data packer with observation mask for Bayesian geomean and others
prepare_jags_data2 <- function(sim, df_mu = 6, obs_var_model = 3) {
  Z <- make_basis(sim$years, df_mu)
  
  # keep ALL observed y
  idx <- which(!is.na(sim$y), arr.ind = TRUE)
  y_obs <- sim$y[idx]
  
  # SE handling: indicator + harmless filler
  se_mat <- sim$se
  has_se_mat <- is.finite(se_mat)
  obs_has_se <- as.integer(has_se_mat[idx])
  se_obs <- se_mat[idx]
  se_obs[obs_has_se == 0L] <- 1.0 # not used when obs_has_se==0
  
  # masks used elsewhere
  obs_mask   <- ifelse(!is.na(sim$y), 1L, 0L)
  n_present  <- as.integer(colSums(obs_mask))
  
  list(
    nsp = nrow(sim$y), nyears = ncol(sim$y), Z = Z, K = ncol(Z),
    n_obs = nrow(idx), y_obs = y_obs, se_obs = as.numeric(se_obs),
    obs_has_se = obs_has_se,
    sp_idx = idx[,1], t_idx = idx[,2],
    obs_var_model = as.integer(obs_var_model),
    obs_mask = obs_mask, n_present = n_present,
    eps = 1.0e-12
  )
}

# 4) JAGS models: obs_var switch for all models ####
# No observation with finite y is dropped from the likelihood
# where SEs are present, ovm=2/3 behaves as intended; where SEs
# are missing, those cells behave as ovm=4 (common theta)
obs_precision_block <- "
  w1 <- equals(obs_var_model, 1)
  w2 <- equals(obs_var_model, 2)
  w3 <- equals(obs_var_model, 3)
  w4 <- equals(obs_var_model, 4)

  c     ~ dunif(0, 10)
  theta ~ dunif(0, 5)

  for (i in 1:n_obs) {
    # ovm=1: effectively 'perfect' (debug)
    tau_w1[i] <- w1 * 1.0E6

    # ovm=2: use SE if present; else fall back to theta
    tau_w2[i] <- w2 * (
                  obs_has_se[i] * (1 / (pow(se_obs[i], 2) + eps)) +
                  (1 - obs_has_se[i]) * (1 / (pow(theta, 2) + eps))
                )

    # ovm=3: scaled SE if present; else fall back to theta
    tau_w3[i] <- w3 * (
                  obs_has_se[i] * (1 / (pow(c * se_obs[i], 2) + eps)) +
                  (1 - obs_has_se[i]) * (1 / (pow(theta, 2) + eps))
                )

    # ovm=4: common theta
    tau_w4[i] <- w4 * (1 / (pow(theta, 2) + eps))

    tau_i[i] <- tau_w1[i] + tau_w2[i] + tau_w3[i] + tau_w4[i]
    y_obs[i] ~ dnorm(l[sp_idx[i], t_idx[i]], tau_i[i])
  }
"

jags_partial <- paste0("
model {

  beta0 ~ dnorm(0, 1.0E-4)
  sigma_b ~ dunif(0, 10)
  tau_b <- pow(sigma_b, -2)
  
  for (k in 1:K) { b[k] ~ dnorm(0, tau_b) }
  
  sigma_delta ~ dunif(0, 2)
  tau_delta <- pow(sigma_delta, -2)
  sigma_sp ~ dunif(0, 5)
  tau_sp <- pow(sigma_sp, -2)
  
  for (t in 2:nyears) {
      mu[t-1] <- beta0 + inprod(b[], Z[t-1, 1:K])
      delta[t-1] ~ dnorm(0, tau_delta)
      mu_full[t-1] <- mu[t-1] + delta[t-1]
  }
    
  for (s in 1:nsp) {
    l[s,1] ~ dnorm(0, 1.0E-2)
    for (t in 2:nyears) {
      l[s,t] ~ dnorm(l[s,t-1] + mu_full[t-1], tau_sp)
    }
  }
",
obs_precision_block,
"
  M_smooth[1] <- 0
  M_full[1]   <- 0
  
  for (t in 2:nyears) {
    M_smooth[t] <- M_smooth[t-1] + mu[t-1]
    M_full[t]   <- M_full[t-1]   + mu_full[t-1]
    G_smooth[t-1] <- exp(mu[t-1])
    G_full[t-1]   <- exp(mu_full[t-1])
  }
  
  for (t in 1:nyears) {
    mean_logI[t] <- mean(l[1:nsp, t])
    MSI[t] <- exp(mean_logI[t] - mean_logI[1])
  }
  
}
")

jags_nopool <- paste0("
model {
  
  # # Shared sigma_b can induce pooling
  # sigma_b ~ dunif(0, 10)
  # tau_b <- pow(sigma_b, -2)
  # for (s in 1:nsp) {
  #   beta0[s] ~ dnorm(0, 1.0E-4)
  #   for (k in 1:K) { b[s,k] ~ dnorm(0, tau_b) }
  # }

  for (s in 1:nsp) {
    beta0[s] ~ dnorm(0, 1.0E-4)
    sigma_b[s] ~ dunif(0, 10)
    tau_b[s] <- pow(sigma_b[s], -2)
    for (k in 1:K) { 
    b[s,k] ~ dnorm(0, tau_b[s]) 
    }
  }

  # sigma_delta ~ dunif(0, 2)
  # tau_delta <- pow(sigma_delta, -2)

  # # Shared delta_t = common shock (can also lead to pooling)
  # for (t in 2:nyears) {
  #   delta[t-1] ~ dnorm(0, tau_delta)
  #   for (s in 1:nsp) {
  #     mu_s[s,t-1] <- beta0[s] + inprod(b[s,1:K], Z[t-1,1:K])
  #     mu_full_s[s,t-1] <- mu_s[s,t-1] + delta[t-1]
  #   }
  # }
  # 
  
  # No common year shock (purely species-specific smoothers)
  for (t in 2:nyears) {
    for (s in 1:nsp) {
      mu_s[s,t-1] <- beta0[s] + inprod(b[s,1:K], Z[t-1,1:K])
      mu_full_s[s,t-1] <- mu_s[s,t-1]
    }
  }
  
  # Species-specific process variance (no pooling via sigma_sp)
  for (s in 1:nsp) {
    sigma_sp[s] ~ dunif(0, 5)
    tau_sp[s] <- pow(sigma_sp[s], -2)
    l[s,1] ~ dnorm(0, 1.0E-2)
    for (t in 2:nyears) {
      l[s,t] ~ dnorm(l[s,t-1] + mu_full_s[s,t-1], tau_sp[s])
    }
  }
",
obs_precision_block,
"
  M_smooth[1] <- 0
  M_full[1]   <- 0
  
  for (t in 2:nyears) {
    mean_mu[t-1]      <- mean(mu_s[1:nsp, t-1]) # log growth (smooth)
    mean_mu_full[t-1] <- mean(mu_full_s[1:nsp, t-1]) # log growth inc. shocks
    M_smooth[t] <- M_smooth[t-1] + mean_mu[t-1]
    M_full[t]   <- M_full[t-1] + mean_mu_full[t-1]
    G_smooth[t-1] <- exp(mean_mu[t-1])
    G_full[t-1]   <- exp(mean_mu_full[t-1])
  }
  
  for (t in 1:nyears) {
    mean_logI[t] <- mean(l[1:nsp, t])
    MSI[t] <- exp(mean_logI[t] - mean_logI[1])
  }
  
}
")

# 4a) JAGS models: Bayesian geomean based on Soldaat et al. 2017
# species-wise independent; no time model, no imputation of missing data (i.e. sample-conditional)
jags_bayes_geomean <- paste0("
model {

  for (s in 1:nsp) {
  sigma_sp[s] ~ dunif(0, 5)
  tau_sp[s] <- pow(sigma_sp[s], -2)
    for (t in 1:nyears) {
      l[s,t] ~ dnorm(0, tau_sp[s])
    }
  }
",
obs_precision_block,
"
  # Geometric-mean index over observed species at each year
  for (t in 1:nyears) {
    for (s in 1:nsp) {
      l_eff[s,t] <- l[s,t] * obs_mask[s,t] # mask unobserved species
    }
    Lbar_num[t] <- sum(l_eff[1:nsp, t]) # single assignment per t
  
    # equals is equals(x,y), LOGICAL test for equality (i.e. returns 1 if true)
    denom[t] <- n_present[t] + equals(n_present[t], 0) # avoid division by zero
  
    Lbar[t] <- Lbar_num[t] / denom[t]
  }
  Lbar_base <- Lbar[1]
  for (t in 1:nyears) {
    MSI[t] <- exp(Lbar[t] - Lbar_base)
  }
  
}
"
)

# with rw1 + drift for imputation
jags_bayes_geomean_impute <- paste0("
model {

  tau_gamma ~ dgamma(2,200) # weak, sd(gamma_s) ~ 0.1
  
  # Species-wise local-level (RW1), no cross-species pooling
  for (s in 1:nsp) {
    sigma_sp[s] ~ dunif(0, 10) # wide, weakly-informative
    tau_sp[s] <- pow(sigma_sp[s], -2)
    
    gamma_s[s] ~ dnorm(0, tau_gamma)
    l[s,1] ~ dnorm(0, 1.0E-6)
    for (t in 2:nyears) {
      l[s,t] ~ dnorm(l[s,t-1] + gamma_s[s], tau_sp[s]) # RW1 increment + spp specific drift
      
      # If you prefer heavy tails
      # would require nu_rw argument to be specified in workflow
      #l[s,t] ~ dt(l[s,t-1], tau_sp[s], nu_rw)
    }
  }

",obs_precision_block,
"
  # All-species mean each year (NOT masked) and MSI baseline at t=1
  for (t in 1:nyears) {
    mean_logI[t] <- mean(l[1:nsp, t])
  }
  for (t in 1:nyears) {
    MSI[t] <- exp(mean_logI[t] - mean_logI[1])
  }
}
")

## Post-hoc ridge smoother for growth ####
# Choose lambda by generalized cross-validation (on a single y, e.g., median growth)
.choose_lambda_gcv <- function(y, Z, grid = exp(seq(log(1e-6), log(1e3), length.out = 60))) {
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  sv <- svd(Zc) # Zc = U D V'
  s2 <- sv$d^2 # eigenvalues of Zc'Zc
  n  <- nrow(Zc)
  gcv <- vapply(grid, function(lam) {
    # Smoother S = U diag(s2/(s2+lam)) U'
    Sy  <- sv$u %*% ( (s2/(s2+lam)) * crossprod(sv$u, yc) )
    yhat <- as.vector(Sy + mean(y))
    rss  <- sum((y - yhat)^2)
    trS  <- sum(s2/(s2+lam))
    rss / ((n - trS)^2 + 1e-12)
  }, numeric(1))
  grid[which.min(gcv)]
}

# Smooth an entire draws x time matrix of growth using a fixed lambda
.smooth_growth_draws <- function(Gmat, Z, lambda) {
  # Gmat: draws x (T-1), Z: (T-1) x K (no intercept), ridge on columns of Z
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  ZZ <- crossprod(Zc)
  P  <- solve(ZZ + diag(lambda, ncol(Zc)), tol = 1e-10)
  Bx <- P %*% t(Zc) # K x (T-1); precompute
  # For each draw: center y, ridge fit b, add back mean
  apply(Gmat, 1, function(y) {
    yc <- y - mean(y)
    b  <- Bx %*% yc
    as.vector(Zc %*% b + mean(y))
  }) |> t() # returns draws x (T-1)
}

# post-hoc smoothing for Bayes geomean
posthoc_smooth_bayes_geomean <- function(fit, years, impute_all,
                                         basis = c("auto","ruppert","ns"),
                                         df_mu = 6, num_knots = 12) {
  basis <- match.arg(basis)
  T <- length(years)
  sims <- fit$BUGSoutput$sims.list
  
  # 1) Choose source of the log-index
  Lbar <- if (isTRUE(impute_all)) sims$mean_logI else sims$Lbar
  if (is.null(Lbar)) {
    stop(sprintf("Required node %s not found in fit$BUGSoutput$sims.list.",
                 if (isTRUE(impute_all)) "'mean_logI'" else "'Lbar'"))
  }
  if (ncol(Lbar) != T) {
    stop(sprintf("Time dimension mismatch: index has %d columns but length(years) = %d.",
                 ncol(Lbar), T))
  }
  
  # 2) Growth draws (draws x (T-1))
  dL <- Lbar[, -1, drop = FALSE] - Lbar[, -T, drop = FALSE]
  
  # 3) Build growth-year basis Z ((T-1) x K), robustly
  Z <- NULL; basis_used <- NULL
  if (basis %in% c("auto","ruppert")) {
    # clamp knot count to available degrees of freedom
    nk_eff <- max(1L, min(num_knots, (T - 1) - 2L))
    Z_try <- try(build_ruppert_basis(t = 1:(T - 1), num_knots = nk_eff)$Z, silent = TRUE)
    if (!inherits(Z_try, "try-error")) { Z <- Z_try; basis_used <- sprintf("ruppert[%d]", nk_eff) }
  }
  if (is.null(Z)) {
    # clamp df for ns to be feasible
    df_eff <- max(1L, min(df_mu, (T - 1) - 1L))
    Z <- make_basis(years, df_mu = df_eff)
    basis_used <- sprintf("ns[df=%d]", df_eff)
  }
  
  # 4) Choose lambda on the posterior-median growth, then smooth all draws
  lam   <- .choose_lambda_gcv(apply(dL, 2, stats::median, na.rm = TRUE), Z)
  mu_sm <- .smooth_growth_draws(dL, Z, lam) # draws x (T-1)
  
  # 5) Re-integrate to cumulative log and rebase to t=1
  M_sm   <- cbind(0, t(apply(mu_sm, 1, cumsum))) # draws x T
  MSI_sm <- exp(M_sm - M_sm[, 1])
  
  list(
    Lambda = posterior_summary(mu_sm), # smoothed log-growth
    Gexp   = posterior_summary(exp(mu_sm)), # multiplicative growth
    M      = posterior_summary(M_sm), # cumulative log-growth
    MSI    = posterior_summary(MSI_sm), # index, baseline t=1
    lambda = lam,
    basis  = basis_used,
    impute_all = impute_all
  )
}

# 5) Freeman-conformant (BRC-compatible) model + data prep ####
# Ruppert basis for t=1..(ny-1)
build_ruppert_basis <- function(t, num_knots = 12, centre = TRUE, rescale01 = TRUE) {
  stopifnot(is.numeric(t), length(t) >= 3)
  t <- as.numeric(t)
  x <- if (rescale01) (t - min(t)) / (max(t) - min(t)) else t
  x <- if (centre) x - mean(x) else x
  knots <- stats::quantile(x, probs = seq(0, 1, length.out = num_knots + 2))[2:(num_knots + 1)]
  Z <- outer(x, as.numeric(knots), function(xx, kk) pmax(xx - kk, 0)^3)
  X <- cbind(Intercept = 1, Linear = x)
  list(X = unname(X), Z = unname(Z), knots = as.numeric(knots), x = x)
}

# For empirical matrices (or simulated y/se re-used)
prepare_jags_data_freeman <- function(mat_index, mat_se = NULL, years = NULL,
                                      num_knots = 12,
                                      obs_var_model = NULL, # preferred
                                      seFromData   = NULL, # legacy
                                      Y1perfect    = TRUE,
                                      tiny_tau     = 1e-12,
                                      m.scale      = c("loge","log10","logit")) { # Note that loge assumes that data already on natural log scale
  m.scale <- match.arg(m.scale)
  stopifnot(is.matrix(mat_index))
  nsp <- nrow(mat_index); ny <- ncol(mat_index)
  if (is.null(years)) years <- seq_len(ny)
  
  # Resolve observation-variance option
  ovm <- if (!is.null(obs_var_model)) {
    as.integer(obs_var_model)
  } else if (!is.null(seFromData)) {
    if (isTRUE(seFromData)) 2L else 4L
  } else 4L
  if (!ovm %in% 1:4) stop("obs_var_model must be 1,2,3 or 4.")
  
  # Transform index to loge and build observation mask
  Y <- .to_loge(mat_index, m.scale) # Note that loge assumes that data already on natural log scale
  obs_mask <- 1L * is.finite(Y)
  
  # Transform SEs to loge scale if provided; allow missing SEs always
  if (!is.null(mat_se)) {
    stopifnot(all(dim(mat_se) == dim(mat_index)))
    SE <- if (m.scale == "loge") {
      mat_se
    } else if (m.scale == "log10") {
      mat_se * log(10)
    } else { # logit -> loge delta-method style scaling
      pmax(1e-12, mat_se / pmax(1e-6, mat_index) / pmax(1e-6, 1 - mat_index))
    }
  } else {
    SE <- matrix(NA_real_, nsp, ny)
  }
  
  # SE-present flag only where Y is observed
  se_present <- 1L * (is.finite(SE) & is.finite(Y))
  
  # Per-cell precision from SE where available; neutral filler elsewhere
  tau_obs <- matrix(tiny_tau, nsp, ny)
  if (any(se_present == 1L)) {
    tau_obs[se_present == 1L] <- 1 / pmax(1e-12, SE[se_present == 1L])^2
  }
  
  # First observed year per species (fallback to 1 if none observed)
  FY <- apply(Y, 1, function(z) { ii <- which(is.finite(z)); if (length(ii) == 0) 1L else min(ii) })
  FY <- as.integer(FY)
  
  # Growth-year basis (Ruppert) on t = 1..(ny-1)
  basis <- build_ruppert_basis(t = 1:(ny - 1), num_knots = max(1, num_knots))
  X <- basis$X; Z <- basis$Z
  
  # Year-1-perfect mask (only at genuinely observed FY cell)
  y1perfect_mask <- matrix(0L, nsp, ny)
  has_FY_obs <- vapply(seq_len(nsp), function(s) is.finite(Y[s, FY[s]]), logical(1))
  y1perfect_mask[cbind(seq_len(nsp), FY)] <- as.integer(Y1perfect) * as.integer(has_FY_obs)
  
  list(
    nsp = nsp, nyears = ny, num_knots = ncol(Z),
    estimate = Y, X = X, Z = Z, FY = FY, years = years,
    # observation components
    obs_var_model = ovm,
    obs_mask = obs_mask,
    tau_obs  = tau_obs,
    se_present = se_present,
    y1perfect_mask = y1perfect_mask,
    big_tau  = 1e12,
    tiny_tau = tiny_tau,
    eps = 1.0e-12,
    # misc
    m_scale = m.scale,
    ones = rep(1, nsp)
  )
}


freeman_brc_model <- ("
model {

  # Fixed effects for linear + spline
  for (j in 1:2) { beta[j] ~ dnorm(0, 1.0E-6) }
  taub ~ dgamma(1.0E-6, 1.0E-6)
  for (k in 1:num_knots) { b[k] ~ dnorm(0, taub) }

  # Species growth variance
  sigma_spi ~ dunif(0.0001, 30)
  tau_spi <- pow(sigma_spi, -2)

  # Observation-variance switch
  w1 <- equals(obs_var_model, 1)  # perfect (debug)
  w2 <- equals(obs_var_model, 2)  # use SEs as given
  w3 <- equals(obs_var_model, 3)  # global scale 'c' on SEs
  w4 <- equals(obs_var_model, 4)  # common theta
  c ~ dunif(0, 10)
  theta ~ dunif(0.0001, 10)
  tau_theta <- pow(theta, -2)

  # Common smoothed growth for t = 1..(nyears-1)
  for (t in 1:(nyears-1)) {
    mfe[t] <- beta[1]*X[t,1] + beta[2]*X[t,2]
    for (k in 1:num_knots) { tmp[t,k] <- b[k] * Z[t,k] }
    mre[t] <- sum(tmp[t,1:num_knots])
    Lambda[t] <- mfe[t] + mre[t]
    Gexp[t] <- exp(Lambda[t]) # multiplicative growth
  }

  # Cumulative indicator M
  M[1] <- 0
  for (t in 2:nyears) { M[t] <- M[t-1] + Lambda[t-1] }

  # Species growth and index reconstruction
  for (s in 1:nsp) {
    for (t in 1:(nyears-1)) {
      g[s,t] ~ dnorm(Lambda[t], tau_spi)
    }
    # prefix sums of growth
    csum_g[s,1] <- 0
    for (t in 2:nyears) { csum_g[s,t] <- csum_g[s,t-1] + g[s,t-1] }

    # Latent anchor (vague prior on log scale)
    anchor[s] ~ dnorm(0, 1.0E-6)

    # reconstruct indices for all years
    for (t in 1:nyears) {
      fwd[s,t]  <- csum_g[s,t] - csum_g[s, FY[s]]
      back[s,t] <- csum_g[s, FY[s]] - csum_g[s,t]
      spindex[s,t] <- anchor[s] + step(t - FY[s]) * fwd[s,t] - step(FY[s] - t) * back[s,t]
    }
  }

  # Observation model (2-D masked, SE switch)
  for (s in 1:nsp) {
    for (t in 1:nyears) {
      # Components (only one w* is 1)
      tau_w1[s,t] <- w1 * ( obs_mask[s,t] *  big_tau + (1 - obs_mask[s,t]) * tiny_tau )
      # Allow some SEs to be missing
      tau_w2[s,t] <- w2 * ( obs_mask[s,t]*(se_present[s,t] * tau_obs[s,t] + (1 - se_present[s,t]) * tau_theta ) + (1 - obs_mask[s,t]) * tiny_tau )
      # Allow some SEs to be missing
      tau_w3[s,t] <- w3 * ( obs_mask[s,t] * ( se_present[s,t] * (tau_obs[s,t] / (pow(c,2) + eps)) + (1 - se_present[s,t]) * tau_theta ) + (1 - obs_mask[s,t]) * tiny_tau )
      tau_w4[s,t] <- w4 * ( obs_mask[s,t] *  tau_theta + (1 - obs_mask[s,t]) * tiny_tau )

      tau_core[s,t] <- tau_w1[s,t] + tau_w2[s,t] + tau_w3[s,t] + tau_w4[s,t]

      # Apply Year-1-perfect only at FY cell if requested
      tau_eff[s,t]  <- y1perfect_mask[s,t] * big_tau
                     + (1 - y1perfect_mask[s,t]) * tau_core[s,t]

      estimate[s,t] ~ dnorm(spindex[s,t], tau_eff[s,t])
    }
  }

  # Mean log index (per year) as a single expression
  for (t in 1:nyears) {
    Mprime[t] <- inprod(ones[1:nsp], spindex[1:nsp, t]) / nsp
  }
  
}
")

run_freeman_brc <- function(data_list, n_iter = 4000, n_burnin = 1000, n_chains = 3, n_thin = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  R2jags::jags(
    data = data_list,
    parameters.to.save = c("beta","b","taub","sigma_spi","theta","Lambda","Gexp",
                           "M","Mprime","spindex","g"),
    model.file = textConnection(freeman_brc_model),
    n.chains = n_chains, n.iter = n_iter, n.burnin = n_burnin, n.thin = n_thin
  )
}

# 6) Fit wrappers and convergence checks ####
fit_partial_or_nopool <- function(model = c("partial","nopool"), data_list,
                                  n_iter = 4000, n_burn = 1000, n_chains = 3, n_thin = 5) {
  model <- match.arg(model)
  model_string <- if (model=="partial") jags_partial else jags_nopool
  params <- if (model=="partial") c("l","delta","G_smooth","G_full","M_smooth","M_full","MSI","mu","mu_full",
                                    "sigma_b","sigma_sp","sigma_delta","b","beta0","c","theta")
  else c("l","G_smooth","G_full","M_smooth","M_full","MSI","mean_mu","mean_mu_full",
         "sigma_b","b","beta0","c","theta","sigma_sp")
  jags(data = data_list, parameters.to.save = params, model.file = textConnection(model_string),
       n.chains = n_chains, n.iter = n_iter, n.burnin = n_burn, n.thin = n_thin)
}

fit_bayes_geomean <- function(data_list, n_iter = 4000, n_burn = 1000, n_chains = 3, n_thin = 5,
                              impute_all = TRUE, nu_rw = NULL) {
  model_string <- if (isTRUE(impute_all)) jags_bayes_geomean_impute else jags_bayes_geomean
  if (!is.null(nu_rw)) data_list$nu_rw <- nu_rw
  jags(data = data_list,
       parameters.to.save = c("MSI","l","theta","c", "Lbar", "mean_logI"), # keep l for diagnostics
       model.file = textConnection(model_string),
       n.chains = n_chains, n.iter = n_iter, n.burnin = n_burn, n.thin = n_thin)
}

rhat_report <- function(fit, threshold = 1.01, top_n = 30) {
  sm <- as.data.frame(fit$BUGSoutput$summary)
  sm$param <- rownames(sm)
  bad <- sm %>% filter(!is.na(Rhat) & Rhat > threshold) %>% arrange(desc(Rhat)) %>% slice_head(n = top_n)
  list(n_bad = nrow(bad), table = bad)
}

# 7) Posterior summaries ####
posterior_summary <- function(arr, probs = c(0.025, 0.5, 0.975)) {
  q <- apply(arr, seq_along(dim(arr))[-1], quantile, probs = probs, na.rm = TRUE)
  list(lower = q[1,, drop=FALSE], median = q[2,, drop=FALSE], upper = q[3,, drop=FALSE])
}
extract_summaries_partial <- function(fit) {
  s <- fit$BUGSoutput$sims.list
  list(MSI = posterior_summary(s$MSI),
       log_smooth = posterior_summary(s$mu), log_full = posterior_summary(s$mu_full),
       G_smooth = posterior_summary(s$G_smooth), G_full = posterior_summary(s$G_full),
       M_smooth = posterior_summary(s$M_smooth), M_full = posterior_summary(s$M_full))
}
extract_summaries_nopool <- function(fit) {
  s <- fit$BUGSoutput$sims.list
  list(MSI = posterior_summary(s$MSI),
       log_smooth = posterior_summary(s$mean_mu), log_full = posterior_summary(s$mean_mu_full), # means across species 
       G_smooth = posterior_summary(s$G_smooth), G_full = posterior_summary(s$G_full),
       M_smooth = posterior_summary(s$M_smooth), M_full = posterior_summary(s$M_full))
}
extract_summaries_bayes_geomean <- function(fit) {
  s <- fit$BUGSoutput$sims.list
  list(MSI = posterior_summary(s$MSI))
}
extract_summaries_freeman <- function(fit) {
  s <- fit$BUGSoutput$sims.list
  list(M = posterior_summary(s$M), Mprime = posterior_summary(s$Mprime), 
       Lambda = posterior_summary(s$Lambda), Gexp = posterior_summary(s$Gexp))
}

# 8) Plotting ####
plot_indicator <- function(years, ms_list, truths = NULL, title = "MSI") {

  stopifnot(is.numeric(years))
  keep <- !vapply(ms_list, is.null, logical(1))
  ms_list <- ms_list[keep]
  if (length(ms_list) == 0) stop("No model summaries to plot.")
  
  labels <- names(ms_list)
  if (is.null(labels) || any(is.na(labels) | labels == "")) {
    labels <- paste0("Model", seq_along(ms_list))
    names(ms_list) <- labels
  }
  
  assembled <- Map(function(s, nm) {
    med <- as.numeric(s$median); lwr <- as.numeric(s$lower); upr <- as.numeric(s$upper)
    if (!length(med) || length(med) != length(years))
      stop("Length mismatch for model '", nm, "': ", length(med), " vs years ", length(years))
    tibble::tibble(model = nm, year = years, med = med, lwr = lwr, upr = upr)
  }, ms_list, labels)
  
  df <- dplyr::bind_rows(assembled)
  df$model <- as.factor(df$model)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = med, colour = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(y = "Index (baseline = 1)", x = "Year", title = title) +
    .msitheme()
  
  if (!is.null(truths)) {
    # Expect data.frame with columns year, truth
    truths <- as.data.frame(truths)
    stopifnot(all(c("year","truth") %in% names(truths)))
    p <- p + ggplot2::geom_line(data = truths, ggplot2::aes(x = year, y = truth, linetype = truth_type),
                                inherit.aes = FALSE, na.rm = TRUE)
  }
  p
}

plot_growth <- function(years, growth_list, truth = NULL, title) {
  stopifnot(is.numeric(years))
  # Keep only non-NULL entries
  keep <- !vapply(growth_list, is.null, logical(1))
  growth_list <- growth_list[keep]
  if (length(growth_list) == 0L) stop("No growth summaries to plot.")
  # Ensure labels
  labels <- names(growth_list)
  if (is.null(labels) || anyNA(labels) || any(labels == "")) {
    labels <- paste0("Model", seq_along(growth_list))
    names(growth_list) <- labels
  }
  # Assemble a tidy data frame
  assembled <- Map(function(s, nm) {
    med <- as.numeric(s$median)
    lwr <- if (!is.null(s$lower))  as.numeric(s$lower)  else rep(NA_real_, length(med))
    upr <- if (!is.null(s$upper))  as.numeric(s$upper)  else rep(NA_real_, length(med))
    if (!length(med) || length(med) != length(years)) {
      stop("Length mismatch for model '", nm, "': median length ",
           length(med), " vs years length ", length(years))
    }
    tibble::tibble(model = nm, year = years, med = med, lwr = lwr, upr = upr)
  }, growth_list, labels)
  
  df <- dplyr::bind_rows(assembled)
  df$model <- as.factor(df$model)
  
  # Core plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = med, colour = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr),
                         alpha = 0.15, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::labs(x = "Year", y = title, title = title)
  
  # Optional truth overlay: data.frame with columns year, truth
  if (!is.null(truth)) {
    truth <- as.data.frame(truth)
    if (all(c("year", "truth", "truth_type") %in% names(truth))) {
      p <- p + ggplot2::geom_line(data = truth,
                                  ggplot2::aes(x = year, y = truth, linetype = truth_type),
                                  inherit.aes = FALSE, na.rm = TRUE)
    } else {
      warning("`truth` must have columns 'year' and 'truth'; layer omitted.")
    }
  }
  
  # Theme (use project theme if present)
  if (exists(".msitheme", mode = "function")) {
    p <- p + .msitheme()
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  p
}

# Cumulative common log-growth (length T; x = years)
# cum_list: named list of posterior summaries, each a list(lower, median, upper),
# where each vector has length(years)
# truths: optional data.frame with columns year, truth  (or NULL)
plot_cumlog <- function(years, cum_list, truths = NULL,
                        title = "Cumulative common log-growth (M)",
                        ylab  = "Cumulative log-growth") {
  stopifnot(is.numeric(years))
  if (!length(cum_list)) stop("No cumulative series supplied.")
  
  # Keep only non-NULL entries and ensure labels
  keep <- !vapply(cum_list, is.null, logical(1))
  cum_list <- cum_list[keep]
  if (!length(cum_list)) stop("All cumulative series were NULL.")
  labels <- names(cum_list)
  if (is.null(labels) || anyNA(labels) || any(labels == "")) {
    labels <- paste0("Model", seq_along(cum_list))
    names(cum_list) <- labels
  }
  
  # Assemble tidy data
  assembled <- Map(function(s, nm) {
    med <- as.numeric(s$median)
    lwr <- if (!is.null(s$lower))  as.numeric(s$lower)  else rep(NA_real_, length(med))
    upr <- if (!is.null(s$upper))  as.numeric(s$upper)  else rep(NA_real_, length(med))
    if (!length(med) || length(med) != length(years)) {
      stop("Length mismatch for model '", nm, "': median length ",
           length(med), " vs years length ", length(years))
    }
    tibble::tibble(model = nm, year = years, med = med, lwr = lwr, upr = upr)
  }, cum_list, labels)
  df <- dplyr::bind_rows(assembled)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = med, colour = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr),
                         alpha = 0.15, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::labs(x = "Year", y = ylab, title = title)
  
  # Optional truth overlay
  if (!is.null(truths)) {
    truths <- as.data.frame(truths)
    if (all(c("year", "truth", "truth_type") %in% names(truths))) {
      p <- p + ggplot2::geom_line(
        data = truths,
        ggplot2::aes(x = year, y = truth, linetype = truth_type),
        inherit.aes = FALSE, na.rm = TRUE
      )
    } else {
      warning("`truths` must have columns 'year' and 'truth'; layer omitted.")
    }
  }
  
  if (exists(".msitheme", mode = "function")) p <- p + .msitheme() else p <- p + ggplot2::theme_minimal()
  p
}

# 9) Empirical input helpers (matrices) ####
# Build a sim-like list from empirical matrices
as_sim_from_empirical <- function(index_mat, se_mat = NULL, years = NULL,
                                  m.scale = c("loge","log10","logit"),
                                  eps_p = 1e-6) {
  m.scale <- match.arg(m.scale)
  stopifnot(is.matrix(index_mat))
  if (!is.null(se_mat)) {
    stopifnot(is.matrix(se_mat), all(dim(se_mat) == dim(index_mat)))
  }
  nsp <- nrow(index_mat); ny <- ncol(index_mat)
  
  # transform index to loge
  idx_loge <- .to_loge(index_mat, m.scale)
  
  # transform SE to loge (allow partial missing)
  if (is.null(se_mat)) {
    se_loge <- matrix(NA_real_, nsp, ny)
  } else if (m.scale == "loge") {
    se_loge <- se_mat
  } else if (m.scale == "log10") {
    se_loge <- se_mat * log(10)
  } else { # "logit": index on probability scale, SE on p-scale
    # clamp p away from 0/1 for numerical stability; preserve NA where index missing
    p <- index_mat
    p_clamped <- pmin(pmax(p, eps_p), 1 - eps_p)
    se_loge <- se_mat / (p_clamped * (1 - p_clamped))
    se_loge[!is.finite(p)] <- NA_real_
  }
  
  # ensure SEs are NA wherever the index is missing (all scales)
  se_loge[!is.finite(idx_loge)] <- NA_real_
  
  list(
    years = years %||% seq_len(ny),
    y     = idx_loge,
    se    = se_loge,
    l_true = NULL, mu_true = NULL,
    delta_true = rep(0, ny - 1),
    gamma_s = rep(0, nsp),
    innov_dist = NA, df_u = NA,
    inclusion_bias = list(enabled = FALSE),
    I = matrix(1L, nsp, ny),
    use_delta = FALSE,
    meta = list(empirical = TRUE)
  )
}

# 10) Full workflow runner ####
run_full_analysis <- function(data_source = c("simulate","empirical"),
                              sim_args = list(),
                              # if there is actual data and data_source = "empirical"
                              empirical_index_mat = NULL, empirical_se_mat = NULL, empirical_years = NULL, empirical_m_scale = "loge",
                              fit_models = c("partial","freeman","nopool","bayes_geomean"),
                              # fitting controls
                              jags_partial = list(df_mu = 6, obs_var_model = 3, n_iter = 4000, n_burn = 1000, n_chains = 3, n_thin = 5),
                              jags_nopool  = list(df_mu = 6, obs_var_model = 3, n_iter = 4000, n_burn = 1000, n_chains = 3, n_thin = 5),
                              jags_freeman = list(obs_var_model = 3, n_iter = 4000, n_burnin = 1000, n_chains = 3, n_thin = 5, seed = NULL),
                              jags_bayes_geomean = list(obs_var_model = 3, n_iter = 4000, n_burn = 1000, n_chains = 3, n_thin = 5),
                              smooth_geomean = list(enable = TRUE, prefer_freeman_basis = TRUE),
                              plot_geomean = TRUE, impute_all_geomean = TRUE,
                              plot_MNAR = TRUE,
                              brc_opts = list(num_knots = 12, seFromData = FALSE, Y1perfect = TRUE, m.scale = "loge"),
                              growth_scale = c("log", "exp")) {
  data_source <- match.arg(data_source)
  if (data_source == "simulate") {
    sim <- do.call(simulate_species_data, sim_args)
  } else {
    stopifnot(!is.null(empirical_index_mat))
    if (data_source == "empirical" && brc_opts$m.scale != "loge") {
      stop("For empirical data, brc_opts$m.scale must be 'loge' because indices are ",
           "already converted to loge in as_sim_from_empirical().")
    }
    sim <- as_sim_from_empirical(empirical_index_mat, empirical_se_mat, empirical_years, m.scale = empirical_m_scale)
  }
  
  results <- list(); fits <- list(); checks <- list()
  
  # Partial
  if ("partial" %in% fit_models) {
    dl <- prepare_jags_data2(sim, df_mu = jags_partial$df_mu %||% 6, obs_var_model = jags_partial$obs_var_model %||% 3)
    fit <- fit_partial_or_nopool("partial", dl, n_iter = jags_partial$n_iter %||% 4000, n_burn = jags_partial$n_burn %||% 1000,
                                 n_chains = jags_partial$n_chains %||% 3, n_thin = jags_partial$n_thin %||% 2)
    fits$partial <- fit
    checks$partial <- rhat_report(fit)
    results$partial <- extract_summaries_partial(fit)
  }
  
  # No pool
  if ("nopool" %in% fit_models) {
    dl <- prepare_jags_data2(sim, df_mu = jags_nopool$df_mu %||% 6, obs_var_model = jags_nopool$obs_var_model %||% 3)
    fit <- fit_partial_or_nopool("nopool", dl, n_iter = jags_nopool$n_iter %||% 4000, n_burn = jags_nopool$n_burn %||% 1000,
                                 n_chains = jags_nopool$n_chains %||% 3, n_thin = jags_nopool$n_thin %||% 5)
    fits$nopool <- fit
    checks$nopool <- rhat_report(fit)
    results$nopool <- extract_summaries_nopool(fit)
  }
  
  # Freeman (BRC)
  if ("freeman" %in% fit_models) {
    # choose option: prefer explicit brc_opts$obs_var_model, else map seFromData -> {2,4}
    freeman_ovm <- jags_freeman$obs_var_model %||% if (isTRUE(brc_opts$seFromData)) 2L else 4L
    
    dat <- prepare_jags_data_freeman(sim$y, sim$se,
                                     years = sim$years,
                                     num_knots = brc_opts$num_knots %||% 12,
                                     obs_var_model = freeman_ovm,
                                     Y1perfect = brc_opts$Y1perfect %||% TRUE,
                                     m.scale = brc_opts$m.scale %||% "loge")
    fit <- run_freeman_brc(dat, n_iter = jags_freeman$n_iter %||% 4000, n_burnin = jags_freeman$n_burnin %||% 1000,
                           n_chains = jags_freeman$n_chains %||% 3, n_thin = jags_freeman$n_thin %||% 5,
                           seed = jags_freeman$seed %||% NULL)
    fits$freeman <- fit
    checks$freeman <- rhat_report(fit)
    results$freeman <- extract_summaries_freeman(fit)
  }
  
  # Bayesian geomean (no pooling across species beyond averaging; no time model)
  if ("bayes_geomean" %in% fit_models) {
    dl <- prepare_jags_data2(sim, df_mu = 6, obs_var_model = jags_bayes_geomean$obs_var_model %||% 4)
    fit <- fit_bayes_geomean(dl, n_iter = jags_bayes_geomean$n_iter %||% 4000, n_burn = jags_bayes_geomean$n_burn %||% 1000,
                             n_chains = jags_bayes_geomean$n_chains %||% 3, n_thin = jags_bayes_geomean$n_thin %||% 5,
                             impute_all = impute_all_geomean, nu_rw = NULL)
    fits$bayes_geomean <- fit
    checks$bayes_geomean <- rhat_report(fit)
    results$bayes_geomean <- extract_summaries_bayes_geomean(fit)
    # smooth Bayesian geomean
    if (isTRUE(smooth_geomean$enable)) {
      want_ruppert <- isTRUE(smooth_geomean$prefer_freeman_basis) && ("freeman" %in% fit_models)
      basis_mode   <- if (want_ruppert) "ruppert" else "ns"
      sm <- posthoc_smooth_bayes_geomean(
        fit = fit,
        years = sim$years,
        impute_all = impute_all_geomean,
        basis = basis_mode,
        df_mu = jags_partial$df_mu %||% 6,
        num_knots = brc_opts$num_knots %||% 12
      )
      results$bayes_geomean_smooth <- sm
    }
  }
  
  ## Plots ####
  growth_scale <- match.arg(growth_scale)
  plots <- list()
  
  ## MSI ####
  MSI_list <- list()
  if (!is.null(results$partial)) MSI_list$Partial <- results$partial$MSI
  if (!is.null(results$nopool))  MSI_list$Nopool  <- results$nopool$MSI
  if (!is.null(results$bayes_geomean) && isTRUE(plot_geomean)) MSI_list$BayesGeomean <- results$bayes_geomean$MSI
  if (!is.null(results$bayes_geomean_smooth)) MSI_list$BayesGeomean_S <- results$bayes_geomean_smooth$MSI
  if (!is.null(results$freeman)) {
    if (!is.null(fits$freeman$BUGSoutput$sims.list$Mprime)) {
      Mdraws <- fits$freeman$BUGSoutput$sims.list$Mprime # iters x T
      MSI_draws <- exp(Mdraws - Mdraws[, 1])
      MSI_list$Freeman <- posterior_summary(MSI_draws)
    } else {
      Mpr <- results$freeman$Mprime
      MSI_list$Freeman <- list(
        lower  = exp(Mpr$lower - Mpr$lower[1]),
        median = exp(Mpr$median - Mpr$median[1]),
        upper  = exp(Mpr$upper - Mpr$upper[1])
      )
    }
  }
  truth_idx <- NULL
  if (!is.null(sim$l_true)) {
    Lbar <- colMeans(sim$l_true)
    truth_idx <- tibble::tibble(year = sim$years, truth = exp(Lbar - Lbar[1]), truth_type = "All spp.")
    if (!is.null(sim$lmnar_selected_true) && isTRUE(plot_MNAR)) {
      Lbar_mnar <- colMeans(sim$lmnar_selected_true, na.rm = TRUE)
      truth_idx <- tibble::add_row(truth_idx, year = sim$years, truth = exp(Lbar_mnar - Lbar_mnar[1]), truth_type = "Sampled spp.")
    }
  }
  if (length(MSI_list) > 0) {
    plots$MSI <- plot_indicator(sim$years, MSI_list, truths = truth_idx,
                                title = "Multi-species indicator (M_prime)")
  }
  
  ## Growth (interval: length T-1; x = years[-1]) ####
  years_m1 <- sim$years[-1]
  growth_list <- list()
  
  grab <- function(x, ...) { nms <- c(...); for (nm in nms) if (!is.null(x[[nm]])) return(x[[nm]]); NULL }
  
  truth_g <- NULL
  if (growth_scale == "log") {
    if (!is.null(results$partial)) growth_list$Partial <- grab(results$partial, "log_smooth")
    if (!is.null(results$nopool))  growth_list$Nopool  <- grab(results$nopool,  "log_smooth")
    if (!is.null(results$freeman)) growth_list$Freeman <- results$freeman$Lambda
    if (!is.null(results$bayes_geomean_smooth)) growth_list$BayesGeomean_S <- results$bayes_geomean_smooth$Lambda
    ylab_g <- "Common log-growth"
    if (!is.null(sim$mu_common_true)) {
      truth_g <- tibble::tibble(year = years_m1, truth = sim$mu_common_true, truth_type = "All spp.")
      if (!is.null(sim$mu_selected_true) && isTRUE(plot_MNAR)) {
        truth_g <- tibble::add_row(truth_g, year = years_m1, truth = sim$mu_selected_true, truth_type = "Sampled spp.")
      }
    }
  } else {
    if (!is.null(results$partial)) growth_list$Partial <- grab(results$partial, "G_smooth", "G")
    if (!is.null(results$nopool))  growth_list$Nopool  <- grab(results$nopool,  "G_smooth", "G")
    if (!is.null(results$bayes_geomean_smooth)) growth_list$BayesGeomean_S <- results$bayes_geomean_smooth$Gexp
    if (!is.null(results$freeman)) growth_list$Freeman <- results$freeman$Gexp
    ylab_g <- "Common multiplicative growth"
    if (!is.null(sim$mu_common_true)) {
      truth_g <- tibble::tibble(year = years_m1, truth = exp(sim$mu_common_true), truth_type = "All spp.")
      if (!is.null(sim$mu_selected_true) && isTRUE(plot_MNAR)) {
        truth_g <- tibble::add_row(truth_g, year = years_m1, truth = exp(sim$mu_selected_true), truth_type = "Sampled spp.")
      }
    }
  }
  
  if (length(growth_list) > 0) {
    ln <- vapply(growth_list, function(s) length(s$median), integer(1))
    stopifnot(all(ln == length(years_m1)))
    plots$Growth <- plot_growth(years_m1, growth_list, truth = truth_g, title = ylab_g)
  }
  
  ## Cumulative M (length T; x = years) ####
  cum_list <- list()
  if (!is.null(results$partial) && !is.null(results$partial$M_smooth)) cum_list$Partial <- results$partial$M_smooth
  if (!is.null(results$nopool)  && !is.null(results$nopool$M_smooth))  cum_list$Nopool  <- results$nopool$M_smooth
  if (!is.null(results$freeman) && !is.null(results$freeman$M)) cum_list$Freeman <- results$freeman$M
  if (!is.null(results$bayes_geomean_smooth)) cum_list$BayesGeomean_S <- results$bayes_geomean_smooth$M
  
  truth_M <- NULL
  if (!is.null(sim$M_common_true)) {
    truth_M <- tibble::tibble(year = sim$years, truth = sim$M_common_true, truth_type = "All spp.")
    if (!is.null(sim$M_selected_true) && isTRUE(plot_MNAR)) {
    truth_M <- tibble::add_row(truth_M, year = sim$years, truth = sim$M_selected_true, truth_type = "Sampled spp.")
    }
  }
  
  if (length(cum_list) > 0) {
    lnM <- vapply(cum_list, function(s) length(s$median), integer(1))
    stopifnot(all(lnM == length(sim$years)))
    plots$CumLog <- plot_cumlog(sim$years, cum_list, truths = truth_M,
                                title = "Cumulative common log-growth (M)")
  }
  
  list(sim = sim, fits = fits, checks = checks, results = results, plots = plots)
}

# 10) Inclusion diagnostics repeated (used externally too) ####
transition_probs <- function(I) {
  Ilag <- I[, -ncol(I), drop=FALSE]
  Icur <- I[, -1, drop=FALSE]
  N11 <- sum(Ilag == 1 & Icur == 1)
  N10 <- sum(Ilag == 1 & Icur == 0)
  N01 <- sum(Ilag == 0 & Icur == 1)
  N00 <- sum(Ilag == 0 & Icur == 0)
  p11 <- ifelse((N11 + N10) > 0, N11 / (N11 + N10), NA_real_)
  p00 <- ifelse((N01 + N00) > 0, N00 / (N01 + N00), NA_real_)
  pi  <- mean(Icur)
  list(p11 = p11, p00 = p00, pi = pi,
       Ein = ifelse(!is.na(p11), 1/(1 - p11), NA_real_),
       Eout= ifelse(!is.na(p00), 1/(1 - p00), NA_real_))
}

run_lengths <- function(I, v = 1L) {
  rl <- integer()
  for (s in 1:nrow(I)) {
    x <- as.integer(I[s,] == v)
    if (all(x == 0)) next
    r <- rle(x)
    rl <- c(rl, r$lengths[r$values == 1])
  }
  rl
}

burstiness_index <- function(x) {
  m <- mean(x); s <- sd(x)
  if (!is.finite(m) || !is.finite(s) || (m + s == 0)) return(NA_real_)
  (s - m) / (s + m)
}

evaluate_inclusion_process <- function(sim, verbose = TRUE) {
  stopifnot(!is.null(sim$I))
  I <- sim$I
  trans <- transition_probs(I)
  rl1 <- run_lengths(I, 1L)
  rl0 <- run_lengths(I, 0L)
  B1 <- burstiness_index(rl1)
  B0 <- burstiness_index(rl0)
  nY <- length(sim$years)
  
  # r_{s,t-1} used in selection
  if (!is.null(sim$mu_true)) {
    r_mat <- outer(rep(1, nrow(I)), sim$mu_true + sim$delta_true) + outer(sim$gamma_s, rep(1, nY-1))
    r_src <- "mu_true + delta_true (+ gamma_s)"
  } else if (!is.null(sim$mu_list)) {
    mu_vec_t <- function(t) sapply(seq_len(nrow(I)), function(s) sim$mu_list[[ sim$guild[s] ]][t])
    r_mat <- t(sapply(1:(nY-1), function(t) mu_vec_t(t) + sim$delta_true[t] + sim$gamma_s))
    r_src <- "mu_list (+ delta_true, gamma_s)"
  } else {
    r_mat <- matrix(NA_real_, nrow = nY-1, ncol = nrow(I))
    r_src <- "unavailable (no mu_true/mu_list)"
  }
  
  years_idx <- 2:nY
  T1 <- nY - 1
  cor_t   <- rep(NA_real_, T1)
  n_pairs <- integer(T1)
  sd_Iv   <- rep(NA_real_, T1)
  sd_rv   <- rep(NA_real_, T1)
  reason  <- rep(NA_character_, T1)
  
  for (tt in seq_len(T1)) {
    t <- years_idx[tt]
    # orient r_mat so rows = time (T-1), cols = species (S)
    r_t <- if (nrow(r_mat) == T1 && ncol(r_mat) == nrow(I)) {
      as.numeric(r_mat[tt, ])
    } else if (ncol(r_mat) == T1 && nrow(r_mat) == nrow(I)) {
      as.numeric(r_mat[, tt])
    } else {
      stop("r_mat must be (T-1)×S or S×(T-1); got ",
           paste(dim(r_mat), collapse = "×"), " vs I dims ",
           paste(dim(I), collapse = "×"))
    }
    I_t <- as.numeric(I[, t])
    ok <- is.finite(I_t) & is.finite(r_t)
    n_ok <- sum(ok)
    n_pairs[tt] <- n_ok
    if (n_ok < 3L) { reason[tt] <- "too_few_pairs"; next }
    
    sdI <- stats::sd(I_t[ok]); sdR <- stats::sd(r_t[ok])
    sd_Iv[tt] <- sdI; sd_rv[tt] <- sdR
    
    if (!is.finite(sdI) || sdI == 0) {
      reason[tt] <- if (!is.finite(sdR) || sdR == 0) "no_variation_both" else "no_variation_I"
      next
    }
    if (!is.finite(sdR) || sdR == 0) { reason[tt] <- "no_variation_r"; next }
    
    cor_t[tt] <- stats::cor(I_t[ok], r_t[ok])
    reason[tt] <- "ok"
  }
  
  if (all(!is.finite(r_mat))) {
    reason[] <- "no_generator_growth"
  }
  
  diag_df <- data.frame(
    year   = years_idx,
    n_pairs = n_pairs,
    sd_I    = sd_Iv,
    sd_r    = sd_rv,
    reason  = reason,
    stringsAsFactors = FALSE
  )
  
  note <- NULL
  if (all(is.na(cor_t))) {
    # summarise reasons (including counts for clarity)
    lvl <- c("too_few_pairs","no_variation_I","no_variation_r","no_variation_both","no_generator_growth","ok")
    freq <- table(factor(reason, levels = lvl))
    # keep nonzero only
    freq <- freq[freq > 0]
    top <- if (length(freq)) paste(paste0(names(freq), "×", as.integer(freq)), collapse = ", ") else "no classified reason"
    note <- paste0(
      "All annual selection–growth correlations are NA for T-1 = ", T1, ". ",
      "This is expected when inclusion or growth lacks cross-species variation (e.g., inclusion all 1s). ",
      "Reason summary: [", top, "]. Source of r: ", r_src, "."
    )
    if (isTRUE(verbose)) message(note)
  }
  
  list(
    transition = trans,
    run_lengths_in = rl1,
    run_lengths_out = rl0,
    burstiness_in = B1,
    burstiness_out = B0,
    cor_I_r_by_year = cor_t,
    diagnostics = diag_df,
    note = note
  )
}

# 11) Example usage ####
sim_args <- list(n_species = 40, n_years = 30, seed = 232680,
                 ## Growth rate DGP
                 dgp_mode = "spline", # "spline","rw","ar1","changepoint","seasonal","mixture"
                 ## DGP options
                 # spline options
                 df_mu = 6, # default
                 # RW / AR1
                 sigma_eta = 0.05, phi = 0.7,
                 # changepoint
                 n_cp = 2, jump_sd = 0.15, piecewise_linear = FALSE,
                 # seasonal
                 A = 0.15, period = 8, phi0 = runif(1, 0, 2*pi),
                 # mixture of random walks
                 K_guilds = 5, # number of groups
                 ## Simulate species options
                 # species/state variation
                 sigma_sp = 0.05, sigma_delta = 0.05, sd_alpha0 = 0.4,
                 innov_dist = "normal", df_u = 3, # "normal" or "student_t" (df_u for t)
                 sp_trend = "random_slope", sigma_gamma = 0.05, log_sd_se = 0.35, # "none" or "random_slope"
                 use_delta = FALSE, # cross-species annual shocks?
                 # observation model
                 sigma_obs_mean = 0.02, prop_missing = 0,
                 # Inclusion MNAR (additional to prop_missing)
                 inclusion_bias = list(enabled=F, # include inclusion bias?
                                       a0=-0.4, 
                                       a1=3.0,
                                       rho1=1.2,
                                       rho0=0.2,
                                       p_init=1))
out <- run_full_analysis(data_source = "simulate", # simulated data or empirical?
                         sim_args = sim_args, # as above
                         #fit_models = c("partial","freeman", "nopool","bayes_geomean"),
                         #fit_models = c("partial","freeman","bayes_geomean"),
                         #fit_models = c("partial","freeman", "nopool"),
                         #fit_models = c("partial"),
                         #fit_models = c("freeman"),
                         #fit_models = c("freeman","bayes_geomean"),
                         fit_models = c("bayes_geomean"),
                         ## Model-specific settings
                         jags_partial = list(df_mu=6, obs_var_model=2, n_iter=500, n_burn=100),
                         jags_nopool = list(df_mu=6, obs_var_model=2, n_iter=500, n_burn=100),
                         jags_freeman = list(obs_var_model=2, n_iter=500, n_burnin=100),
                         jags_bayes_geomean = list(obs_var_model=2, n_iter=500, n_burnin=100),
                         smooth_geomean = list(enable = TRUE, prefer_freeman_basis = FALSE),
                         plot_geomean = TRUE, # plot unsmoothed Bayes geomean MSI
                         # impute missing spp/year combinations
                         impute_all_geomean = T, # if false, then under MNAR the estimand differs from other models
                         plot_MNAR = TRUE, # plot sampled species "truth" alongside actual truth
                         ## Cross-model settings
                         # legacy seFromData overridden by model-specific obs_var_model settings, seFromData = T (or F) sets obs_var_model = 2 (or 4
                         # but only if obs_var_model not specified
                         # m.scale = loge assumes data already on natural log scale
                         brc_opts = list(num_knots=12, seFromData=TRUE, Y1perfect=TRUE, m.scale = "loge"),
                         ## Growth-rate presentation scale (model returns both anyway, this is for plot)
                         growth_scale = "log")
# Convergence report
out$checks
# Plots
print(out$plots$MSI) # M_prime equivalents
print(out$plots$CumLog) # M equivalents where applicable
print(out$plots$Growth) # annual growth on scale specified in run_full_analysis()
# Inclusion diagnostics
evaluate_inclusion_process(out$sim, verbose = TRUE)

#### TO DO ####
## Fix back-transformation of logit-scale empirical data (currently uses exp(), which results in different estimand)
## Shrinkage evaluation!
## Visualise (sample of) species' trends

# End of file #
