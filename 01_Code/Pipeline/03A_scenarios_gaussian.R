# ==============================================================================
# 03A_scenarios_gaussian.R
# ------------------------------------------------------------------------------
# MODEL A: Default Gaussian assumption -- no deviation from factorstochvol.
#
# Generates N_SCENARIOS one-step-ahead return scenarios from the POSTERIOR
# PREDICTIVE distribution under the model's Gaussian innovation assumption:
#
#   f_{j,T+1} | h_{m+j,T+1} ~ N(0, exp(h_{m+j,T+1}))
#   eps_{i,T+1}| h_{i,T+1}  ~ N(0, exp(h_{i,T+1}))
#
# Procedure (one scenario per post-burnin MCMC draw):
#   Step 1. For draw s, retrieve h_T^(s) from logvol(fit)[, T, s].
#   Step 2. Advance each log-variance one step via the posterior-mean AR(1):
#             h_{T+1,i}^(s) ~ N(mu_i + phi_i*(h_{T,i}^(s) - mu_i), sigma_i^2)
#           where (mu_i, phi_i, sigma_i) are the posterior-mean AR(1) parameters
#           estimated in 02_model_fit.R.
#           [DEVIATION NOTE: A fully Bayesian approach would use draw-specific
#            (mu^s, phi^s, sigma^s). Here we use posterior means, which
#            underestimates SV parameter uncertainty in the tails.]
#   Step 3. Draw factor innovations  f^(s) ~ N(0, diag(exp(h_{factor,T+1}^(s)))).
#   Step 4. Draw idio innovations  eps^(s) ~ N(0, diag(exp(h_{idio,T+1}^(s)))).
#   Step 5. Form return scenario: y^(s) = Lambda^(s) %*% f^(s) + eps^(s).
#
# Output (saved to PATH_DATA):
#   scenarios_gaussian.rds  -- N_SCENARIOS x m matrix of return scenarios
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

fsv_fit      <- readRDS(file.path(PATH_DATA, "fsv_fit.rds"))
sv_params    <- readRDS(file.path(PATH_DATA, "sv_ar1_params.rds"))  # [m+r, 3]
returns_dem  <- readRDS(file.path(PATH_DATA, "returns_demeaned.rds"))

m     <- length(TICKERS)
r     <- N_FACTORS
T_obs <- nrow(returns_dem)
S     <- N_SCENARIOS       # = MCMC_DRAWS - MCMC_BURNIN = 40 000

cat(sprintf("Model A (Gaussian): generating %d scenarios for m=%d assets.\n", S, m))

# ── Extract Posterior Draws ────────────────────────────────────────────────────
# logvol(fit)[i, t, s] = log-VOLATILITY of series i at time t under draw s.
# Volatility  = exp(logvol_val)
# Log-variance (paper's h) = 2 * logvol_val
lv_draws     <- logvol(fsv_fit)   # [m+r, T, S]
lambda_draws <- loadings(fsv_fit) # [m, r, S]

# Log-variance at the LAST observation time T for each draw
# Shape: [m+r, S]  (log-variance at T for all series and all draws)
logvar_T <- 2 * lv_draws[, T_obs, ]   # 2 * log-vol = log-var

# ── AR(1) One-Step-Ahead Advance ──────────────────────────────────────────────
# For each draw s, advance h_{T+1, i}^(s) using posterior-mean AR(1) parameters.
# h_{T+1} = mu + phi * (h_T^(s) - mu) + sigma * N(0,1)
#
# Vectorised form: logvar_T1[i, s] = mu[i] + phi[i]*(logvar_T[i,s] - mu[i]) + sigma[i]*z
mu_vec    <- sv_params[, "mu"]     # [m+r] posterior-mean log-variance levels
phi_vec   <- sv_params[, "phi"]    # [m+r] persistence parameters
sigma_vec <- sv_params[, "sigma"]  # [m+r] vol-of-vol parameters

set.seed(42)
# Draw standard normals for AR(1) innovation: [m+r, S]
z_sv <- matrix(rnorm((m + r) * S), nrow = m + r, ncol = S)

# Advance: logvar_T1[i,s] = mu_i + phi_i*(logvar_T[i,s] - mu_i) + sigma_i * z[i,s]
logvar_T1 <- mu_vec + phi_vec * (logvar_T - mu_vec) + sigma_vec * z_sv  # [m+r, S]

# Convert forecast log-variance back to volatility (conditional std dev)
vol_T1 <- exp(0.5 * logvar_T1)   # [m+r, S]  vol = exp(h/2) = exp(logvar/2)

vol_idio   <- vol_T1[1:m, ]          # [m, S]  idiosyncratic volatilities
vol_factor <- vol_T1[(m+1):(m+r), ]  # [r, S]  factor volatilities

# ── Draw Gaussian Innovations and Form Return Scenarios ───────────────────────
# For each draw s:
#   f^(s) ~ N(0, diag(vol_factor[,s]^2))  -->  f = vol_factor * z_f
#   e^(s) ~ N(0, diag(vol_idio[,s]^2))    -->  e = vol_idio   * z_e
#   y^(s) = Lambda^(s) %*% f^(s) + e^(s)

z_factor <- matrix(rnorm(r * S), nrow = r, ncol = S)   # [r, S] standard normals
z_idio   <- matrix(rnorm(m * S), nrow = m, ncol = S)   # [m, S] standard normals

# Scale by forecast volatility: element-wise multiply
f_scenarios   <- vol_factor * z_factor   # [r, S]  factor scenario draws
eps_scenarios <- vol_idio   * z_idio     # [m, S]  idio  scenario draws

# Apply factor loadings: Lambda^(s) %*% f^(s) for each s
# lambda_draws is [m, r, S]; f_scenarios is [r, S]
# sapply loops over draws, giving [m, S]
common_component <- sapply(seq_len(S), function(s) {
  lambda_draws[, , s] %*% f_scenarios[, s]
})   # [m, S]

# Return scenarios: y = common + idiosyncratic
scenarios_A <- t(common_component + eps_scenarios)   # [S, m]
colnames(scenarios_A) <- TICKERS

# ── Sanity Checks ─────────────────────────────────────────────────────────────
cat("Scenario summary (Gaussian, model A):\n")
cat("  Dimensions:", dim(scenarios_A), "\n")
cat("  Column means (should be ~0, demeaned):", round(colMeans(scenarios_A), 5), "\n")
cat("  Column sds (annualised %):\n")
print(round(apply(scenarios_A, 2, sd) * sqrt(252) * 100, 2))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(scenarios_A, file.path(PATH_DATA, "scenarios_gaussian.rds"))
cat("\nScenarios saved: scenarios_gaussian.rds\n")
cat("03A_scenarios_gaussian.R complete.\n")
