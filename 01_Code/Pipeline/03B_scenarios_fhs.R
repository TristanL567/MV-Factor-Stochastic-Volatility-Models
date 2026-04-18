# ==============================================================================
# 03B_scenarios_fhs.R
# ------------------------------------------------------------------------------
# MODEL B: Filtered Historical Simulation (FHS).
#
# DEVIATION FROM GAUSSIAN ASSUMPTION (Model A):
#   In Model A, factor and idiosyncratic innovations are drawn from N(0,1).
#   Here we REPLACE that parametric draw with a BOOTSTRAP RESAMPLE from the
#   empirical distribution of standardised residuals.
#
#   The factor SV model is used for what it does well: estimating time-varying
#   volatility structure (Lambda, h_t). The INNOVATION DISTRIBUTION is taken
#   directly from the data, making no parametric tail assumption.
#
#   This approach (Barone-Adesi, Giannopoulos & Vosper 1999) addresses the
#   documented excess kurtosis in the standardised residuals (see 02_model_fit.R)
#   and produces heavier-tailed ES estimates than Model A.
#
# Procedure:
#   Steps 1-2: identical to Model A (advance log-variance AR(1), get vol_T1).
#   Step 3: instead of N(0,1), RESAMPLE (with replacement) from the empirical
#           distribution of the standardised residuals computed in 02_model_fit.R.
#   Steps 4-5: identical to Model A (scale by forecast vol, form y^(s)).
#
# Output (saved to PATH_DATA):
#   scenarios_fhs.rds             -- N_SCENARIOS x m matrix of return scenarios
#   fhs_residual_diagnostics.rds  -- diagnostics on standardised residual iid-ness
#
# Figures (saved to PATH_CHART):
#   qq_std_residuals.png  -- Q-Q plots documenting departure from N(0,1)
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

fsv_fit      <- readRDS(file.path(PATH_DATA, "fsv_fit.rds"))
sv_params    <- readRDS(file.path(PATH_DATA, "sv_ar1_params.rds"))
std_resids   <- readRDS(file.path(PATH_DATA, "std_residuals.rds"))  # [m+r, T]
returns_dem  <- readRDS(file.path(PATH_DATA, "returns_demeaned.rds"))

m     <- length(TICKERS)
r     <- N_FACTORS
T_obs <- nrow(returns_dem)
S     <- N_SCENARIOS

cat(sprintf("Model B (FHS): generating %d scenarios for m=%d assets.\n", S, m))
cat("DEVIATION: Gaussian N(0,1) innovation draws replaced by bootstrap resampling\n")
cat("           from empirical standardised residuals.\n")

# ── 1. Document Departure from Gaussianity (Q-Q Plots) ────────────────────────
# This figure motivates why we deviate from Model A: the standardised residuals
# have heavier tails than N(0,1), visible as S-curve deviations in the Q-Q plot.
series_labels <- rownames(std_resids)   # TICKERS + Factor names

png(file.path(PATH_CHART, "qq_std_residuals.png"),
    width = 1600, height = 1200, res = 120)
par(mfrow = c(ceiling((m + r) / 3), 3), mar = c(4, 4, 3, 1))
for (i in seq_len(m + r)) {
  z <- std_resids[i, ]
  qqnorm(z, main = paste0("Q-Q: ", series_labels[i]),
         col = adjustcolor("steelblue", alpha.f = 0.4), pch = 16, cex = 0.5)
  qqline(z, col = "tomato", lwd = 1.5)   # red line = perfect Gaussian
}
dev.off()
cat("Q-Q plots saved: qq_std_residuals.png\n")

# ── 2. Log-Variance Forecasts (identical to Model A) ──────────────────────────
# Steps 1-2 of Model A are reused exactly.
lv_draws  <- logvol(fsv_fit)               # [m+r, T, S]
lambda_draws <- loadings(fsv_fit)          # [m, r, S]
logvar_T  <- 2 * lv_draws[, T_obs, ]      # [m+r, S]  log-var at time T

mu_vec    <- sv_params[, "mu"]
phi_vec   <- sv_params[, "phi"]
sigma_vec <- sv_params[, "sigma"]

set.seed(43)
z_sv      <- matrix(rnorm((m + r) * S), nrow = m + r, ncol = S)
logvar_T1 <- mu_vec + phi_vec * (logvar_T - mu_vec) + sigma_vec * z_sv  # [m+r, S]
vol_T1    <- exp(0.5 * logvar_T1)   # [m+r, S]

vol_idio   <- vol_T1[1:m, ]
vol_factor <- vol_T1[(m+1):(m+r), ]

# ── 3. Bootstrap Resample from Empirical Standardised Residuals ───────────────
# DEVIATION FROM MODEL A STARTS HERE.
# We resample from the EMPIRICAL joint distribution of standardised residuals.
# Using joint (block) resampling preserves contemporaneous cross-sectional
# dependence that exists beyond the factor structure.
#
# std_resids is [m+r, T]: each column is one observation in time.
# We sample T indices WITH REPLACEMENT to get S bootstrap scenarios.
T_hist   <- ncol(std_resids)              # historical observations available
boot_idx <- sample(seq_len(T_hist), size = S, replace = TRUE)  # [S] time indices

# Resampled standardised innovations: [m+r, S]
z_boot   <- std_resids[, boot_idx]

# Split into factor and idiosyncratic parts
z_factor_boot <- z_boot[(m+1):(m+r), ]   # [r, S]  resampled factor std innovations
z_idio_boot   <- z_boot[1:m, ]           # [m, S]  resampled idio  std innovations

# ── 4. Scale by Forecast Volatility (identical to Model A) ───────────────────
# Element-wise scaling: innovation * model-forecast vol gives the actual shock
f_scenarios   <- vol_factor * z_factor_boot   # [r, S]
eps_scenarios <- vol_idio   * z_idio_boot     # [m, S]

# ── 5. Form Return Scenarios ──────────────────────────────────────────────────
common_component <- sapply(seq_len(S), function(s) {
  lambda_draws[, , s] %*% f_scenarios[, s]
})   # [m, S]

scenarios_B <- t(common_component + eps_scenarios)   # [S, m]
colnames(scenarios_B) <- TICKERS

# ── Sanity Check and Kurtosis Comparison ──────────────────────────────────────
compute_ex_kurt <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^4) / (s^4) - 3
}

# Ljung-Box diagnostics on levels and squared residuals.
# Small p-values indicate remaining serial dependence not captured by the filter.
lb_diag <- do.call(rbind, lapply(seq_len(m + r), function(i) {
  z <- as.numeric(std_resids[i, ])
  data.frame(
    series       = series_labels[i],
    kurtosis_ex  = compute_ex_kurt(z),
    lb10_p_raw   = Box.test(z, lag = 10, type = "Ljung-Box")$p.value,
    lb20_p_raw   = Box.test(z, lag = 20, type = "Ljung-Box")$p.value,
    lb10_p_sq    = Box.test(z^2, lag = 10, type = "Ljung-Box")$p.value,
    lb20_p_sq    = Box.test(z^2, lag = 20, type = "Ljung-Box")$p.value,
    stringsAsFactors = FALSE
  )
}))

cat("Ljung-Box diagnostics (p-values) on standardised residuals:\n")
print(within(lb_diag, {
  kurtosis_ex <- round(kurtosis_ex, 3)
  lb10_p_raw  <- round(lb10_p_raw, 4)
  lb20_p_raw  <- round(lb20_p_raw, 4)
  lb10_p_sq   <- round(lb10_p_sq, 4)
  lb20_p_sq   <- round(lb20_p_sq, 4)
}))
saveRDS(lb_diag, file.path(PATH_DATA, "fhs_residual_diagnostics.rds"))

cat("\nScenario summary (FHS, model B):\n")
cat("  Dimensions:", dim(scenarios_B), "\n")
cat("  Excess kurtosis comparison (FHS vs Gaussian):\n")
for (i in seq_len(m)) {
  ek_B <- compute_ex_kurt(scenarios_B[, i])
  cat(sprintf("    %-4s: FHS kurt = %.2f\n", TICKERS[i], ek_B))
}
cat("  FHS kurtosis > 0 confirms heavier tails than Gaussian model A.\n")

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(scenarios_B, file.path(PATH_DATA, "scenarios_fhs.rds"))
cat("\nScenarios saved: scenarios_fhs.rds\n")
cat("03B_scenarios_fhs.R complete.\n")
