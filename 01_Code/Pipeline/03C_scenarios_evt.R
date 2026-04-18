# ==============================================================================
# 03C_scenarios_evt.R
# ------------------------------------------------------------------------------
# MODEL C: Semi-parametric EVT tail model (robustness check).
#
# DEVIATION FROM GAUSSIAN ASSUMPTION (Model A):
#   Extreme Value Theory provides a principled tail model for rare events.
#   Rather than resampling the entire empirical distribution (Model B / FHS),
#   we fit a GENERALISED PARETO DISTRIBUTION (GPD) to the LOWER TAIL of each
#   standardised residual series and use a SEMI-PARAMETRIC distribution:
#
#     - Body  (> threshold u): empirical distribution (as in FHS)
#     - Tails (< threshold u): GPD extrapolation
#
#   This allows smooth extrapolation beyond the historical sample for extreme
#   quantiles, making ES estimates more robust to limited tail data.
#   The GPD is appropriate by the Pickands-Balkema-de Haan theorem:
#   exceedances above a high threshold converge in distribution to GPD.
#
#   Reference: Embrechts, Klüppelberg & Mikosch (1997); McNeil & Frey (2000).
#
# DEVIATION NOTE:
#   GPD is fitted independently to each of the m+r standardised residual series.
#   Cross-sectional dependence in the tails (co-crashes) is preserved via JOINT
#   resampling of the body draws (same as FHS), with the tail part handled
#   marginally. A full copula model would be the next step beyond this.
#
# Output (saved to PATH_DATA):
#   scenarios_evt.rds     -- N_SCENARIOS x m matrix of return scenarios
#   gpd_params.rds        -- (m+r) x 4 data frame of GPD fit parameters
#
# Figures (saved to PATH_CHART):
#   mean_excess_plots.png  -- mean-excess plots for threshold validation
#   gpd_tail_fits.png      -- empirical vs GPD tail probability
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

cat(sprintf("Model C (EVT): generating %d scenarios for m=%d assets.\n", S, m))
cat("DEVIATION: Lower tail replaced by GPD extrapolation; body resampled empirically.\n")

# ── Helper: Fit GPD to Lower Tail ─────────────────────────────────────────────
fit_gpd_lower_tail <- function(z, q_thresh = EVT_THRESHOLD_Q) {
  # Fits GPD to the left tail of z (negated, so we work with upper exceedances).
  # z:        standardised residual vector (zero-mean, unit-scale)
  # q_thresh: empirical quantile below which GPD is fitted
  # Returns: list(threshold, xi, beta, n_exceedances, cdf_at_threshold)

  neg_z    <- -z                                   # flip sign: lower tail → upper
  u        <- quantile(neg_z, 1 - q_thresh)        # threshold on the negated series
  exceed   <- neg_z[neg_z > u] - u                 # exceedances above threshold

  if (length(exceed) < 10) {
    warning("Too few exceedances for GPD fit; reverting to empirical tail.")
    return(NULL)
  }

  # evir::gpd fits by maximum likelihood; shape xi and scale beta are returned
  gpd_fit  <- evir::gpd(neg_z, threshold = u)
  xi       <- gpd_fit$par.ests["xi"]    # shape: xi > 0 → heavy tail (Pareto-like)
  beta     <- gpd_fit$par.ests["beta"]  # scale

  list(
    threshold   = u,
    xi          = xi,
    beta        = beta,
    n_exceed    = length(exceed),
    p_exceed    = mean(neg_z > u),   # empirical exceedance probability
    u_original  = -u                 # threshold on original (un-negated) scale
  )
}

# ── Helper: Draw from GPD Lower Tail ──────────────────────────────────────────
# Transforms uniform(0, p_exceed) to GPD quantiles, then negates back.
draw_gpd_tail <- function(n_draws, gpd_params) {
  # Inverse CDF of GPD (Pickands' representation):
  # Q(p) = u + (beta/xi) * ((1-p)^(-xi) - 1)  for xi != 0
  # Q(p) = u - beta * log(1-p)                 for xi == 0

  u    <- gpd_params$threshold
  xi   <- gpd_params$xi
  beta <- gpd_params$beta
  p_ex <- gpd_params$p_exceed

  # Draw uniform on [0, p_exceed] then map through conditional GPD quantile
  u_unif <- runif(n_draws, min = 0, max = p_ex)

  if (abs(xi) < 1e-6) {
    # xi ≈ 0: exponential tail
    tail_vals <- u + beta * (-log(1 - u_unif / p_ex))
  } else {
    tail_vals <- u + (beta / xi) * ((1 - u_unif / p_ex)^(-xi) - 1)
  }

  # Convert back to original scale (flip sign: was upper tail of -z)
  -tail_vals
}

# ── 1. Fit GPD to Lower Tail of Each Standardised Residual Series ─────────────
series_labels <- rownames(std_resids)
gpd_fits      <- vector("list", m + r)
names(gpd_fits) <- series_labels

gpd_param_df <- data.frame(
  series   = series_labels,
  threshold = NA_real_,
  xi        = NA_real_,   # shape (> 0 = heavy tail, 0 = exponential, < 0 = bounded)
  beta      = NA_real_,   # scale
  n_exceed  = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(m + r)) {
  fit_i <- fit_gpd_lower_tail(std_resids[i, ], q_thresh = EVT_THRESHOLD_Q)
  if (!is.null(fit_i)) {
    gpd_fits[[i]]            <- fit_i
    gpd_param_df$threshold[i] <- fit_i$u_original   # on original scale
    gpd_param_df$xi[i]        <- fit_i$xi
    gpd_param_df$beta[i]      <- fit_i$beta
    gpd_param_df$n_exceed[i]  <- fit_i$n_exceed
    cat(sprintf("  %-12s: xi = %6.3f, beta = %6.3f, n_exceed = %d%s\n",
                series_labels[i], fit_i$xi, fit_i$beta, fit_i$n_exceed,
                ifelse(fit_i$xi > 0, "  [heavy tail]", "")))
  } else {
    cat(sprintf("  %-12s: GPD fit failed; will use empirical tail.\n",
                series_labels[i]))
  }
}

cat("\nGPD summary -- xi > 0 confirms heavy tails (validates EVT approach):\n")
print(round(gpd_param_df[, c("series", "xi", "beta", "n_exceed")], 4))
saveRDS(gpd_param_df, file.path(PATH_DATA, "gpd_params.rds"))

# ── 2. Mean-Excess Plots (Threshold Validation) ────────────────────────────────
# The mean excess function E[X - u | X > u] should be LINEAR in u under GPD.
# A straight line starting from our chosen threshold u confirms GPD is appropriate.
png(file.path(PATH_CHART, "mean_excess_plots.png"),
    width = 1600, height = 1200, res = 120)
par(mfrow = c(ceiling((m + r) / 3), 3), mar = c(4, 4, 3, 1))
for (i in seq_len(m + r)) {
  neg_z     <- -std_resids[i, ]         # negate: work with upper tail
  sorted_z  <- sort(neg_z)
  u_grid    <- sorted_z[seq(1, floor(0.90 * length(sorted_z)), by = 5)]
  me_vals   <- sapply(u_grid, function(u) mean(neg_z[neg_z > u] - u))
  plot(u_grid, me_vals, type = "l", col = "steelblue",
       main = paste0("Mean Excess: ", series_labels[i]),
       xlab = "Threshold u", ylab = "Mean Excess")
  # Mark the chosen threshold
  thresh_chosen <- -gpd_param_df$threshold[i]   # back to negated scale
  abline(v = thresh_chosen, col = "tomato", lty = 2, lwd = 1.5)
  legend("topleft", legend = "Chosen threshold", col = "tomato",
         lty = 2, bty = "n", cex = 0.7)
}
dev.off()
cat("Mean-excess plots saved: mean_excess_plots.png\n")

# GPD tail fit diagnostics: empirical vs fitted survival in the exceedance domain.
png(file.path(PATH_CHART, "gpd_tail_fits.png"),
    width = 1600, height = 1200, res = 120)
par(mfrow = c(ceiling((m + r) / 3), 3), mar = c(4, 4, 3, 1))
for (i in seq_len(m + r)) {
  fit_i <- gpd_fits[[i]]
  if (is.null(fit_i)) {
    plot.new()
    title(main = paste0("GPD Tail Fit: ", series_labels[i]))
    text(0.5, 0.5, "Fit unavailable\n(too few exceedances)")
    next
  }

  neg_z <- -std_resids[i, ]
  u     <- fit_i$threshold
  xi    <- fit_i$xi
  beta  <- fit_i$beta
  y_exc <- sort(neg_z[neg_z > u] - u)

  n_exc    <- length(y_exc)
  emp_surv <- rev(seq_len(n_exc)) / n_exc
  fit_surv <- if (abs(xi) < 1e-6) {
    exp(-y_exc / beta)
  } else {
    pmax((1 + xi * y_exc / beta)^(-1 / xi), 1e-10)
  }

  plot(y_exc, emp_surv, log = "y", type = "p",
       pch = 16, cex = 0.55, col = adjustcolor("steelblue", alpha.f = 0.65),
       xlab = "Exceedance above threshold", ylab = "Tail survival P(Y > y)",
       main = paste0("GPD Tail Fit: ", series_labels[i]))
  lines(y_exc, fit_surv, col = "tomato", lwd = 2)
  legend("topright", legend = c("Empirical tail", "Fitted GPD"),
         col = c("steelblue", "tomato"), pch = c(16, NA), lty = c(NA, 1),
         bty = "n", cex = 0.75)
}
dev.off()
cat("GPD tail-fit diagnostics saved: gpd_tail_fits.png\n")

# ── 3. Log-Variance Forecasts (identical to Models A and B) ───────────────────
lv_draws     <- logvol(fsv_fit)
lambda_draws <- loadings(fsv_fit)
logvar_T     <- 2 * lv_draws[, T_obs, ]

mu_vec    <- sv_params[, "mu"]
phi_vec   <- sv_params[, "phi"]
sigma_vec <- sv_params[, "sigma"]

set.seed(44)
z_sv      <- matrix(rnorm((m + r) * S), nrow = m + r, ncol = S)
logvar_T1 <- mu_vec + phi_vec * (logvar_T - mu_vec) + sigma_vec * z_sv
vol_T1    <- exp(0.5 * logvar_T1)

vol_idio   <- vol_T1[1:m, ]
vol_factor <- vol_T1[(m+1):(m+r), ]

# ── 4. Semi-Parametric Draw: Body Resampled, Tail from GPD ───────────────────
# DEVIATION FROM MODEL B:
#   In FHS we resample from the full empirical distribution.
#   Here the body draw is the same, but draws that FALL INTO THE TAIL REGION
#   are REPLACED by GPD draws. This is the semi-parametric combination.
#
#   Implementation:
#   (a) Draw body index from historical dates (same as FHS).
#   (b) For each series i and draw s, check if the historical draw is in the tail.
#       If z_boot[i,s] < threshold_i: replace with a GPD draw instead.
#   This gives a hybrid distribution: empirical for moderate observations,
#   GPD-parametric for extremes.

T_hist   <- ncol(std_resids)
boot_idx <- sample(seq_len(T_hist), size = S, replace = TRUE)
z_boot   <- std_resids[, boot_idx]   # [m+r, S]  initial bootstrap draws

# For each series, replace tail observations with GPD draws
z_evt <- z_boot   # start from FHS draws, then overwrite tail elements

for (i in seq_len(m + r)) {
  if (is.null(gpd_fits[[i]])) next    # no GPD fit; keep empirical draw
  thresh_i  <- gpd_fits[[i]]$u_original   # threshold on original (un-negated) scale
  tail_mask <- z_evt[i, ] < thresh_i      # TRUE where the draw is in the left tail

  n_tail_draws <- sum(tail_mask)
  if (n_tail_draws > 0) {
    z_evt[i, tail_mask] <- draw_gpd_tail(n_tail_draws, gpd_fits[[i]])
  }
}

# Split into factor and idiosyncratic parts
z_factor_evt <- z_evt[(m+1):(m+r), ]
z_idio_evt   <- z_evt[1:m, ]

# ── 5. Scale and Form Return Scenarios ────────────────────────────────────────
f_scenarios   <- vol_factor * z_factor_evt
eps_scenarios <- vol_idio   * z_idio_evt

common_component <- sapply(seq_len(S), function(s) {
  lambda_draws[, , s] %*% f_scenarios[, s]
})   # [m, S]

scenarios_C <- t(common_component + eps_scenarios)   # [S, m]
colnames(scenarios_C) <- TICKERS

# ── 6. Tail Comparison: GPD vs FHS vs Gaussian ────────────────────────────────
scenarios_A <- readRDS(file.path(PATH_DATA, "scenarios_gaussian.rds"))
scenarios_B <- readRDS(file.path(PATH_DATA, "scenarios_fhs.rds"))

cat("\nLeft-tail comparison (5th percentile of each ETF's scenario distribution):\n")
cat(sprintf("  %-6s  %10s  %10s  %10s\n", "ETF", "Gaussian", "FHS", "EVT"))
for (i in seq_len(m)) {
  q5_A <- quantile(scenarios_A[, i], 0.05)
  q5_B <- quantile(scenarios_B[, i], 0.05)
  q5_C <- quantile(scenarios_C[, i], 0.05)
  cat(sprintf("  %-6s  %10.5f  %10.5f  %10.5f\n",
              TICKERS[i], q5_A, q5_B, q5_C))
}
cat("  More negative 5th percentile in EVT/FHS vs Gaussian confirms heavier tails.\n")

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(scenarios_C, file.path(PATH_DATA, "scenarios_evt.rds"))
cat("\nScenarios saved: scenarios_evt.rds\n")
cat("03C_scenarios_evt.R complete.\n")
