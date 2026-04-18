# ==============================================================================
# 02_model_fit.R
# ------------------------------------------------------------------------------
# Fits the multivariate factor stochastic volatility model to the demeaned
# ETF log-return matrix using factorstochvol::fsvsample().
#
# Model (Kastner, Frühwirth-Schnatter & Lopes 2017):
#   y_t = Lambda f_t + eps_t,   t = 1,...,T
#   f_{jt} | h_{m+j,t} ~ N(0, exp(h_{m+j,t}))   [factor innovations]
#   eps_{it} | h_{it}  ~ N(0, exp(h_{it}))       [idiosyncratic innovations]
#   h_{it} = mu_i + phi_i*(h_{i,t-1} - mu_i) + sigma_i * eta_{it}  [SV AR(1)]
#
# Identification: lower-triangular Lambda with positive diagonal.
# Sampler: deep ASIS interweaving (Kastner et al. 2017, up to 400x efficiency
#          gain over the standard Gibbs sampler).
#
# Outputs (saved to PATH_DATA):
#   fsv_fit.rds           -- fsvdraws MCMC object (all post-burnin draws)
#   sv_ar1_params.rds     -- m+r x 3 matrix of posterior-mean AR(1) parameters
#                            (used by 03A/03B/03C for one-step-ahead forecasting)
#   std_residuals.rds     -- (m+r) x T matrix of standardised innovations
#                            (used by 03B and 03C)
#
# Figures (saved to PATH_CHART):
#   diag_trace_loadings.png   -- trace plots for diagonal loading parameters
#   factor_logvol_paths.png   -- posterior mean factor log-variance paths +-2 sd
#   loadings_heatmap.png      -- posterior mean loading matrix heatmap
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

returns_demeaned <- readRDS(file.path(PATH_DATA, "returns_demeaned.rds"))
T_obs <- nrow(returns_demeaned)
m     <- ncol(returns_demeaned)   # number of assets = 5

cat(sprintf("Fitting factor SV: m = %d, r = %d, T = %d\n", m, N_FACTORS, T_obs))
cat(sprintf("MCMC: %d draws, %d burn-in, thin = %d\n",
            MCMC_DRAWS, MCMC_BURNIN, MCMC_THIN))
cat("Asset order:", colnames(returns_demeaned), "\n")
cat("Lower triangular identification: Factor 1 anchored to SPY (equity),",
    "Factor 2 anchored to QQQ (growth).\n")

# ── 1. MCMC Sampling ───────────────────────────────────────────────────────────
# fsvsample() implements the full Bayesian factor SV sampler.
# interweaving = "deep": activates ASIS deep interweaving, which also
#   reparameterises the factor log-variance path, yielding the largest
#   efficiency gains (Kastner et al. 2017, Section 3.3).
# Default priors: Normal(0, B_Lambda) for loadings, standard SV priors for
#   (mu, phi, sigma) of each log-variance AR(1) process.
set.seed(2024)
# factorstochvol API compatibility:
# newer versions expect numeric interweaving codes (deep = 4),
# older versions accept character labels (e.g., "deep").
interweaving_arg <- {
  default_iw <- formals(fsvsample)[["interweaving"]]
  if (is.numeric(default_iw)) 4 else "deep"
}

fsv_fit <- fsvsample(
  y            = returns_demeaned,  # T x m demeaned return matrix
  factors      = N_FACTORS,         # number of latent factors r
  draws        = MCMC_DRAWS,
  burnin       = MCMC_BURNIN,
  thin         = MCMC_THIN,
  interweaving = interweaving_arg,  # deep ASIS: most efficient sampler
  quiet        = FALSE              # print progress to console
)

S <- MCMC_DRAWS - MCMC_BURNIN   # number of retained post-burnin draws
cat(sprintf("\nSampling complete. %d post-burnin draws retained.\n", S))

# ── 2. Save Raw MCMC Object ────────────────────────────────────────────────────
# The fsvdraws object is large (~several hundred MB for 40,000 draws).
# Downstream scripts load it on demand.
saveRDS(fsv_fit, file.path(PATH_DATA, "fsv_fit.rds"))
cat("MCMC object saved.\n")

# ── 3. Extract Posterior Summaries ────────────────────────────────────────────
# factorstochvol accessor functions:
#   loadings(fit)  -> array [m, r, S]:    factor loading draws
#   logvol(fit)    -> array [m+r, T, S]:  log-VOLATILITY draws (= h/2 in paper notation)
#                     First m series = idiosyncratic, last r = factor log-vols.
#   factors(fit)   -> array [r, T, S]:    latent factor draws

lambda_draws <- loadings(fsv_fit)           # [m, r, S]
lv_draws     <- logvol(fsv_fit)             # [m+r, T, S]  (log-volatility = h/2)
factor_draws <- factors(fsv_fit)            # [r, T, S]

# Posterior means
lambda_mean <- apply(lambda_draws, c(1, 2), mean)  # [m, r]
lv_mean     <- apply(lv_draws,     c(1, 2), mean)  # [m+r, T]
factor_mean <- apply(factor_draws, c(1, 2), mean)  # [r, T]

rownames(lambda_mean) <- TICKERS
colnames(lambda_mean) <- paste0("F", seq_len(N_FACTORS))

cat("\nPosterior mean loading matrix (Lambda):\n")
print(round(lambda_mean, 4))
cat("Interpretation: positive Factor 1 loading on SPY/QQQ = equity risk factor;\n")
cat("                negative TLT loading confirms flight-to-quality pattern.\n")

# ── 4. Compute Posterior-Mean AR(1) SV Parameters ─────────────────────────────
# Each of the m+r log-variance processes follows: h_{i,t+1} = mu_i + phi_i*(h_{i,t}-mu_i) + sigma_i*eta
# where h_{i,t} = 2 * lv_{i,t} (log-VARIANCE = 2 * log-VOLATILITY).
#
# NOTE: factorstochvol does not expose (mu, phi, sigma) draws directly via
# a standard accessor at all package versions. We therefore estimate these
# parameters from the POSTERIOR MEAN log-variance path using OLS AR(1).
# This is a practical approximation: it uses the smoothed posterior mean rather
# than the full joint posterior of (h, theta). Comment is explicit where this
# deviates from full Bayesian propagation.
#
# DEVIATION FROM FULL BAYESIAN APPROACH:
#   A fully Bayesian one-step-ahead predictive distribution would use each
#   MCMC draw s jointly: (Lambda^s, h_T^s, mu^s, phi^s, sigma^s).
#   Here we use the POSTERIOR MEAN of mu, phi, sigma (estimated below) while
#   still using per-draw h_T^s. This underestimates parameter uncertainty in
#   the predictive distribution but is tractable without access to per-draw
#   SV parameter draws.

estimate_ar1 <- function(h_path) {
  # OLS AR(1) regression on a single log-variance time series.
  # Returns named vector c(mu, phi, sigma).
  T_h   <- length(h_path)
  mu    <- mean(h_path)                              # unconditional mean
  y_lag <- h_path[1:(T_h - 1)] - mu
  y_cur <- h_path[2:T_h]       - mu
  phi   <- sum(y_cur * y_lag) / sum(y_lag^2)        # OLS slope
  phi   <- max(-0.999, min(0.999, phi))              # clamp within stationarity
  resid <- y_cur - phi * y_lag                       # AR(1) residuals
  sigma <- sd(resid)                                 # innovation std dev
  c(mu = mu, phi = phi, sigma = sigma)
}

# Log-variance = 2 * log-volatility (paper convention: h = log(sigma^2))
logvar_mean <- 2 * lv_mean    # [m+r, T] posterior mean log-variance paths

# Fit AR(1) to each of the m+r series using the posterior mean path
sv_ar1_params <- t(apply(logvar_mean, 1, estimate_ar1))   # [m+r, 3]
rownames(sv_ar1_params) <- c(TICKERS, paste0("Factor", seq_len(N_FACTORS)))
colnames(sv_ar1_params) <- c("mu", "phi", "sigma")

cat("\nPosterior-mean AR(1) SV parameters:\n")
print(round(sv_ar1_params, 4))
cat("Note: mu = log-variance level; phi = persistence (close to 1 = high clustering);\n")
cat("      sigma = vol-of-vol (innovation std dev of log-variance process).\n")

saveRDS(sv_ar1_params, file.path(PATH_DATA, "sv_ar1_params.rds"))

# ── 5. Standardised Residuals ──────────────────────────────────────────────────
# Standardised factor innovations: f_{j,t} / exp(lv_{m+j, t})
#   Under the model these should be iid N(0,1).
# Standardised idiosyncratic: eps_{i,t} / exp(lv_{i,t})
#   eps_{i,t} = y_{i,t} - Lambda_mean %*% factor_mean[,t]
#
# These are used by 03B (FHS resampling) and 03C (EVT fitting).
# We use POSTERIOR MEAN Lambda and factors (not full posterior) for tractability.

idio_mean <- t(returns_demeaned) - lambda_mean %*% factor_mean   # [m, T]

# Posterior mean log-volatilities, separated into idiosyncratic and factor parts
lv_idio   <- lv_mean[1:m, ]          # [m, T] idiosyncratic log-vols
lv_factor <- lv_mean[(m+1):(m+N_FACTORS), ]  # [r, T] factor log-vols

# Standardise: divide by posterior mean conditional volatility
std_idio   <- idio_mean   / exp(lv_idio)    # [m, T]
std_factor <- factor_mean / exp(lv_factor)  # [r, T]

# Stack: first m rows = idiosyncratic series, last r rows = factor series
std_residuals <- rbind(std_idio, std_factor)   # [m+r, T]
rownames(std_residuals) <- rownames(sv_ar1_params)

cat("\nStandardised residual summary (should be near N(0,1) under the model):\n")
cat("  Column means:", round(rowMeans(std_residuals), 3), "\n")
cat("  Column SDs:  ", round(apply(std_residuals, 1, sd), 3), "\n")
cat("  Excess kurtosis:\n")
ex_kurt <- apply(std_residuals, 1, function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^4) / (s^4) - 3
})
print(round(ex_kurt, 3))
cat("  Excess kurtosis >> 0 confirms fat-tailed innovations: motivates FHS and EVT.\n")

saveRDS(std_residuals, file.path(PATH_DATA, "std_residuals.rds"))

# ── 6. Inefficiency Factors ────────────────────────────────────────────────────
# Inefficiency Factor (IF) = integrated autocorrelation time.
# IF ~ 1 means efficient sampling; IF >> 1 means slow mixing.
# Deep interweaving should give IF << 100 for diagonal loadings.
if (requireNamespace("coda", quietly = TRUE)) {
  cat("\nInefficiency Factors (diagonal loadings) -- lower is better:\n")
  for (j in seq_len(N_FACTORS)) {
    chain  <- lambda_draws[j, j, ]                      # diagonal element [j,j]
    n_eff  <- coda::effectiveSize(chain)
    IF_val <- S / n_eff
    cat(sprintf("  lambda[%d,%d] (%-3s, Factor %d): IF = %.1f\n",
                j, j, TICKERS[j], j, IF_val))
  }
}

# ── 7. Diagnostic Figures ─────────────────────────────────────────────────────

# --- (a) Trace plots: diagonal loadings ---
png(file.path(PATH_CHART, "diag_trace_loadings.png"),
    width = 1400, height = 500, res = 120)
par(mfrow = c(1, N_FACTORS), mar = c(4, 4, 3, 1))
for (j in seq_len(N_FACTORS)) {
  chain <- lambda_draws[j, j, ]
  plot(chain, type = "l", col = "steelblue",
       main = bquote(paste("Trace: ", lambda[.(j)*","*.(j)],
                           " (", .(TICKERS[j]), ", Factor ", .(j), ")")),
       xlab = "Draw (post burn-in)", ylab = "Loading value", lwd = 0.4)
  abline(h = mean(chain), col = "tomato", lty = 2, lwd = 1.5)  # posterior mean
}
dev.off()

# --- (b) Factor log-variance paths: posterior mean +-2 sd ---
dates_ret  <- as.Date(rownames(returns_demeaned))
lv_sd      <- apply(lv_draws, c(1, 2), sd)   # [m+r, T] posterior sd of log-vol

png(file.path(PATH_CHART, "factor_logvol_paths.png"),
    width = 1400, height = 450 * N_FACTORS, res = 120)
par(mfrow = c(N_FACTORS, 1), mar = c(4, 4, 3, 1))
for (j in seq_len(N_FACTORS)) {
  idx   <- m + j                                  # factor series index in [m+r]
  h_mu  <- 2 * lv_mean[idx, ]                     # posterior mean log-variance
  h_sd  <- 2 * lv_sd[idx, ]                       # posterior sd  log-variance
  ylim  <- range(h_mu + 2*h_sd, h_mu - 2*h_sd)
  plot(dates_ret, h_mu, type = "l", col = "steelblue", ylim = ylim,
       main = paste0("Factor ", j, " log-variance: posterior mean \u00b1 2 sd"),
       xlab = "Date", ylab = "log(variance)", lwd = 1.2)
  polygon(c(dates_ret, rev(dates_ret)),
          c(h_mu + 2*h_sd, rev(h_mu - 2*h_sd)),
          col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)
}
dev.off()

# --- (c) Posterior mean loading matrix heatmap ---
lm_df <- melt(lambda_mean)
colnames(lm_df) <- c("Asset", "Factor", "Loading")
p_heat <- ggplot(lm_df, aes(x = Factor, y = Asset, fill = Loading)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Loading, 3)), size = 3.5) +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "steelblue",
                       midpoint = 0, name = "Loading") +
  labs(title = "Posterior Mean Factor Loading Matrix",
       subtitle = "Lower triangular identification; positive diagonal enforced") +
  theme_minimal(base_size = 12)
ggsave(file.path(PATH_CHART, "loadings_heatmap.png"),
       p_heat, width = 6, height = 4, dpi = 150)

cat("\nDiagnostic figures saved to:", PATH_CHART, "\n")
cat("02_model_fit.R complete.\n")
