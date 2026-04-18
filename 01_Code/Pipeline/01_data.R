# ==============================================================================
# 01_data.R
# ------------------------------------------------------------------------------
# Downloads adjusted ETF prices from Yahoo Finance, computes log returns,
# demeanes each series, and saves the cleaned data for downstream scripts.
#
# Outputs (saved to PATH_DATA):
#   prices.rds         -- aligned adjusted-price xts object
#   log_returns.rds    -- raw log-return xts object
#   returns_demeaned.rds -- T x m matrix, mean-removed (model input)
#   col_means.rds      -- m-vector of removed daily means
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

# ── 1. Download Adjusted Prices ────────────────────────────────────────────────
# getSymbols() writes each ticker into env_prices as an xts object.
# Ad() extracts the adjusted-close column (accounts for dividends and splits).
cat("Downloading data for:", paste(TICKERS, collapse = ", "), "\n")

env_prices <- new.env()
getSymbols(
  TICKERS,
  src         = "yahoo",
  from        = DATE_START,
  to          = DATE_END,
  auto.assign = TRUE,
  env         = env_prices
)

# Pull adjusted closes into a list then merge on the common date index
prices_list <- lapply(TICKERS, function(tk) Ad(get(tk, envir = env_prices)))
names(prices_list) <- TICKERS

# merge(..., join = "inner") keeps only days where ALL five ETFs traded
prices_xts <- do.call(merge, c(prices_list, join = "inner"))
colnames(prices_xts) <- TICKERS

cat(sprintf("Aligned price series: %d obs  [%s  to  %s]\n",
            nrow(prices_xts),
            format(start(prices_xts), "%Y-%m-%d"),
            format(end(prices_xts),   "%Y-%m-%d")))

# ── 2. Log Returns ─────────────────────────────────────────────────────────────
# diff(log(p_t)) = log(p_t / p_{t-1}): continuously compounded daily returns.
log_returns_xts <- diff(log(prices_xts))
log_returns_xts <- na.omit(log_returns_xts)   # drop the NA in the first row

cat(sprintf("Log-return series: %d obs\n", nrow(log_returns_xts)))

# Convert to plain T x m matrix for factorstochvol (which expects a matrix)
returns_mat <- as.matrix(log_returns_xts)

# ── 3. Demean ──────────────────────────────────────────────────────────────────
# The factor SV model is specified for CENTRED returns.
# Equity ETFs carry a positive drift (equity risk premium ~5-8% p.a.) which,
# if left in, inflates the log-variance level estimates (mu parameter).
# We remove the in-sample column mean, which is the minimum-assumption demean.
#
# NOTE: col_means.rds is saved for reference / future use.
# The backtest in 04_risk_budgeting.R uses the RAW log_returns_xts (not
# demeaned), so the drift is automatically included in backtest returns.
col_means        <- colMeans(returns_mat)
returns_demeaned <- sweep(returns_mat, 2, col_means, "-")

cat("\nSample means removed (annualised %):\n")
print(round(col_means * 252 * 100, 2))

# Sanity check: verify the column means of the demeaned matrix are ~0
stopifnot(all(abs(colMeans(returns_demeaned)) < 1e-12))

# ── 4. Descriptive Statistics ─────────────────────────────────────────────────
# Summary table: annualised vol, skewness, excess kurtosis.
# Excess kurtosis >> 0 for all series confirms non-Gaussian tails (motivates
# the FHS and EVT approaches in 03B and 03C).
moments_skew <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^3) / (s^3)
}

moments_ex_kurt <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^4) / (s^4) - 3
}

desc_stats <- data.frame(
  ticker   = TICKERS,
  ann_vol  = round(apply(returns_mat, 2, sd) * sqrt(252) * 100, 2),
  skewness = round(apply(returns_mat, 2, moments_skew), 3),
  ex_kurt  = round(apply(returns_mat, 2, moments_ex_kurt), 3)
)
cat("\nDescriptive statistics (raw log returns):\n")
print(desc_stats)
cat("  Excess kurtosis > 0 for all series: Gaussian innovations will underestimate tail risk.\n")

# ── 5. Save ────────────────────────────────────────────────────────────────────
saveRDS(prices_xts,        file.path(PATH_DATA, "prices.rds"))
saveRDS(log_returns_xts,   file.path(PATH_DATA, "log_returns.rds"))
saveRDS(returns_demeaned,  file.path(PATH_DATA, "returns_demeaned.rds"))
saveRDS(col_means,         file.path(PATH_DATA, "col_means.rds"))
saveRDS(desc_stats,        file.path(PATH_DATA, "desc_stats.rds"))

cat("\nData saved to:", PATH_DATA, "\n")
