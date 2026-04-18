# ==============================================================================
# 00_config.R
# ------------------------------------------------------------------------------
# Single source of truth for all pipeline parameters.
# Every downstream script begins with: source("00_config.R")
# ==============================================================================

# ── Packages ───────────────────────────────────────────────────────────────────
REQUIRED_PACKAGES <- c(
  "quantmod",       # Yahoo Finance data download
  "factorstochvol", # Multivariate factor SV model (Hosszejni & Kastner 2021)
  "evir",           # GPD fitting for EVT (Embrechts et al.)
  "ggplot2",        # Plotting
  "patchwork",      # Combining ggplot panels side by side
  "reshape2",       # melt() for long-format plotting
  "lubridate",      # Date arithmetic (annual rebalancing logic)
  "xts"             # Time-series alignment and indexing
)

missing_pkgs <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
      "Install with:\n",
      "install.packages(c(", paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))"
    ),
    call. = FALSE
  )
}

invisible(lapply(REQUIRED_PACKAGES, library, character.only = TRUE))

# ── Directory Paths ────────────────────────────────────────────────────────────
PATH_ROOT  <- "C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models"
PATH_CODE  <- file.path(PATH_ROOT, "01_Code/Pipeline")
PATH_DATA  <- file.path(PATH_ROOT, "02_Data")
PATH_CHART <- file.path(PATH_ROOT, "03_Charts")

dir.create(PATH_DATA,  showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_CHART, showWarnings = FALSE, recursive = TRUE)

# ── ETF Universe ───────────────────────────────────────────────────────────────
# Five ETFs spanning equity, fixed income, and gold.
TICKERS <- c("SPY", "QQQ", "TLT", "SHY", "GLD")

# Fixed intra-bucket allocation ratios (must sum to 1.0 within each bucket).
# These are HELD CONSTANT across all strategies; only the bucket-level split varies.
BUCKET_WEIGHTS <- list(
  equity = c(SPY = 0.80, QQQ = 0.20),  # US large-cap 80%, US tech/growth 20%
  bonds  = c(TLT = 0.50, SHY = 0.50),  # Long Treasury 50%, Short Treasury 50%
  gold   = c(GLD = 1.00)               # Gold 100%
)

# ── Sample Period ──────────────────────────────────────────────────────────────
# GLD launched Nov 2004; 2006 gives a clean full-year start with all ETFs live.
DATE_START <- "2006-01-03"
DATE_END   <- Sys.Date()

# ── Factor SV Model Settings ───────────────────────────────────────────────────
# Asset ordering matters: lower triangular identification anchors factor j to
# asset j. Placing SPY first makes Factor 1 the equity risk factor; QQQ second
# makes Factor 2 the residual growth/tech factor.
# With m = 5 assets and r = 2 factors the model is well-identified.
N_FACTORS   <- 2L

MCMC_DRAWS  <- 50000L   # Total MCMC iterations (includes burn-in)
MCMC_BURNIN <- 10000L   # Iterations discarded as burn-in
MCMC_THIN   <- 1L       # Thinning interval (1 = keep every draw post burn-in)

# ── ES Parameters ─────────────────────────────────────────────────────────────
# ES_95: Expected Shortfall at 95% confidence = expected loss in worst 5% of days
ES_LEVEL <- 0.05   # Tail probability alpha (left tail)

# ── Risk Budget Targets ────────────────────────────────────────────────────────
# Each budget vector specifies the TARGET fraction of portfolio ES to be borne
# by each asset-class bucket. Must sum to 1.
#
# These same fractions are ALSO used as the static strategy weights, so that
# the two approaches are directly comparable.
#
# Example: pf1 targets equal equity/bond risk with a 10% gold allocation.
BUDGETS <- list(
  pf1 = c(equity = 0.45, bonds = 0.45, gold = 0.10),
  pf2 = c(equity = 0.50, bonds = 0.50, gold = 0.00)
)

# ── Static Strategy ────────────────────────────────────────────────────────────
# The static strategy uses BUDGET PERCENTAGES DIRECTLY as portfolio weights
# (i.e. 45% equity bucket, 45% bond bucket, 10% gold) and rebalances annually.
# This is the benchmark against the risk-contribution-targeted dynamic strategy.
REBALANCE_FREQ <- "annual"   # "annual" supported; extend to "quarterly" if needed

# Dynamic strategy: rebalance only when realised ES risk contributions move
# outside target bands (plus a cool-down to avoid churn).
RC_LOOKBACK_DAYS    <- 252L  # Rolling window used to estimate realised ES RC
RC_BAND_LOWER       <- 0.03  # Allowed downside deviation from target RC share
RC_BAND_UPPER       <- 0.03  # Allowed upside   deviation from target RC share
RC_MIN_DAYS_BETWEEN <- 21L   # Minimum gap between consecutive rebalances

# ── EVT Settings ──────────────────────────────────────────────────────────────
# Threshold quantile for GPD fitting: fit the tail below this empirical quantile
# of the standardised residuals. Mean-excess plot should be consulted to validate.
EVT_THRESHOLD_Q <- 0.10   # Lower 10% of standardised residual distribution

# ── Optimiser Settings ─────────────────────────────────────────────────────────
OPT_TOL     <- 1e-9    # Convergence tolerance for optim()
OPT_MAXITER <- 5000L   # Maximum number of optimiser iterations

# ── Scenario Count ─────────────────────────────────────────────────────────────
# Number of one-step-ahead scenarios generated for ES estimation.
# Using all post-burnin MCMC draws (one scenario per draw) ensures full
# propagation of parameter uncertainty into the ES estimate.
N_SCENARIOS <- MCMC_DRAWS - MCMC_BURNIN   # = 40 000

# ── Helper: Expand Bucket Weights to Asset Weights ─────────────────────────────
# Converts a named bucket-weight vector (length 3) into a named asset-weight
# vector (length 5) using the fixed intra-bucket ratios in BUCKET_WEIGHTS.
bucket_to_asset_weights <- function(w_buckets) {
  w_assets <- numeric(length(TICKERS))
  names(w_assets) <- TICKERS
  for (b in names(BUCKET_WEIGHTS)) {
    for (asset in names(BUCKET_WEIGHTS[[b]])) {
      w_assets[asset] <- w_buckets[b] * BUCKET_WEIGHTS[[b]][asset]
    }
  }
  w_assets
}

cat("Config loaded. Pipeline root:", PATH_ROOT, "\n")
