# ==============================================================================
# 04_risk_budgeting.R
# ------------------------------------------------------------------------------
# Shared optimization engine. Runs for all 3 scenario sets x 2 portfolios = 6
# dynamic allocations, plus 2 static allocations as the benchmark.
#
# DYNAMIC STRATEGY (risk-contribution targeted):
#   Finds bucket weights w* = (w_equity, w_bonds, w_gold) such that the ES
#   risk contribution of each bucket matches the target budget:
#     RC_bucket_j^ES(w*) / ES(w*) = b_j   for j in {equity, bonds, gold}
#   where intra-bucket weights are FIXED (SPY/QQQ split, TLT/SHY split).
#
#   Objective (minimised):  sum_j (RC_j/ES - b_j)^2
#   Method: optim() L-BFGS-B on softmax-reparameterised bucket weights.
#   Long-only: enforced via softmax (weights always > 0) and zero-budget
#              buckets forced to 0 before renormalisation.
#
# STATIC STRATEGY (comparison benchmark):
#   Uses BUDGET PERCENTAGES DIRECTLY as portfolio weights (no optimisation).
#   E.g. b = (0.45, 0.45, 0.10) → w_equity = 45%, w_bonds = 45%, w_gold = 10%.
#   Rebalanced annually: on the first trading day of each calendar year, weights
#   are reset to the static target.
#
# ES DECOMPOSITION (Euler / Tasche 1999):
#   ES_alpha(r_p) = sum_i RC_i^ES,   where RC_i = w_i * E[-r_i | r_p <= VaR_alpha]
#   The Euler decomposition holds because ES is positively homogeneous of degree 1.
#   At the optimal weights, each bucket's contribution equals b_j * ES.
#
# Output (saved to PATH_DATA):
#   results_dynamic.rds  -- list: one entry per model variant x portfolio
#   results_static.rds   -- list: one entry per portfolio
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

log_returns_xts <- readRDS(file.path(PATH_DATA, "log_returns.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Core ES and Risk Contribution Functions
# ══════════════════════════════════════════════════════════════════════════════

compute_es_contributions <- function(scenarios, w_assets, alpha = ES_LEVEL) {
  # ---------------------------------------------------------------------------
  # Compute ES and Euler asset-level risk contributions from a scenario matrix.
  #
  # Args:
  #   scenarios : N x m matrix of return scenarios (rows = scenarios, cols = assets)
  #   w_assets  : m-vector of ASSET-level portfolio weights (must sum to 1)
  #   alpha     : left-tail probability (default ES_LEVEL = 0.05)
  #
  # Returns: list(es, var, rc_asset, n_tail)
  #   es       : scalar, Expected Shortfall (positive = loss)
  #   var      : scalar, Value at Risk (positive = loss)
  #   rc_asset : m-vector, Euler risk contribution per asset (sum = es)
  # ---------------------------------------------------------------------------
  pf_ret    <- as.numeric(scenarios %*% w_assets)   # N-vector of portfolio returns
  var_val   <- quantile(pf_ret, alpha)               # alpha-quantile (negative loss)
  tail_mask <- pf_ret <= var_val                     # indicator: in-tail scenarios
  n_tail    <- sum(tail_mask)

  es <- -mean(pf_ret[tail_mask])   # ES = -E[r_p | r_p <= VaR]; positive = loss

  # Euler decomposition: RC_i = w_i * E[-r_i | r_p <= VaR_alpha]
  tail_ret  <- scenarios[tail_mask, , drop = FALSE]   # n_tail x m
  rc_asset  <- -w_assets * colMeans(tail_ret)         # m-vector

  list(es = es, var = -var_val, rc_asset = rc_asset, n_tail = n_tail)
}


aggregate_to_buckets <- function(rc_asset) {
  # Sums asset-level risk contributions to the three bucket levels.
  # Uses BUCKET_WEIGHTS from config to identify which assets belong to each bucket.
  sapply(names(BUCKET_WEIGHTS), function(b) {
    sum(rc_asset[names(BUCKET_WEIGHTS[[b]])])
  })
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Dynamic Risk Budget Optimiser
# ══════════════════════════════════════════════════════════════════════════════

risk_budget_objective <- function(theta, scenarios, target_b, alpha) {
  # ---------------------------------------------------------------------------
  # Objective function for the risk budgeting optimiser.
  # theta: unconstrained R^3 vector; softmax maps to the 2-simplex.
  # ---------------------------------------------------------------------------
  # Softmax ensures sum(w) = 1 and w > 0 for all buckets
  w_buckets       <- exp(theta) / sum(exp(theta))
  names(w_buckets) <- names(target_b)

  # Force zero-budget buckets to exactly 0, then renormalise
  zero_b           <- target_b == 0
  w_buckets[zero_b] <- 0
  total_w          <- sum(w_buckets)
  if (total_w < 1e-10) return(1e6)
  w_buckets        <- w_buckets / total_w

  # Expand to asset weights using fixed intra-bucket ratios
  w_assets <- bucket_to_asset_weights(w_buckets)

  # Compute ES and Euler decomposition
  res       <- compute_es_contributions(scenarios, w_assets, alpha)
  rc_bucket <- aggregate_to_buckets(res$rc_asset)

  # Relative contributions (fraction of total ES borne by each bucket)
  rel_rc    <- rc_bucket / res$es

  # Penalise deviation from target only for non-zero-budget buckets
  active    <- !zero_b
  sum((rel_rc[active] - target_b[active])^2)
}


optimise_risk_budget <- function(scenarios, target_b, alpha = ES_LEVEL,
                                 n_restarts = 5L) {
  # ---------------------------------------------------------------------------
  # Finds optimal bucket weights minimising the risk-budget objective.
  # Multiple random restarts reduce sensitivity to local minima.
  #
  # Returns: list(w_buckets, w_assets, es, rc_asset, rc_bucket, rel_rc, conv)
  # ---------------------------------------------------------------------------
  best_val <- Inf
  best_res <- NULL

  for (k in seq_len(n_restarts)) {
    theta0 <- if (k == 1) rep(0, 3) else rnorm(3, 0, 0.5)  # 1st: equal start

    opt <- tryCatch(
      optim(
        par     = theta0,
        fn      = risk_budget_objective,
        scenarios = scenarios,
        target_b  = target_b,
        alpha     = alpha,
        method  = "L-BFGS-B",
        lower   = rep(-20, 3),  # keeps exp(theta) numerically stable
        upper   = rep( 20, 3),
        control = list(maxit = OPT_MAXITER, reltol = OPT_TOL)
      ),
      error = function(e) list(value = Inf, par = theta0, convergence = 99)
    )

    if (opt$value < best_val) {
      best_val <- opt$value
      best_res <- opt
    }
  }

  # Recover final weights from best solution
  theta_star    <- best_res$par
  w_buckets     <- exp(theta_star) / sum(exp(theta_star))
  names(w_buckets) <- names(target_b)

  zero_b           <- target_b == 0
  w_buckets[zero_b] <- 0
  w_buckets        <- w_buckets / sum(w_buckets)

  w_assets  <- bucket_to_asset_weights(w_buckets)
  res_final <- compute_es_contributions(scenarios, w_assets)
  rc_bucket <- aggregate_to_buckets(res_final$rc_asset)
  rel_rc    <- rc_bucket / res_final$es

  # Warn if optimiser did not converge or result is far from target
  max_dev <- max(abs(rel_rc[!zero_b] - target_b[!zero_b]))
  if (max_dev > 0.02) {
    warning(sprintf("Largest RC deviation from target: %.1f pct. ",
                    max_dev * 100),
            "Long-only constraint may prevent exact risk parity.")
  }

  list(
    w_buckets  = round(w_buckets, 6),
    w_assets   = round(w_assets,  6),
    es         = res_final$es,
    rc_asset   = res_final$rc_asset,
    rc_bucket  = rc_bucket,
    rel_rc     = rel_rc,
    target_b   = target_b,
    max_rc_dev = max_dev,
    obj_val    = best_val,
    converged  = best_res$convergence == 0
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Run Dynamic Optimiser (3 models x 2 portfolios)
# ══════════════════════════════════════════════════════════════════════════════

scenario_files <- c(
  gaussian = "scenarios_gaussian.rds",
  fhs      = "scenarios_fhs.rds",
  evt      = "scenarios_evt.rds"
)

results_dynamic <- list()

for (model_name in names(scenario_files)) {
  cat(sprintf("\n── Model %s ──────────────────────────────────────\n",
              toupper(model_name)))
  scenarios <- readRDS(file.path(PATH_DATA, scenario_files[model_name]))

  for (pf_name in names(BUDGETS)) {
    target_b <- BUDGETS[[pf_name]]
    key      <- paste0(model_name, "_", pf_name)

    cat(sprintf("  Optimising %s | Target: %s\n", key,
                paste(sprintf("%s=%.0f%%", names(target_b), target_b*100),
                      collapse = ", ")))

    res <- optimise_risk_budget(scenarios, target_b)
    results_dynamic[[key]] <- res

    cat(sprintf("    Weights (buckets): equity=%.1f%%  bonds=%.1f%%  gold=%.1f%%\n",
                res$w_buckets["equity"]*100,
                res$w_buckets["bonds"]*100,
                res$w_buckets["gold"]*100))
    cat(sprintf("    Achieved RC:       equity=%.1f%%  bonds=%.1f%%  gold=%.1f%%\n",
                res$rel_rc["equity"]*100,
                res$rel_rc["bonds"]*100,
                res$rel_rc["gold"]*100))
    cat(sprintf("    Daily ES (95%%): %.4f  |  Max RC deviation: %.2f pct\n",
                res$es, res$max_rc_dev*100))
    if (!res$converged) cat("    WARNING: optimiser did not converge.\n")
  }
}

saveRDS(results_dynamic, file.path(PATH_DATA, "results_dynamic.rds"))
cat("\nDynamic results saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Static Strategy (benchmark)
# ══════════════════════════════════════════════════════════════════════════════
# Static weights = target budget percentages applied directly to bucket levels.
# No model or optimisation involved: just hold the target percentage in each
# asset-class bucket, distributed within buckets by fixed intra-bucket ratios.

cat("\n── Static Strategy ──────────────────────────────────────────────────────\n")

results_static <- list()

for (pf_name in names(BUDGETS)) {
  target_b <- BUDGETS[[pf_name]]

  # Static weights are simply the budget targets
  w_static_assets <- bucket_to_asset_weights(target_b)   # m-vector of asset weights

  cat(sprintf("\n  Portfolio %s | Static weights:\n", pf_name))
  print(round(w_static_assets * 100, 2))

  # Evaluate ES and risk contributions for each model's scenario set
  static_by_model <- list()
  for (model_name in names(scenario_files)) {
    scenarios <- readRDS(file.path(PATH_DATA, scenario_files[model_name]))
    res_s     <- compute_es_contributions(scenarios, w_static_assets)
    rc_b      <- aggregate_to_buckets(res_s$rc_asset)
    rel_rc_s  <- rc_b / res_s$es

    static_by_model[[model_name]] <- list(
      w_buckets = target_b,      # static: bucket weights == target budgets
      w_assets  = w_static_assets,
      es        = res_s$es,
      rc_asset  = res_s$rc_asset,
      rc_bucket = rc_b,
      rel_rc    = rel_rc_s,
      target_b  = target_b
    )

    cat(sprintf("    [%s]  ES=%.4f  RC: equity=%.1f%%  bonds=%.1f%%  gold=%.1f%%\n",
                model_name, res_s$es,
                rel_rc_s["equity"]*100,
                rel_rc_s["bonds"]*100,
                rel_rc_s["gold"]*100))
  }
  results_static[[pf_name]] <- list(
    w_assets     = w_static_assets,
    by_model     = static_by_model
  )
}

saveRDS(results_static, file.path(PATH_DATA, "results_static.rds"))
cat("\nStatic results saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Historical Backtest
# ══════════════════════════════════════════════════════════════════════════════
# Applies both dynamic and static weights over the full historical sample.
# Static strategy: annual calendar rebalancing to fixed budget weights.
# Dynamic strategy: rebalance only when realised ES risk-contribution shares
# (rolling window) breach target tolerance bands.
#
# NOTE: This is an IN-SAMPLE backtest since the factor SV model was fitted on
# the full sample. It evaluates whether the risk-contribution targets were
# approximately satisfied ex-post in terms of REALISED risk contributions.

# Convert log returns to arithmetic returns for correct portfolio compounding.
# If r_log = log(P_t / P_{t-1}), then simple return is exp(r_log) - 1.
returns_simple <- exp(as.matrix(log_returns_xts)) - 1   # T x m, raw simple returns

rc_band_breached <- function(window_returns, w_assets, target_b,
                             alpha      = ES_LEVEL,
                             band_lower = RC_BAND_LOWER,   # from 00_config.R
                             band_upper = RC_BAND_UPPER) { # from 00_config.R
  # -------------------------------------------------------------------------
  # Checks whether REALISED ES risk-contribution shares have drifted outside
  # the target tolerance bands using a rolling historical return window.
  #
  # DESIGN NOTE — historical returns, not model scenarios:
  #   The optimisation in Section 2 uses MODEL-BASED posterior-predictive
  #   scenarios (from 03A/B/C) to estimate forward-looking ES and set weights.
  #   This trigger check uses REALISED returns in a rolling lookback window.
  #   Rationale: (a) computationally feasible at daily frequency without
  #   re-running MCMC; (b) measures whether the portfolio's ACTUAL risk
  #   concentration has drifted, regardless of model assumptions.
  #   The two ES estimates will differ; the trigger is intentionally backward-
  #   looking while the optimisation is forward-looking.
  # -------------------------------------------------------------------------
  es_obj <- compute_es_contributions(window_returns, w_assets, alpha)
  if (!is.finite(es_obj$es) || es_obj$es <= 0) return(FALSE)
  rel_rc <- aggregate_to_buckets(es_obj$rc_asset) / es_obj$es
  active <- target_b > 0
  lower  <- target_b - band_lower
  upper  <- target_b + band_upper
  any(rel_rc[active] < lower[active] | rel_rc[active] > upper[active])
}

compute_backtest_static <- function(w_target, returns_simple, rebal_freq = "annual") {
  # Benchmark: annual calendar rebalancing to fixed target weights.
  T_bt     <- nrow(returns_simple)
  dates_bt <- as.Date(rownames(returns_simple))
  years_bt <- year(dates_bt)

  w_current      <- w_target
  pf_returns     <- numeric(T_bt)
  rebalance_flag <- logical(T_bt)
  weights_ts     <- matrix(NA, nrow = T_bt, ncol = ncol(returns_simple),
                           dimnames = list(NULL, TICKERS))

  for (t in seq_len(T_bt)) {
    if (rebal_freq == "annual" && t > 1 && years_bt[t] > years_bt[t - 1]) {
      w_current <- w_target
      rebalance_flag[t] <- TRUE
    }

    weights_ts[t, ] <- w_current
    pf_returns[t]   <- sum(w_current * returns_simple[t, ])

    # Buy-and-hold weight drift after realised arithmetic returns.
    w_current <- w_current * (1 + returns_simple[t, ])
    w_current <- w_current / sum(w_current)
  }

  data.frame(
    date      = dates_bt,
    pf_return = pf_returns,
    rebalance = rebalance_flag,
    weights_ts,
    check.names = FALSE
  )
}

compute_backtest_dynamic <- function(w_target, target_b, returns_simple,
                                     lookback   = RC_LOOKBACK_DAYS,    # from 00_config.R
                                     min_gap    = RC_MIN_DAYS_BETWEEN, # from 00_config.R
                                     band_lower = RC_BAND_LOWER,       # from 00_config.R
                                     band_upper = RC_BAND_UPPER) {     # from 00_config.R
  # Dynamic: trigger rebalance only when realised ES RC shares breach bands.
  # lookback  : rolling window length (days) for realised ES RC estimation
  # min_gap   : minimum trading days between consecutive rebalances (cool-down)
  # band_lower/upper: tolerance around target RC shares (e.g. 0.03 = ±3 pp)
  T_bt     <- nrow(returns_simple)
  dates_bt <- as.Date(rownames(returns_simple))

  w_current      <- w_target
  pf_returns     <- numeric(T_bt)
  rebalance_flag <- logical(T_bt)
  weights_ts     <- matrix(NA, nrow = T_bt, ncol = ncol(returns_simple),
                           dimnames = list(NULL, TICKERS))
  last_rebal_t   <- 1L

  for (t in seq_len(T_bt)) {
    can_check <- t > lookback
    cooled    <- (t - last_rebal_t) >= min_gap

    if (can_check && cooled) {
      win_idx <- (t - lookback):(t - 1)          # rolling window of realized returns
      win_ret <- returns_simple[win_idx, , drop = FALSE]
      if (rc_band_breached(win_ret, w_current, target_b,
                           band_lower = band_lower,
                           band_upper = band_upper)) {
        w_current <- w_target
        rebalance_flag[t] <- TRUE
        last_rebal_t <- t
      }
    }

    weights_ts[t, ] <- w_current
    pf_returns[t]   <- sum(w_current * returns_simple[t, ])

    # Buy-and-hold weight drift after realised arithmetic returns.
    w_current <- w_current * (1 + returns_simple[t, ])
    w_current <- w_current / sum(w_current)
  }

  data.frame(
    date      = dates_bt,
    pf_return = pf_returns,
    rebalance = rebalance_flag,
    weights_ts,
    check.names = FALSE
  )
}

# Run backtest for each strategy (2 portfolios x [3 dynamic + 1 static] = 8)
cat("\nRunning historical backtests...\n")
backtests <- list()

for (pf_name in names(BUDGETS)) {
  # Static strategy backtest
  w_static <- results_static[[pf_name]]$w_assets
  bt_key   <- paste0("static_", pf_name)
  backtests[[bt_key]] <- compute_backtest_static(w_static, returns_simple)
  cat(sprintf("  Backtest complete: %s\n", bt_key))

  # Dynamic strategy backtest (trigger-based ES RC band rebalancing)
  for (model_name in names(scenario_files)) {
    key      <- paste0(model_name, "_", pf_name)
    w_dyn    <- results_dynamic[[key]]$w_assets
    bt_key2  <- paste0("dynamic_", key)
    backtests[[bt_key2]] <- compute_backtest_dynamic(w_dyn, BUDGETS[[pf_name]], returns_simple)
    cat(sprintf("  Backtest complete: %s (rebalances=%d)\n",
                bt_key2, sum(backtests[[bt_key2]]$rebalance)))
  }
}

# Compute summary statistics for each backtest
backtest_summary <- do.call(rbind, lapply(names(backtests), function(key) {
  bt      <- backtests[[key]]
  pf_r    <- bt$pf_return
  ann_ret <- (prod(1 + pf_r)^(252 / length(pf_r))) - 1
  ann_vol <- sd(pf_r)   * sqrt(252)
  sharpe  <- ann_ret / ann_vol
  # Historical ES (empirical): conditional mean of losses in the alpha tail.
  var_cut <- as.numeric(quantile(pf_r, ES_LEVEL))
  hist_es <- -mean(pf_r[pf_r <= var_cut])
  # Max drawdown
  cum_r   <- cumprod(1 + pf_r)
  roll_max <- cummax(cum_r)
  max_dd  <- max((roll_max - cum_r) / roll_max)
  n_rebal <- if ("rebalance" %in% names(bt)) sum(bt$rebalance) else NA_integer_

  data.frame(
    strategy = key,
    ann_ret  = round(ann_ret  * 100, 2),
    ann_vol  = round(ann_vol  * 100, 2),
    sharpe   = round(sharpe,  3),
    hist_es  = round(hist_es  * 100, 4),   # percent
    max_dd   = round(max_dd   * 100, 2),
    n_rebal  = as.integer(n_rebal),
    stringsAsFactors = FALSE
  )
}))

cat("\nBacktest summary:\n")
print(backtest_summary)

saveRDS(backtests,         file.path(PATH_DATA, "backtests.rds"))
saveRDS(backtest_summary,  file.path(PATH_DATA, "backtest_summary.rds"))
cat("\nBacktest results saved.\n")
cat("04_risk_budgeting.R complete.\n")
