# ==============================================================================
# 05_results.R
# ------------------------------------------------------------------------------
# Loads all results and produces comparison tables and figures.
#
# Comparison axes:
#   1. MODEL:     Gaussian (A) vs FHS (B) vs EVT (C)
#   2. STRATEGY:  Dynamic (risk-contribution targeted) vs Static (budget weights)
#   3. PORTFOLIO: pf1 (45/45/10) vs pf2 (50/50/0)
#
# Figures (saved to PATH_CHART):
#   01_factor_logvol_paths.png      -- Factor log-variance paths (from model fit)
#   02_loadings_heatmap.png         -- Already saved in 02_model_fit.R
#   03_qq_std_residuals.png         -- Already saved in 03B_scenarios_fhs.R
#   04_es_model_comparison.png      -- ES estimates: Gaussian vs FHS vs EVT
#   05_rc_bars_pf1.png              -- Risk contribution bars, portfolio 1
#   06_rc_bars_pf2.png              -- Risk contribution bars, portfolio 2
#   07_weights_comparison_pf1.png   -- Optimal weights: static vs dynamic (pf1)
#   08_weights_comparison_pf2.png   -- Optimal weights: static vs dynamic (pf2)
#   09_cumulative_returns_pf1.png   -- Cumulative return backtest, pf1
#   10_cumulative_returns_pf2.png   -- Cumulative return backtest, pf2
#
# Tables (printed to console):
#   T1: Portfolio weights (static vs dynamic, by model)
#   T2: ES estimates and risk contributions at optimal weights
#   T3: Backtest performance summary
# ==============================================================================

source(file.path("C:/Users/Tristan Leiter/Documents/MV-Factor-Stochastic-Volatility-Models",
                 "01_Code/Pipeline/00_config.R"))

results_dynamic  <- readRDS(file.path(PATH_DATA, "results_dynamic.rds"))
results_static   <- readRDS(file.path(PATH_DATA, "results_static.rds"))
backtests        <- readRDS(file.path(PATH_DATA, "backtests.rds"))
backtest_summary <- readRDS(file.path(PATH_DATA, "backtest_summary.rds"))
log_returns_xts  <- readRDS(file.path(PATH_DATA, "log_returns.rds"))

model_names <- c("gaussian", "fhs", "evt")
model_labels <- c(gaussian = "Gaussian (A)", fhs = "FHS (B)", evt = "EVT (C)")
bucket_names <- c("equity", "bonds", "gold")
bucket_cols  <- c(equity = "#E07B54", bonds = "#5B8DB8", gold = "#D4AF37")


# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1: Portfolio Weights Comparison
# ══════════════════════════════════════════════════════════════════════════════
cat("\n═══ TABLE 1: Portfolio Weights (%) ═══════════════════════════════════\n")

for (pf_name in names(BUDGETS)) {
  cat(sprintf("\nPortfolio %s | Target budgets: %s\n",
              toupper(pf_name),
              paste(sprintf("%s=%.0f%%", names(BUDGETS[[pf_name]]),
                            BUDGETS[[pf_name]]*100), collapse = ", ")))

  # Build a row per strategy variant
  wt_table <- data.frame(
    strategy = character(),
    SPY = numeric(), QQQ = numeric(), TLT = numeric(),
    SHY = numeric(), GLD = numeric(),
    stringsAsFactors = FALSE
  )

  # Static row
  w_s <- results_static[[pf_name]]$w_assets * 100
  wt_table <- rbind(wt_table,
    data.frame(strategy = "Static", t(w_s), stringsAsFactors = FALSE))

  # Dynamic rows (one per model)
  for (mn in model_names) {
    key <- paste0(mn, "_", pf_name)
    w_d <- results_dynamic[[key]]$w_assets * 100
    wt_table <- rbind(wt_table,
      data.frame(strategy = model_labels[mn], t(w_d), stringsAsFactors = FALSE))
  }

  print(wt_table, row.names = FALSE, digits = 3)
}


# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2: ES and Risk Contributions at Optimal Weights
# ══════════════════════════════════════════════════════════════════════════════
cat("\n═══ TABLE 2: ES and Realised Risk Contributions ═══════════════════════\n")

for (pf_name in names(BUDGETS)) {
  cat(sprintf("\nPortfolio %s\n", toupper(pf_name)))
  target_b <- BUDGETS[[pf_name]]

  rc_table <- data.frame(
    strategy  = character(),
    model     = character(),
    es_pct    = numeric(),
    rc_equity = numeric(),
    rc_bonds  = numeric(),
    rc_gold   = numeric(),
    max_dev   = numeric(),
    stringsAsFactors = FALSE
  )

  # Static rows
  for (mn in model_names) {
    s_res   <- results_static[[pf_name]]$by_model[[mn]]
    rel_rc  <- s_res$rel_rc * 100
    max_dev <- max(abs(s_res$rel_rc[target_b > 0] - target_b[target_b > 0])) * 100
    rc_table <- rbind(rc_table, data.frame(
      strategy  = "Static",
      model     = model_labels[mn],
      es_pct    = round(s_res$es * 100, 4),
      rc_equity = round(rel_rc["equity"], 1),
      rc_bonds  = round(rel_rc["bonds"],  1),
      rc_gold   = round(rel_rc["gold"],   1),
      max_dev   = round(max_dev, 2),
      stringsAsFactors = FALSE
    ))
  }

  # Dynamic rows
  for (mn in model_names) {
    key    <- paste0(mn, "_", pf_name)
    d_res  <- results_dynamic[[key]]
    rel_rc <- d_res$rel_rc * 100
    rc_table <- rbind(rc_table, data.frame(
      strategy  = "Dynamic",
      model     = model_labels[mn],
      es_pct    = round(d_res$es * 100, 4),
      rc_equity = round(rel_rc["equity"], 1),
      rc_bonds  = round(rel_rc["bonds"],  1),
      rc_gold   = round(rel_rc["gold"],   1),
      max_dev   = round(d_res$max_rc_dev * 100, 2),
      stringsAsFactors = FALSE
    ))
  }

  colnames(rc_table)[4:6] <- paste0("RC ", c("Eq%", "Bd%", "Au%"))
  colnames(rc_table)[7]   <- "MaxDev%"
  print(rc_table, row.names = FALSE)
  cat("  Target RC:  equity=", BUDGETS[[pf_name]]["equity"]*100,
      "  bonds=",  BUDGETS[[pf_name]]["bonds"]*100,
      "  gold=",   BUDGETS[[pf_name]]["gold"]*100, "\n")
  cat("  ES is daily 95% Expected Shortfall (positive = expected loss in worst 5%).\n")
  cat("  MaxDev = max absolute deviation of realised RC from target budget.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3: Backtest Performance
# ══════════════════════════════════════════════════════════════════════════════
cat("\n═══ TABLE 3: Backtest Performance Summary ══════════════════════════════\n")
print(backtest_summary, row.names = FALSE)
cat("  ann_ret/ann_vol in %;  sharpe = ann_ret / ann_vol;\n",
    " hist_es = empirical 95% ES (daily %);  max_dd = maximum drawdown (%);\n",
    " n_rebal = number of rebalancing events over the full sample.\n",
    " Static: annual calendar rebalancing (~17 events).\n",
    " Dynamic: trigger-based rebalancing when rolling-window ES RC breaches ±",
    RC_BAND_LOWER * 100, "pp bands (min gap:", RC_MIN_DAYS_BETWEEN, "days).\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: ES Model Comparison (same weights, different tail models)
# ══════════════════════════════════════════════════════════════════════════════
# For each portfolio, compare the ES estimate of the STATIC weights
# across the three models. This isolates the model effect on tail risk estimation.

es_compare <- do.call(rbind, lapply(names(BUDGETS), function(pf_name) {
  do.call(rbind, lapply(model_names, function(mn) {
    data.frame(
      portfolio = pf_name,
      model     = model_labels[mn],
      es_pct    = results_static[[pf_name]]$by_model[[mn]]$es * 100,
      stringsAsFactors = FALSE
    )
  }))
}))

p_es <- ggplot(es_compare, aes(x = model, y = es_pct, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f%%", es_pct)), vjust = -0.4, size = 3.5) +
  facet_wrap(~ portfolio, labeller = as_labeller(
    c(pf1 = "Portfolio 1 (45/45/10)", pf2 = "Portfolio 2 (50/50/0)"))) +
  scale_fill_manual(values = c("Gaussian (A)" = "#5B8DB8",
                               "FHS (B)"      = "#E07B54",
                               "EVT (C)"      = "#6B9E6B"),
                    name = "Tail model") +
  labs(title = "Daily ES (95%) by Tail Model",
       subtitle = "Static portfolio weights; heavier bar = more conservative tail estimate",
       x = NULL, y = "Daily ES (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(PATH_CHART, "04_es_model_comparison.png"),
       p_es, width = 9, height = 5, dpi = 150)
cat("\nFigure 4 saved: es_model_comparison.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURES 5-6: Risk Contribution Bar Charts (per portfolio)
# ══════════════════════════════════════════════════════════════════════════════

make_rc_bar_chart <- function(pf_name, results_dynamic, results_static, model_names) {
  target_b <- BUDGETS[[pf_name]]

  rc_rows <- list()

  # Static entries
  for (mn in model_names) {
    rel_rc <- results_static[[pf_name]]$by_model[[mn]]$rel_rc * 100
    for (b in bucket_names) {
      rc_rows[[length(rc_rows)+1]] <- data.frame(
        strategy = "Static", model = model_labels[mn],
        bucket = b, rc_pct = rel_rc[b], stringsAsFactors = FALSE)
    }
  }
  # Dynamic entries
  for (mn in model_names) {
    rel_rc <- results_dynamic[[paste0(mn, "_", pf_name)]]$rel_rc * 100
    for (b in bucket_names) {
      rc_rows[[length(rc_rows)+1]] <- data.frame(
        strategy = "Dynamic", model = model_labels[mn],
        bucket = b, rc_pct = rel_rc[b], stringsAsFactors = FALSE)
    }
  }

  rc_df <- do.call(rbind, rc_rows)
  rc_df$bucket   <- factor(rc_df$bucket, levels = bucket_names)
  rc_df$x_label  <- paste0(rc_df$strategy, "\n", rc_df$model)
  rc_df$x_label  <- factor(rc_df$x_label,
                            levels = unique(rc_df$x_label))

  # Target lines
  target_lines <- data.frame(
    bucket  = names(target_b),
    target  = as.numeric(target_b) * 100
  )
  target_lines <- target_lines[target_lines$target > 0, ]
  target_lines$bucket <- factor(target_lines$bucket, levels = bucket_names)

  ggplot(rc_df, aes(x = x_label, y = rc_pct, fill = bucket)) +
    geom_col(position = "stack", width = 0.7) +
    geom_hline(data = target_lines,
               aes(yintercept = target, colour = bucket),
               linetype = "dashed", linewidth = 0.8) +
    scale_fill_manual(values  = bucket_cols, name = "Bucket") +
    scale_colour_manual(values = bucket_cols, guide = "none") +
    labs(
      title    = sprintf("Risk Contributions by Bucket — Portfolio %s", toupper(pf_name)),
      subtitle = sprintf("Target: %s  |  Dashed lines = target levels",
                         paste(sprintf("%s=%.0f%%", names(target_b), target_b*100),
                               collapse = ", ")),
      x = NULL, y = "Risk Contribution (% of total ES)"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 8),
          legend.position = "right")
}

for (pf_name in names(BUDGETS)) {
  p_rc <- make_rc_bar_chart(pf_name, results_dynamic, results_static, model_names)
  fname <- sprintf("0%d_rc_bars_%s.png",
                   4 + which(names(BUDGETS) == pf_name), pf_name)
  ggsave(file.path(PATH_CHART, fname), p_rc, width = 11, height = 6, dpi = 150)
  cat(sprintf("Figure saved: %s\n", fname))
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURES 7-8: Asset Weight Comparison (static vs dynamic, per portfolio)
# ══════════════════════════════════════════════════════════════════════════════

etf_cols <- c(SPY = "#1f77b4", QQQ = "#aec7e8",
              TLT = "#ff7f0e", SHY = "#ffbb78", GLD = "#D4AF37")

make_weight_chart <- function(pf_name) {
  rows <- list()

  # Static weights
  w_s <- results_static[[pf_name]]$w_assets * 100
  for (tk in TICKERS)
    rows[[length(rows)+1]] <- data.frame(
      strategy = "Static", model = "—", ticker = tk, weight = w_s[tk],
      stringsAsFactors = FALSE)

  # Dynamic weights
  for (mn in model_names) {
    w_d <- results_dynamic[[paste0(mn, "_", pf_name)]]$w_assets * 100
    for (tk in TICKERS)
      rows[[length(rows)+1]] <- data.frame(
        strategy = "Dynamic", model = model_labels[mn],
        ticker = tk, weight = w_d[tk], stringsAsFactors = FALSE)
  }

  wt_df <- do.call(rbind, rows)
  wt_df$x_label <- ifelse(wt_df$strategy == "Static", "Static",
                           paste0("Dynamic\n", wt_df$model))
  wt_df$x_label <- factor(wt_df$x_label, levels = unique(wt_df$x_label))
  wt_df$ticker  <- factor(wt_df$ticker, levels = TICKERS)

  ggplot(wt_df, aes(x = x_label, y = weight, fill = ticker)) +
    geom_col(position = "stack", width = 0.65) +
    geom_text(aes(label = ifelse(weight > 2, sprintf("%.1f%%", weight), "")),
              position = position_stack(vjust = 0.5), size = 3, colour = "white") +
    scale_fill_manual(values = etf_cols, name = "ETF") +
    labs(
      title    = sprintf("Portfolio Weights — Portfolio %s", toupper(pf_name)),
      subtitle = "Static = target budget percentages; Dynamic = risk-budget optimised",
      x = NULL, y = "Weight (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}

for (pf_name in names(BUDGETS)) {
  p_wt  <- make_weight_chart(pf_name)
  fname <- sprintf("0%d_weights_comparison_%s.png",
                   6 + which(names(BUDGETS) == pf_name), pf_name)
  ggsave(file.path(PATH_CHART, fname), p_wt, width = 9, height = 5, dpi = 150)
  cat(sprintf("Figure saved: %s\n", fname))
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURES 9-10: Cumulative Returns Backtest
# ══════════════════════════════════════════════════════════════════════════════

make_cumret_chart <- function(pf_name) {
  # Pull static and all dynamic backtest series for this portfolio
  keys_plot <- c(
    paste0("static_", pf_name),
    paste0("dynamic_gaussian_",  pf_name),
    paste0("dynamic_fhs_",       pf_name),
    paste0("dynamic_evt_",       pf_name)
  )
  labels_plot <- c(
    "Static",
    "Dynamic — Gaussian (A)",
    "Dynamic — FHS (B)",
    "Dynamic — EVT (C)"
  )
  colours_plot <- c("grey40", "#5B8DB8", "#E07B54", "#6B9E6B")
  names(colours_plot) <- labels_plot

  # Combine into long-format data frame
  cum_rows <- mapply(function(key, lab) {
    bt <- backtests[[key]]
    data.frame(
      date     = bt$date,
      cum_ret  = cumprod(1 + bt$pf_return) - 1,
      strategy = lab,
      stringsAsFactors = FALSE
    )
  }, keys_plot, labels_plot, SIMPLIFY = FALSE)

  cum_df <- do.call(rbind, cum_rows)
  cum_df$strategy <- factor(cum_df$strategy, levels = labels_plot)

  # Rebalancing event markers: vertical ticks at dates when the dynamic
  # strategy triggered a rebalance (one set of markers per model variant).
  # Drawn at the bottom of the chart to avoid obscuring the return lines.
  rebal_rows <- do.call(rbind, lapply(seq_along(keys_plot), function(i) {
    key <- keys_plot[i]
    bt  <- backtests[[key]]
    if (!"rebalance" %in% names(bt)) return(NULL)
    rdates <- bt$date[bt$rebalance]
    if (length(rdates) == 0) return(NULL)
    data.frame(date = rdates, strategy = labels_plot[i],
               stringsAsFactors = FALSE)
  }))

  p <- ggplot(cum_df, aes(x = date, y = cum_ret * 100, colour = strategy)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = colours_plot, name = "Strategy") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("Cumulative Returns — Portfolio %s", toupper(pf_name)),
      subtitle = "Static (annual rebalance) vs Dynamic (rebalance on ES RC band breaches)",
      x = NULL, y = "Cumulative Return (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))

  # Overlay rebalancing rug marks on the x-axis if data is available
  if (!is.null(rebal_rows) && nrow(rebal_rows) > 0) {
    rebal_rows$strategy <- factor(rebal_rows$strategy, levels = labels_plot)
    p <- p +
      geom_rug(data = rebal_rows,
               aes(x = date, colour = strategy),
               sides = "b",          # bottom axis only
               linewidth = 0.4,
               alpha = 0.6,
               inherit.aes = FALSE)
  }
  p
}

for (pf_name in names(BUDGETS)) {
  p_cum  <- make_cumret_chart(pf_name)
  fname  <- sprintf("%02d_cumulative_returns_%s.png",
                    8 + which(names(BUDGETS) == pf_name), pf_name)
  ggsave(file.path(PATH_CHART, fname), p_cum, width = 11, height = 5.5, dpi = 150)
  cat(sprintf("Figure saved: %s\n", fname))
}

cat("\n05_results.R complete. All tables printed and figures saved to:\n")
cat("  ", PATH_CHART, "\n")
