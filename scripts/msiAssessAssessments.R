## Generic interval-change assessment for any index-like series
## draws: iterations x T matrix of posterior draws for X_t
## years: length-T vector of time labels (e.g. years)
## year1, year2: the two time points to compare (must be in 'years')
## scale: "log" if draws are on log scale; "index" if on index scale
## prob: credible level (default 0.95)
## estimand_label: optional textual description of X_t
assess_interval_change_from_draws <- function(draws,
                                            years,
                                            year1,
                                            year2,
                                            scale = c("log", "index"),
                                            prob = 0.95,
                                            estimand_label = NULL) {
  scale <- match.arg(scale)
  
  if (length(dim(draws)) != 2L) {
    stop("`draws` must be a 2D matrix: iterations x time.")
  }
  n_iter <- nrow(draws)
  T <- ncol(draws)
  
  if (length(years) != T) {
    stop("length(years) must equal ncol(draws).")
  }
  
  id1 <- match(year1, years)
  id2 <- match(year2, years)
  if (any(is.na(c(id1, id2)))) {
    stop("Both `year1` and `year2` must appear in `years`.")
  }
  if (id1 == id2) stop("year1 and year2 must be different.")
  
  # keep the (year1, year2) ordering as given by the user.
  if (scale == "log") {
    ## draws are X_t on log scale; difference is log-ratio
    log_diff <- draws[, id2] - draws[, id1]
    ratio    <- exp(log_diff)
  } else {
    ## draws are X_t on index scale; work via ratio then log
    ratio    <- draws[, id2] / draws[, id1]
    log_diff <- log(ratio)
  }
  
  alpha <- 1 - prob
  q_ratio    <- stats::quantile(ratio,    probs = c(alpha/2, 0.5, 1 - alpha/2), na.rm = TRUE)
  q_log_diff <- stats::quantile(log_diff, probs = c(alpha/2, 0.5, 1 - alpha/2), na.rm = TRUE)
  
  p_gt1 <- mean(ratio > 1, na.rm = TRUE)
  p_lt1 <- mean(ratio < 1, na.rm = TRUE)
  
  if (is.null(estimand_label)) {
    if (scale == "index") {
      estimand_label <- paste(
        "Ratio of the index at year2 to year1,",
        "X(year2) / X(year1), on the original index scale."
      )
    } else {
      estimand_label <- paste(
        "Difference in log-index between year2 and year1,",
        "log X(year2) - log X(year1), i.e. the log of the index ratio."
      )
    }
  }
  
  out <- list(
    year1 = year1,
    year2 = year2,
    estimand = estimand_label,
    ratio = list(
      lower  = unname(q_ratio[1]),
      median = unname(q_ratio[2]),
      upper  = unname(q_ratio[3]),
      p_gt1  = p_gt1,
      p_lt1  = p_lt1
    ),
    log_diff = list(
      lower  = unname(q_log_diff[1]),
      median = unname(q_log_diff[2]),
      upper  = unname(q_log_diff[3]),
      p_gt0  = mean(log_diff > 0, na.rm = TRUE),
      p_lt0  = mean(log_diff < 0, na.rm = TRUE)
    )
  )
  class(out) <- "msi_interval_assess"
  out
}

## Optional print method for nice console output
print.msi_interval_assess <- function(x, digits = 3, ...) {
  cat("Interval assessment for change between", x$year1, "and", x$year2, "\n")
  cat("Estimand:\n  ", x$estimand, "\n\n")
  
  r <- x$ratio
  cat("Ratio X(year2)/X(year1):\n")
  cat("  median =", round(r$median, digits),
      " (", round(r$lower, digits), ", ", round(r$upper, digits), ")",
      "  [", format(signif(100 * (r$median - 1), 3)), "% change]\n", sep = "")
  cat("  P(ratio > 1) =", round(r$p_gt1, digits),
      "; P(ratio < 1) =", round(r$p_lt1, digits), "\n\n")
  
  d <- x$log_diff
  cat("Log-difference log X(year2) - log X(year1):\n")
  cat("  median =", round(d$median, digits),
      " (", round(d$lower, digits), ", ", round(d$upper, digits), ")\n", sep = "")
  cat("  P(diff > 0) =", round(d$p_gt0, digits),
      "; P(diff < 0) =", round(d$p_lt0, digits), "\n")
  
  invisible(x)
}

## Convenience wrapper: interval assessment for Freeman M_prime
## out: object returned by run_full_analysis(...), with a Freeman fit
## year1, year2: years to compare (must match out$sim$years)
assess_freeman_Mprime_change <- function(out, year1, year2, prob = 0.95) {
  if (is.null(out$fits$freeman)) {
    stop("No Freeman fit found in `out$fits$freeman`.")
  }
  fit <- out$fits$freeman
  sims <- fit$BUGSoutput$sims.list
  
  if (is.null(sims$Mprime)) {
    stop("Freeman fit does not contain `Mprime` in sims.list.")
  }
  
  draws <- sims$Mprime # iterations x T (log mean index across spp)
  years <- out$sim$years # time labels
  
  estimand_label <- paste(
    "Change in the mean log-index M_prime across species between year1 and year2,",
    "equivalently the ratio of the geometric-mean abundance index",
    "MSI(year2) / MSI(year1) implied by the Freeman model."
  )
  
  assess_interval_change_from_draws(
    draws = draws,
    years = years,
    year1 = year1,
    year2 = year2,
    scale = "log",
    prob = prob,
    estimand_label = estimand_label
  )
}

# ## Suppose you have run:
# ## out <- run_full_analysis(..., fit_models = c("freeman", ...))
# 
# # Assess net change from baseline year to last year
# res1 <- assess_freeman_Mprime_change(out, year1 = min(out$sim$years), year2 = max(out$sim$years))
# print(res1)
# 
# # Assess net change over a shorter interval, say 2000–2010
# res2 <- assess_freeman_Mprime_change(out, year1 = 2000, year2 = 2010)
# print(res2)

## Partial JAGS model (two options, M_smooth and M_full)
# Structural trend-only:
draws_trend <- out$fits$partial$BUGSoutput$sims.list$M_smooth  # log scale
assess_interval_change_from_draws(draws_trend, out$sim$years, 2, 28,
                                scale = "log",
                                estimand_label = "Change in common log-index M_smooth (trend only)"
)

# Trend + shocks:
draws_full <- out$fits$partial$BUGSoutput$sims.list$M_full
assess_interval_change_from_draws(draws_full, out$sim$years, 2, 28,
                                scale = "log",
                                estimand_label = "Change in common log-index M_full (trend + common shocks)"
)
# 
# ## Bayesian geomean MSI (Soldaat-type)
# draws_msi <- out$fits$bayes_geomean$BUGSoutput$sims.list$MSI  # index scale
# assess_interval_change_from_draws(draws_msi, out$sim$years, 1990, 2020,
#                                 scale = "index",
#                                 estimand_label = "Ratio of Bayesian geomean MSI between the two years"
# )
# 
# ## Post-smoothed MSI or M'
# # (schematic example)
# # Suppose 'Mprime_sm_draws' is iterations x T of smoothed M' (log scale)
freeman_Mprime_smooth <- posthoc_smooth_Mprime(
  fit = out$fits$freeman,
  years = out$sim$years,
  basis = "ruppert",
  num_knots = 12
)
out$results$freeman_Mprime_smooth <- freeman_Mprime_smooth
MSI_list <- list(
  Freeman_raw = out$results$freeman$MSI,
  Freeman_Mp_S = freeman_Mprime_smooth$MSI
)
p_MSI <- plot_indicator(out$sim$years, MSI_list, truths = NULL,
                        title = "Freeman MSI: raw vs post-smoothed M'")
print(p_MSI)
res_sm <- assess_interval_change_from_draws(
  draws = freeman_Mprime_smooth$MSI_sm_draws,
  years = out$sim$years,
  year1 = 1,
  year2 = 30,
  scale = "log",
  estimand_label = "Change in post-smoothed mean log-index M_prime"
)
print(res_sm)
