# Goodness-of-fit diagnostics

#' Goodness-of-fit summary
#'
#' @param fit A fitted \code{pmrm_nlme} object.
#' @param bins Number of bins for residual quantile grouping.
#' @return List with population residuals, conditional (EB) residuals,
#'   RMSE, and information criteria.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' g <- goodness_of_fit(fit)
#' g$rmse
#' g$AIC
#' @export
goodness_of_fit <- function(fit, bins = 10) {
  pres <- stats::residuals(fit, type = "response", level = 0)
  cres <- stats::residuals(fit, type = "normalized")
  list(
    population_residuals = as.numeric(pres),
    conditional_residuals = as.numeric(cres),
    rmse = sqrt(mean(as.numeric(cres)^2, na.rm = TRUE)),
    AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    BIC = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
    residual_bins = cut(as.numeric(cres), breaks = bins)
  )
}

#' Visual predictive check (VPC)
#'
#' @param fit Fitted \code{pmrm_nlme} object.
#' @param n_sims Number of simulated datasets for the VPC.
#' @param probs Quantiles to report.
#' @param seed Random seed.
#' @return Data frame of observed vs simulated quantiles by time bin.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20,
#'                                   times = seq(0, 24, by = 6), seed = 2)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' vpc(fit, n_sims = 5, seed = 1)
#' @export
vpc <- function(fit, n_sims = 50, probs = c(0.05, 0.5, 0.95), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat <- tryCatch(nlme::getData(fit), error = function(e) NULL)
  if (is.null(dat)) dat <- fit$data
  dat <- as.data.frame(dat)
  yhat_pop <- stats::fitted(fit, level = 0)
  sd_res <- stats::sd(stats::residuals(fit, level = 0), na.rm = TRUE)

  sim_q <- t(vapply(seq_len(n_sims), function(s) {
    ysim <- rnorm(nrow(dat), yhat_pop, sd_res)
    vapply(probs, function(p) stats::quantile(ysim, p, na.rm = TRUE),
           numeric(1))
  }, numeric(length(probs))))

  resp <- all.vars(stats::formula(fit))[1]
  if (!resp %in% names(dat)) resp <- "y"
  obs_q <- vapply(probs, function(p)
    stats::quantile(dat[[resp]], p, na.rm = TRUE), numeric(1))

  data.frame(
    prob = paste0("q", probs),
    observed = obs_q,
    sim_median = apply(sim_q, 2, median, na.rm = TRUE),
    sim_lo = apply(sim_q, 2, function(x) stats::quantile(x, 0.05, na.rm = TRUE)),
    sim_hi = apply(sim_q, 2, function(x) stats::quantile(x, 0.95, na.rm = TRUE))
  )
}
