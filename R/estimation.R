# Inference wrappers: unified estimation entry point and parametric bootstrap

#' Unified estimation wrapper
#'
#' @param model Model family: "nlme", "delayed", or "joint".
#' @param ... Arguments forwarded to the corresponding fitting function.
#' @return The fitted model object.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- fit_pmrm("nlme", y ~ time * trt, data = sim)
#' class(fit)
#' @export
fit_pmrm <- function(model = c("nlme", "delayed", "joint"), ...) {
  switch(match.arg(model),
         "nlme"    = pmrm_nlme(...),
         "delayed" = pmrm_delayed(...),
         "joint"   = pmrm_joint(...)
  )
}

#' Parametric bootstrap for confidence intervals
#'
#' Refits the supplied model on \code{n_boot} datasets simulated from the
#' fitted model's estimated parameters (resampling residuals conditional on
#' the fixed-effects mean structure).
#'
#' @param fit A \code{pmrm_nlme} object.
#' @param statistic Function of the refit returning the scalar of interest.
#' @param n_boot Number of bootstrap replicates.
#' @param seed Random seed.
#' @param level Confidence level.
#' @return List with \code{estimate}, \code{ci}, and the bootstrap draws.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' boot <- parametric_bootstrap(fit, n_boot = 5, seed = 1)
#' boot$estimate
#' boot$ci
#' @export
parametric_bootstrap <- function(fit, statistic = function(f) stats::coef(f)[1],
                                 n_boot = 200, seed = NULL, level = 0.95) {
  if (!is.null(seed)) set.seed(seed)
  dat <- tryCatch(nlme::getData(fit), error = function(e) NULL)
  if (is.null(dat)) dat <- fit$data
  dat <- as.data.frame(dat)
  yhat <- stats::fitted(fit, level = 0)
  resid_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)

  boots <- replicate(n_boot, {
    d <- dat
    d$y_boot <- rnorm(nrow(d), yhat[match(seq_len(nrow(d)), seq_len(nrow(d)))],
                      resid_sd)
    # overwrite response column used by the formula
    resp <- all.vars(stats::formula(fit))[1]
    if (!resp %in% names(d)) resp <- "y"
    d[[resp]] <- d$y_boot
    rf <- tryCatch(suppressWarnings(
      update(fit, data = d)), error = function(e) NULL)
    if (is.null(rf)) NA_real_ else tryCatch(statistic(rf),
                                            error = function(e) NA_real_)
  })

  est <- tryCatch(statistic(fit), error = function(e) NA_real_)
  list(estimate = est,
       ci = stats::quantile(boots, c((1 - level) / 2, 1 - (1 - level) / 2),
                            na.rm = TRUE),
       boots = boots)
}
