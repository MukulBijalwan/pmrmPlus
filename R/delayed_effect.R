# Delayed treatment effect models

#' Delayed treatment effect function
#'
#' Computes the fractional treatment effect at time \code{time} given
#' treatment initiation time \code{trt_time} and a delay \code{delay}.
#'
#' @param time Numeric vector of observation times.
#' @param trt_time Treatment initiation time (default 0).
#' @param delay Delay before onset begins (tau).
#' @param shape One of "exponential", "step", "linear".
#' @param lambda Onset rate for the exponential shape.
#' @return Numeric vector in [0, 1] giving the fraction of the full effect
#'   realized at each time.
#' @export
delayed_effect <- function(time, trt_time = 0, delay = 0,
                           shape = c("exponential", "step", "linear"),
                           lambda = 0.1) {
  shape <- match.arg(shape)
  t_eff <- pmax(0, time - trt_time - delay)
  switch(shape,
    "step"        = as.numeric(t_eff > 0),
    "linear"      = {
      denom <- max(t_eff, na.rm = TRUE)
      if (!is.finite(denom) || denom <= 0) rep(1, length(t_eff))
      else pmin(t_eff / denom, 1)
    },
    "exponential" = 1 - exp(-lambda * t_eff)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Fit a progression model with delayed treatment onset
#'
#' The delayed effect enters as \code{delta * Trt_i * d(t_ij)} where d() is a
#' fractional onset function. When \code{delay_estimation = "profile"}, the
#' delay is estimated by grid search over \code{delay_grid}, selecting the
#' value maximizing log-likelihood; an approximate confidence set is derived
#' from the chi-squared threshold (log-lik drop of qchisq(0.95, 1)/2).
#'
#' @param formula LME-style formula, e.g. \code{y ~ time * trt}.
#' @param data Long-format data frame.
#' @param delay_estimation One of "fixed", "profile", or "bayesian"
#'   ("bayesian" currently falls back to profile with a warning).
#' @param delay_grid Grid of candidate delays (in the units of \code{time}).
#' @param shape Onset shape passed to \code{delayed_effect()}.
#' @param lambda Onset rate for exponential onset.
#' @param progression Progression shape passed to \code{pmrm_nlme()}.
#' @param conf_level Confidence level for the delay confidence interval.
#' @param ... Additional arguments to \code{pmrm_nlme()}.
#' @return Object of class \code{"pmrm_delayed"} containing the best fit,
#'   the estimated delay, the profile likelihood, and the delay CI.
#' @export
pmrm_delayed <- function(formula, data,
                         delay_estimation = c("fixed", "profile", "bayesian"),
                         delay_grid = seq(0, 24, by = 3),
                         shape = c("exponential", "step", "linear"),
                         lambda = 0.3,
                         progression = "linear",
                         conf_level = 0.95,
                         ...) {
  delay_estimation <- match.arg(delay_estimation)
  shape <- match.arg(shape)

  if (delay_estimation == "bayesian") {
    warning("Bayesian delayed-effect fitting requires rstan; using profile likelihood.")
    delay_estimation <- "profile"
  }

  fit_one <- function(delay) {
    d <- as.data.frame(data)
    trt_col <- grep("^trt$", names(d), ignore.case = TRUE)[1]
    d$trt_eff <- as.numeric(d[[trt_col]]) *
      delayed_effect(d$time, delay = delay, shape = shape, lambda = lambda)
    f <- tryCatch(
      suppressWarnings(pmrm_nlme(update(formula, . ~ . + trt_eff), data = d,
                                 progression = progression, ...)),
      error = function(e) NULL
    )
    if (!inherits(f, "lme")) NULL else f
  }

  if (delay_estimation == "fixed") {
    fit <- fit_one(0)
    fit$delay_estimate <- 0
    return(structure(fit, class = c("pmrm_delayed", class(fit))))
  }

  fits <- lapply(delay_grid, fit_one)
  logliks <- vapply(fits, function(f) if (is.null(f)) -Inf else as.numeric(logLik(f)),
                    numeric(1))
  ok <- is.finite(logliks)
  if (!any(ok)) stop("No model fits converged on the delay grid.")

  best_idx <- which.max(logliks)
  fit <- fits[[best_idx]]
  threshold <- max(logliks[ok]) - qchisq(conf_level, 1) / 2

  fit$delay_estimate <- delay_grid[best_idx]
  fit$shape <- shape
  fit$lambda <- lambda
  fit$conf_level <- conf_level
  fit$profile_loglik <- data.frame(delay = delay_grid, loglik = logliks)
  fit$delay_ci <- range(delay_grid[ok & logliks >= threshold])

  structure(fit, class = c("pmrm_delayed", class(fit)))
}

#' Summarize a delayed-effect fit
#'
#' @param object A \code{pmrm_delayed} object.
#' @param ... Unused.
#' @return List with fixed effects, delay estimate/CI, and information criteria.
#' @export
summary.pmrm_delayed <- function(object, ...) {
  fx <- tryCatch(nlme::fixed.effects(object),
                 error = function(e) tryCatch(stats::coef(object),
                                              error = function(e2) NULL))
  tt <- tryCatch(summary.lme(object)$tTable, error = function(e) NULL)
  out <- list(
    fixed_effects = fx,
    delay_estimate = object$delay_estimate,
    delay_ci = object$delay_ci,
    AIC = tryCatch(stats::AIC(object), error = function(e) NA_real_),
    BIC = tryCatch(stats::BIC(object), error = function(e) NA_real_)
  )
  if (!is.null(tt)) out$t_table <- tt
  out
}

#' Print method for delayed-effect fits
#' @param x A \code{pmrm_delayed} object.
#' @param ... Passed to print.
#' @export
print.pmrm_delayed <- function(x, ...) {
  cat("pmrmPlus delayed-effect progression model\n")
  cat("Progression shape :", attr(x, "progression"), "\n")
  cat("Estimated delay   :", round(x$delay_estimate, 3), "\n")
  if (!is.null(x$delay_ci)) {
    cat("Delay CI          : [", round(x$delay_ci[1], 3), ", ",
        round(x$delay_ci[2], 3), "]\n", sep = "")
  }
  invisible(x)
}

#' Plot a delayed-effect fit
#'
#' @param x A \code{pmrm_delayed} object.
#' @param type "trajectory" plots fitted population means by arm with the
#'   delay marked; "profile_likelihood" plots the log-likelihood vs delay
#'   with the chi-squared threshold line.
#' @param ... Unused.
#' @return A ggplot object (invisibly).
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline
#'   labs theme_minimal
#' @export
plot.pmrm_delayed <- function(x, type = c("trajectory", "profile_likelihood"), ...) {
  type <- match.arg(type)
  if (type == "profile_likelihood") {
    pl <- x$profile_loglik
    thr <- max(pl$loglik, na.rm = TRUE) - qchisq(x$conf_level %||% 0.95, 1) / 2
    p <- ggplot2::ggplot(pl, ggplot2::aes(x = delay, y = loglik)) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_point(color = "steelblue") +
      ggplot2::geom_vline(xintercept = x$delay_estimate, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = thr, linetype = "dotted") +
      ggplot2::labs(title = "Profile log-likelihood for delay",
                    subtitle = "Dashed: estimate; dotted: chi-square threshold",
                    x = "Delay", y = "Log-likelihood") +
      ggplot2::theme_minimal()
  } else {
    dat <- tryCatch(nlme::getData(x), error = function(e) NULL)
    if (is.null(dat)) dat <- x$data
    dat <- as.data.frame(dat)
    trt_col <- grep("^trt$", names(dat), ignore.case = TRUE)[1]
    dat$.pred <- stats::predict(x, newdata = dat, level = "population")
    p <- ggplot2::ggplot(dat,
                         ggplot2::aes(x = time, y = .pred,
                                      group = factor(.data[[trt_col]]),
                                      color = factor(.data[[trt_col]]))) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::geom_vline(xintercept = x$delay_estimate, linetype = "dashed") +
      ggplot2::labs(title = "Population trajectories by arm",
                    color = "Arm", y = "Outcome") +
      ggplot2::theme_minimal()
  }
  print(p)
  invisible(p)
}

