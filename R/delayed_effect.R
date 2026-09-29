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
#' @examples
#' t <- seq(0, 12, by = 2)
#' delayed_effect(t, delay = 4, shape = "step")
#' delayed_effect(t, delay = 0, shape = "exponential", lambda = 0.1)
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

# Drop the treatment terms from a user-supplied model formula so that the
# treatment effect is carried exclusively by the delayed term `trt_eff`.
# Leaving unrestricted treatment terms in place alongside the delayed term
# makes the delay unidentifiable: the profile likelihood then always peaks at
# delay = 0, whatever the data.
strip_treatment_terms <- function(formula, trt_col) {
  tl <- attr(stats::terms(formula), "term.labels")
  if (length(tl)) {
    keeps <- !vapply(tl, function(tm) {
      trt_col %in% all.vars(stats::as.formula(paste("~", tm)))
    }, logical(1))
    tl <- tl[keeps]
  }
  resp <- all.vars(stats::update(formula, . ~ 1))[1]
  stats::as.formula(paste(resp, "~", paste(c(tl, "trt_eff"), collapse = " + ")),
                    env = environment(formula))
}

# Report why fitting failed; the individual grid fits are wrapped in
# tryCatch(), so without this the user only ever sees a generic message.
stop_fit_failed <- function(what, first_error = NULL) {
  msg <- what
  if (!is.null(first_error)) msg <- paste0(msg, " First error: ", first_error)
  stop(msg, call. = FALSE)
}

#' Fit a progression model with delayed treatment onset
#'
#' The treatment effect enters the model exclusively through a delayed term,
#' \code{delta * Trt_i * t_ij * d(t_ij)}, where \code{d()} is a fractional
#' onset function: treatment terms already present in \code{formula} (such as
#' the \code{trt} and \code{time:trt} terms of \code{y ~ time * trt}) are
#' replaced by this term. Keeping them alongside the delayed term would let an
#' unrestricted treatment effect absorb the delay and make it unidentifiable,
#' in which case the profile likelihood always peaks at zero delay. Because
#' the delayed term is proportional to \code{time}, \code{delta} is the
#' treatment effect on the progression slope at full onset.
#'
#' The data must contain a treatment indicator column named \code{trt} (any
#' case), which is what identifies the treatment terms to replace. Pass other
#' \code{pmrm_nlme()} options -- most importantly \code{subject_var} and
#' \code{time_var} -- through \code{...}.
#'
#' When \code{delay_estimation = "profile"}, the delay is estimated by grid
#' search over \code{delay_grid}, selecting the value maximizing
#' log-likelihood; an approximate confidence set is derived from the
#' chi-squared threshold (log-lik drop of qchisq(0.95, 1)/2). All candidates
#' use the same fixed-effects structure, so the REML log-likelihoods are
#' comparable across the grid.
#'
#' @param formula LME-style formula, e.g. \code{y ~ time * trt}. Terms
#'   involving the treatment indicator are replaced by the delayed treatment
#'   term; the remaining terms specify the common progression model.
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
#' @examples
#' sim <- simulate_progression_trial(n_control = 40, n_treatment = 40,
#'                                   effect_size = -0.8, delay = 6,
#'                                   times = seq(0, 24, by = 6),
#'                                   residual_sd = 1.5, random_sd = 0.2, seed = 7)
#' fit <- pmrm_delayed(y ~ time * trt, data = sim,
#'                     delay_estimation = "profile",
#'                     delay_grid = c(0, 6, 12), lambda = 0.2)
#' fit$delay_estimate
#' summary(fit)$delay_ci
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

  data_names <- names(as.data.frame(data))
  trt_hit <- grep("^trt$", data_names, ignore.case = TRUE)
  if (!length(trt_hit)) {
    stop("Could not find a treatment column named 'trt' in 'data'.")
  }
  trt_col <- data_names[trt_hit[1]]
  # treatment terms are stripped so that the delay stays identifiable
  model_formula <- strip_treatment_terms(formula, trt_col)
  first_error <- NULL

  fit_one <- function(delay) {
    d <- as.data.frame(data)
    # effect on the progression slope, modulated by the onset fraction
    d$trt_eff <- as.numeric(d[[trt_col]]) * d$time *
      delayed_effect(d$time, delay = delay, shape = shape, lambda = lambda)
    f <- tryCatch(
      suppressWarnings(pmrm_nlme(model_formula, data = d,
                                 progression = progression, ...)),
      error = function(e) {
        if (is.null(first_error)) first_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!inherits(f, "lme")) NULL else f
  }

  if (delay_estimation == "fixed") {
    fit <- fit_one(0)
    if (is.null(fit)) {
      stop_fit_failed("The delayed-effect model failed to fit.", first_error)
    }
    fit$delay_estimate <- 0
    return(structure(fit, class = c("pmrm_delayed", class(fit))))
  }

  fits <- lapply(delay_grid, fit_one)
  logliks <- vapply(fits, function(f) if (is.null(f)) -Inf else as.numeric(logLik(f)),
                    numeric(1))
  ok <- is.finite(logliks)
  if (!any(ok)) {
    stop_fit_failed("No model fits converged on the delay grid.", first_error)
  }

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
#' @examples
#' sim <- simulate_progression_trial(n_control = 40, n_treatment = 40,
#'                                   effect_size = -0.8, delay = 6,
#'                                   times = seq(0, 24, by = 6), seed = 7)
#' fit <- pmrm_delayed(y ~ time * trt, data = sim,
#'                     delay_estimation = "profile",
#'                     delay_grid = c(0, 6), lambda = 0.2)
#' s <- summary(fit)
#' s$delay_estimate
#' s$delay_ci
#' @export
summary.pmrm_delayed <- function(object, ...) {
  fx <- tryCatch(nlme::fixed.effects(object),
                 error = function(e) tryCatch(stats::coef(object),
                                              error = function(e2) NULL))
  # `summary()` would dispatch straight back to this method, so drop the
  # "pmrm_delayed" class first to reach the underlying summary.lme() result.
  tt <- tryCatch({
    plain <- object
    class(plain) <- setdiff(class(plain), "pmrm_delayed")
    summary(plain)$tTable
  }, error = function(e) NULL)
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
#' @examples
#' sim <- simulate_progression_trial(n_control = 40, n_treatment = 40,
#'                                   effect_size = -0.8, delay = 6,
#'                                   times = seq(0, 24, by = 6), seed = 7)
#' fit <- pmrm_delayed(y ~ time * trt, data = sim,
#'                     delay_estimation = "profile",
#'                     delay_grid = c(0, 6), lambda = 0.2)
#' print(fit)
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
#' @examples
#' sim <- simulate_progression_trial(n_control = 40, n_treatment = 40,
#'                                   effect_size = -0.8, delay = 6,
#'                                   times = seq(0, 24, by = 6), seed = 7)
#' fit <- pmrm_delayed(y ~ time * trt, data = sim,
#'                     delay_estimation = "profile",
#'                     delay_grid = c(0, 6), lambda = 0.2)
#' plot(fit, type = "profile_likelihood")
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
    trt_hit <- grep("^trt$", names(dat), ignore.case = TRUE)
    if (!length(trt_hit)) {
      stop("Could not find a treatment column named 'trt' in the data.")
    }
    trt_col <- names(dat)[trt_hit[1]]
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

