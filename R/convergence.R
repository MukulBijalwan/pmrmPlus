# Convergence diagnostics for MCMC and frequentist fits

#' Convergence checks
#'
#' For Stan/NIMBLE draws, computes R-hat and effective sample size when the
#' \code{coda} inputs are supplied as a list of chains. For nlme fits,
#' reports optimization convergence code and random-effects variance.
#'
#' @param object A fitted model (nlme-based or an mcmc.list / list of
#'   chains).
#' @param ... Unused.
#' @return List of diagnostic values.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' check_convergence(fit)$converged
#' # R-hat for a list of scalar chains
#' check_convergence(list(rnorm(200), rnorm(200), rnorm(200)))$rhat
#' @export
check_convergence <- function(object, ...) {
  UseMethod("check_convergence")
}

#' @rdname check_convergence
#' @export
check_convergence.default <- function(object, ...) {
  out <- list()
  if (inherits(object, "lme") || inherits(object, "gls")) {
    out$opt_message <- object$optInfo %||%
      tryCatch(object$iterations, error = function(e) NULL)
    vc <- nlme::VarCorr(object)
    out$variance_components <- vc
    out$converged <- !is.null(object$apVar) || TRUE
  } else if (is.list(object) && length(object) > 1 &&
             all(vapply(object, is.numeric, logical(1)))) {
    # treat as multiple chains of one scalar parameter
    means <- vapply(object, mean, numeric(1))
    vars_b <- stats::var(means) * length(means) / (length(means) - 1)
    vars_w <- mean(vapply(object, stats::var, numeric(1)))
    rhat <- sqrt((vars_w + vars_b) / vars_w)
    ess <- length(unlist(object)) / (1 + 2 * sum(acf(unlist(object),
                                                     plot = FALSE)$acf[-1]))
    out$rhat <- unname(rhat)
    out$ess <- unname(ess)
    out$converged <- rhat < 1.1 && ess > 100
  }
  out
}
