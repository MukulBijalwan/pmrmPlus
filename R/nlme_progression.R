# Nonlinear mixed-effects progression models
#
# y_ij = f(t_ij, beta, b_i) + eps_ij,  b_i ~ N(0, Omega), eps_ij ~ N(0, sigma^2)

#' Fit a nonlinear progression mixed-effects model
#'
#' Fits repeated-measures progression data using one of four common shapes
#' (linear, plateau, emax, exponential). Linear models use \code{nlme::lme};
#' nonlinear shapes use \code{nlme::nlme}.
#'
#' @param formula Formula such as \code{y ~ time + trt} (linear shape).
#' @param data Data frame in long format.
#' @param progression Progression shape: one of "linear", "plateau", "emax",
#'   "exponential".
#' @param random Random effects specification; the grouping variable is taken
#'   from \code{subject_var}.
#' @param treatment_var Name of the treatment indicator column.
#' @param time_var Name of the time column.
#' @param subject_var Name of the subject identifier column.
#' @param heteroscedastic If TRUE allow residual variance to increase with
#'   time via a power variance function.
#' @param start Optional named starting values for nonlinear shapes.
#' @param ... Additional arguments passed on to the fitting engine.
#' @return An object of class \code{"pmrm_nlme"}.
#' @references Pinheiro JC, Bates DM (2000). Mixed-Effects Models in S and
#'   S-PLUS. Springer Statistics and Computing. ISBN 978-0-387-98957-0.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' summary(fit)$tTable
#' @export
pmrm_nlme <- function(formula,
                      data,
                      progression = c("linear", "plateau", "emax", "exponential"),
                      random = NULL,
                      treatment_var = "trt",
                      time_var = "time",
                      subject_var = "id",
                      heteroscedastic = FALSE,
                      start = NULL,
                      ...) {
  progression <- match.arg(progression)
  stopifnot(is.data.frame(data))

  data <- as.data.frame(data)
  names(data)[names(data) == time_var] <- "time"
  names(data)[names(data) == subject_var] <- "subject"

  weights <- if (isTRUE(heteroscedastic)) nlme::varPower(form = ~ time) else NULL

  fit <- switch(
    progression,
    "linear" = {
      resp <- all.vars(update(formula, . ~ 1))[1]
      rhs <- paste(deparse(formula[[3]]), collapse = " ")
      rhs_clean <- gsub(time_var, "time", rhs)
      rhs_clean <- gsub(paste0("\\b", subject_var, "\\b"), "subject", rhs_clean)
      f <- as.formula(paste(resp, "~", rhs_clean))
      rf <- as.formula("~ 1 + time | subject")
      ctl <- nlme::lmeControl(maxIter = 200, msMaxIter = 200,
                              msVerbose = FALSE, optMaxiter = 200)
      f_fit <- tryCatch(
        nlme::lme(f, data = data, random = rf, weights = weights,
                  na.action = na.omit, method = "REML", control = ctl, ...),
        error = function(e1) {
          # retry with simplified random effects (random intercept only)
          tryCatch(
            nlme::lme(f, data = data, random = ~ 1 | subject,
                      weights = weights, na.action = na.omit,
                      method = "REML", control = ctl, ...),
            error = function(e2) NULL
          )
        }
      )
      if (is.null(f_fit)) {
        stop("lme failed to converge for this dataset")
      }
      f_fit
    },
    {
      if (is.null(start)) start <- default_start(progression, data)
      nl_form <- switch(
        progression,
        "plateau"     = as.formula(y ~ b0 + b1 * pmin(time, tau)),
        "emax"        = as.formula(y ~ b0 + bmax * time / (ed50 + time)),
        "exponential" = as.formula(y ~ b0 + (binf - b0) * (1 - exp(-lambda * time)))
      )
      params <- setdiff(all.vars(nl_form[[3]]), "time")
      fixed_f <- lapply(params, function(p) as.formula(paste(p, "~ 1")))
      names(fixed_f) <- params
      dat <- data
      dat$y <- as.numeric(dat[[all.vars(formula[[2]])[1]]])
      start_vec <- unlist(lapply(params, function(p) as.numeric(start[[p]])))
      # nlme() expects a two-sided "params ~ 1" random formula (unlike lme()'s
      # one-sided "~ x | group" form); pdDiag() gives an uncorrelated
      # diagonal random-effects covariance matrix.
      ran_full <- nlme::pdDiag(
        as.formula(paste(paste(params, collapse = " + "), "~ 1")))
      ran_intercept <- nlme::pdDiag(as.formula(paste(params[1L], "~ 1")))
      tryCatch({
        nlme::nlme(nl_form, data = dat, fixed = fixed_f,
                   random = ran_full, groups = ~ subject,
                   start = start_vec, weights = weights,
                   na.action = na.omit, ...)
      }, error = function(e) {
        # retry with a random intercept on the baseline parameter only
        nlme::nlme(nl_form, data = dat, fixed = fixed_f,
                   random = ran_intercept, groups = ~ subject,
                   start = start_vec, weights = weights,
                   na.action = na.omit, ...)
      })
    }
  )

  structure(fit,
            class = c("pmrm_nlme", class(fit)),
            progression = progression,
            time_var = time_var,
            subject_var = subject_var,
            treatment_var = treatment_var)
}

default_start <- function(progression, data) {
  ycol <- intersect(c("y", grep("response|adas|updrs", names(data),
                                ignore.case = TRUE, value = TRUE))[1], names(data))
  y <- stats::na.omit(as.numeric(if (length(ycol)) data[[ycol[1]]] else rep(0, nrow(data))))
  b0 <- mean(y, na.rm = TRUE); s <- max(sd(y, na.rm = TRUE), 1e-3)
  switch(progression,
         "plateau"     = list(b0 = b0, b1 = -s, tau = max(data$time, na.rm = TRUE) / 2),
         "emax"        = list(b0 = b0, bmax = -s, ed50 = 1),
         "exponential" = list(b0 = b0, binf = b0 - s, lambda = 0.1),
         list())
}

#' Predict from a fitted progression model
#'
#' @param object A \code{pmrm_nlme} object.
#' @param newdata New data frame containing the predictor columns.
#' @param level Either "population" (fixed effects only) or "individual"
#'   (includes empirical-Bayes random effects).
#' @param ... Unused.
#' @return Numeric vector of predictions.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' nd <- data.frame(time = c(0, 6, 12, 24), trt = c(0, 0, 1, 1))
#' predict(fit, newdata = nd)                 # population level
#' head(predict(fit, level = "individual"))   # subject-specific
#' @export
predict.pmrm_nlme <- function(object, newdata, level = c("population", "individual"), ...) {
  level <- match.arg(level)
  prog <- attr(object, "progression")

  if (prog == "linear") {
    if (level == "individual" && missing(newdata)) {
      return(as.vector(stats::fitted(object)))
    }
    if (missing(newdata)) newdata <- tryCatch(nlme::getData(object), error = function(e) NULL)
    if (is.null(newdata)) stop("newdata required")
    newdata <- as.data.frame(newdata)
    mm <- stats::model.matrix(delete.response(terms(object)), data = newdata)
    # summary() may dispatch to a pmrm_* method that returns a plain list
    # (e.g. summary.pmrm_delayed()), so take the fixed effects straight from
    # the fitted lme object instead.
    fe <- nlme::fixed.effects(object)
    if (!identical(colnames(mm), names(fe))) fe <- fe[colnames(mm)]
    return(as.vector(mm %*% fe))
  }

  # nonlinear shapes
  if (missing(newdata)) newdata <- tryCatch(nlme::getData(object), error = function(e) NULL)
  if (is.null(newdata)) stop("newdata required")
  newdata <- as.data.frame(newdata)
  tt <- if ("time" %in% names(newdata)) newdata$time else rep(0, nrow(newdata))

  if (level == "population") {
    fe <- nlme::fixed.effects(object)
    g <- function(nm) if (nm %in% names(fe)) as.numeric(fe[[nm]]) else NA_real_
    pred <- switch(prog,
      "plateau"     = g("b0") + g("b1") * pmin(tt, g("tau")),
      "emax"        = g("b0") + g("bmax") * tt / (g("ed50") + tt),
      "exponential" = g("b0") + (g("binf") - g("b0")) * (1 - exp(-g("lambda") * tt))
    )
  } else {
    pred <- as.vector(stats::predict(object, newdata = newdata))
  }
  as.vector(pred)
}

#' Population-level trajectory over a time grid
#'
#' @param object A \code{pmrm_nlme} or \code{pmrm_delayed} object.
#' @param times Numeric vector of times.
#' @param trt Treatment arm value (default 0).
#' @return Data frame with columns \code{time}, \code{trt}, \code{predicted}.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20, seed = 1)
#' fit <- pmrm_nlme(y ~ time * trt, data = sim)
#' population_trajectory(fit, times = seq(0, 24, by = 6), trt = c(0, 1))
#' @export
population_trajectory <- function(object, times, trt = 0) {
  nd <- expand.grid(time = times, trt = trt)
  nd$predicted <- predict(object, newdata = nd, level = "population")
  nd
}

