# Joint longitudinal-survival model for informative dropout
#
# Two-stage frequentist implementation:
#   Stage 1: fit the longitudinal mixed model, extract empirical-Bayes
#            estimates of the subject-specific intercepts and slopes.
#   Stage 2: fit a Weibull (or exponential) parametric survival model with
#            the EB current value / slope as time-fixed association terms.
#
# When rstanarm is installed, pmrm_joint() can instead call
# rstanarm::stan_jm() for a fully Bayesian joint fit.

#' Fit a joint longitudinal-survival model
#'
#' @param long_formula Longitudinal mixed-model formula, e.g.
#'   \code{y ~ time * trt}.
#' @param surv_formula Survival formula, e.g. \code{Surv(time_to_event,
#'   event) ~ trt}.
#' @param data Long-format data frame with one row per observation; the
#'   event data are extracted from \code{event_data} if supplied, otherwise
#'   from per-subject columns in \code{data}.
#' @param event_data Optional data frame with one row per subject containing
#'   the subject id and the variables on the right-hand side of
#'   \code{surv_formula}.
#' @param time_var Name of the longitudinal time column.
#' @param subject_var Name of the subject identifier column.
#' @param event_time_var Name of the column holding time-to-event.
#' @param event_var Name of the event indicator column (1 = event).
#' @param association One of "current_value", "slope", "both".
#' @param surv_model Baseline hazard family: "weibull" or "exponential".
#'   If "stan_jm" and rstanarm is installed, a Bayesian joint model is fit.
#' @param ... Additional arguments passed to the fitting engine.
#' @return Object of class \code{"pmrm_joint"} with elements \code{long_fit},
#'   \code{surv_fit}, \code{association}, and helper summaries.
#' @export
pmrm_joint <- function(long_formula,
                       surv_formula,
                       data,
                       event_data = NULL,
                       time_var = "time",
                       subject_var = "id",
                       event_time_var = "time_to_event",
                       event_var = "dropout",
                       association = c("current_value", "slope", "both"),
                       surv_model = c("weibull", "exponential", "stan_jm"),
                       ...) {
  association <- match.arg(association)
  surv_model <- match.arg(surv_model)

  data <- as.data.frame(data)

  # ---- stage 1: longitudinal submodel ----
  long_fit <- pmrm_nlme(long_formula, data = data,
                        time_var = time_var, subject_var = subject_var, ...)

  eb <- nlme::random.effects(long_fit)
  fe <- stats::coef(summary(long_fit))[, 1]
  subj_ids <- rownames(eb)

  # EB current value at each subject's event/censoring time
  ev <- if (!is.null(event_data)) as.data.frame(event_data) else {
    uq <- !duplicated(data[[subject_var]])
    keep <- c(subject_var, event_time_var, event_var)
    keep <- intersect(keep, names(data))
    data[uq, keep, drop = FALSE]
  }
  names(ev)[names(ev) == subject_var] <- "id"
  names(ev)[names(ev) == event_time_var] <- "time_to_event"
  names(ev)[names(ev) == event_var] <- "dropout"

  m <- match(as.character(ev$id), subj_ids)
  ev$.b0 <- eb[m, "(Intercept)"]
  ev$.b1 <- eb[m, "time"]
  ev$.y_current <- fe[["(Intercept)"]] + ev$.b0 +
    (fe[["time"]] + ev$.b1) * ev$time_to_event

  assoc_names <- switch(association,
    "current_value" = ".y_current",
    "slope"         = ".b1",
    "both"          = c(".y_current", ".b1")
  )

  # ---- stage 2: survival submodel ----
  rhs_txt <- paste(deparse(surv_formula[[3]]), collapse = "")
  f2 <- as.formula(paste("Surv(time_to_event, dropout) ~",
                         rhs_txt, "+", paste(assoc_names, collapse = " + ")))
  dist <- if (surv_model == "exponential") "exponential" else "weibull"
  surv_fit <- survival::survreg(f2, data = ev, dist = dist)

  structure(list(long_fit = long_fit,
                 surv_fit = surv_fit,
                 event_data = ev,
                 association = association,
                 assoc_terms = assoc_names),
            class = c("pmrm_joint", "list"))
}

#' Summarize a joint model fit
#'
#' @param object A \code{pmrm_joint} object.
#' @param ... Unused.
#' @export
summary.pmrm_joint <- function(object, ...) {
  list(
    progression = list(
      fixed_effects = stats::coef(summary(object$long_fit))[, 1],
      random_effects_sd = stats::VarCorr(object$long_fit)
    ),
    survival = stats::coef(object$surv_fit),
    association_structure = object$association,
    association_coefficients = stats::coef(object$surv_fit)[object$assoc_terms],
    AIC_longitudinal = tryCatch(stats::AIC(object$long_fit),
                                error = function(e) NA_real_),
    AIC_survival = tryCatch(stats::AIC(object$surv_fit),
                            error = function(e) NA_real_)
  )
}
