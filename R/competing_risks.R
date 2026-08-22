# Competing risks extension: cause-specific hazards for multiple terminal
# events (dropout, death, discontinuation) with a shared longitudinal
# submodel.

#' Fit a competing-risks progression model
#'
#' Fits cause-specific Cox proportional-hazards models for each terminal
#' event type on top of the shared longitudinal mixed-effects submodel.
#' Association with the longitudinal process enters through empirical-Bayes
#' current value and slope at each subject's event/censoring time.
#'
#' @param long_formula Longitudinal formula.
#' @param surv_formulas Either a single survival formula or a named list of
#'   formulas, one per cause (right-hand side only; Surv response is added).
#' @param data Long-format data frame.
#' @param time_var Name of longitudinal time column.
#' @param subject_var Name of subject id column.
#' @param event_time_var Column holding time-to-event.
#' @param event_var Column holding the event code: 0 = censored,
#'   1..K = cause.
#' @param causes Character vector of cause labels in code order.
#' @param association "current_value", "slope", or "both".
#' @param ... Additional arguments to \code{pmrm_nlme()}.
#' @return Object of class \code{"pmrm_competing"}.
#' @export
pmrm_competing <- function(long_formula,
                           surv_formulas = y ~ 1,
                           data,
                           time_var = "time",
                           subject_var = "id",
                           event_time_var = "time_to_event",
                           event_var = "event",
                           causes = c("dropout", "death"),
                           association = c("current_value", "slope", "both"),
                           ...) {
  association <- match.arg(association)
  data <- as.data.frame(data)

  # ---- shared longitudinal submodel ----
  long_fit <- pmrm_nlme(long_formula, data = data,
                        time_var = time_var, subject_var = subject_var, ...)
  eb <- nlme::random.effects(long_fit)
  fe <- stats::coef(summary(long_fit))[, 1]
  subj_ids <- rownames(eb)

  ev <- data[!duplicated(data[[subject_var]]),
             intersect(c(subject_var, event_time_var, event_var), names(data)),
             drop = FALSE]
  names(ev) <- c("id", "time_to_event", "event")
  ev$event <- as.integer(ev$event)

  m <- match(as.character(ev$id), subj_ids)
  ev$.b0 <- eb[m, "(Intercept)"]
  ev$.b1 <- eb[m, "time"]
  ev$.y_current <- fe[["(Intercept)"]] + ev$.b0 +
    (fe[["time"]] + ev$.b1) * ev$time_to_event

  assoc_names <- switch(association,
    "current_value" = ".y_current",
    "slope"         = ".b1",
    "both"          = c(".y_current", ".b1"))

  n_causes <- length(causes)

  if (!is.list(surv_formulas)) {
    base_rhs <- paste(deparse(surv_formulas[[3]]), collapse = "")
    surv_formulas <- setNames(replicate(n_causes, as.formula(paste("~", base_rhs)),
                                        simplify = FALSE), causes)
  }
  stopifnot(length(surv_formulas) == n_causes)

  # ---- cause-specific hazards ----
  cause_fits <- lapply(seq_len(n_causes), function(k) {
    ek <- as.integer(ev$event == k)
    rhs <- paste(deparse(surv_formulas[[k]][[2]]), collapse = "") # may be NULL
    rhs_txt <- tryCatch(paste(deparse(surv_formulas[[k]])[[2]], collapse = ""),
                        error = function(e) NULL)
    if (is.null(rhs_txt)) {
      ff <- surv_formulas[[k]]
      rhs_txt <- if (length(ff) == 3) paste(deparse(ff[[3]]), collapse = "") else "1"
    }
    f_k <- as.formula(paste("Surv(time_to_event,", paste0("ev_", k), ") ~",
                            rhs_txt, "+", paste(assoc_names, collapse = " + ")))
    env <- new.env(parent = environment())
    assign(paste0("ev_", k), ek, envir = env)
    environment(f_k) <- env
    survival::coxph(f_k, data = ev)
  })
  names(cause_fits) <- causes

  structure(list(long_fit = long_fit,
                 cause_fits = cause_fits,
                 event_data = ev,
                 causes = causes,
                 association = association),
            class = c("pmrm_competing", "list"))
}

#' Summarize a competing-risks fit
#'
#' @param object A \code{pmrm_competing} object.
#' @param ... Unused.
#' @export
summary.pmrm_competing <- function(object, ...) {
  list(
    fixed_effects = stats::coef(summary(object$long_fit))[, 1],
    causes = object$causes,
    hazard_ratios = lapply(object$cause_fits, function(f) {
      exp(stats::coef(f))
    }),
    concordance = vapply(object$cause_fits, function(f)
      tryCatch(unname(survival::concordance(f)$concordance),
               error = function(e) NA_real_), numeric(1))
  )
}
