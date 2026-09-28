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
#' @param event_time_var Column holding time-to-event (defaults to the column
#'   name used by the shipped data sets).
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
                           event_time_var = "time_to_dropout",
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

  assoc_names <- switch(association,
    "current_value" = ".y_current",
    "slope"         = ".b1",
    "both"          = c(".y_current", ".b1"))

  # ---- normalise the survival formulas ----
  n_causes <- length(causes)
  if (!is.list(surv_formulas)) {
    base_rhs <- surv_rhs_text(surv_formulas)
    surv_formulas <- setNames(replicate(n_causes, as.formula(paste("~", base_rhs)),
                                        simplify = FALSE), causes)
  }
  stopifnot(length(surv_formulas) == n_causes)
  rhs_txts <- vapply(surv_formulas, surv_rhs_text, character(1))

  # ---- per-subject event data (including survival covariates) ----
  covars <- unique(unlist(lapply(surv_formulas, all.vars)))
  ev <- data[!duplicated(data[[subject_var]]),
             unique(intersect(c(subject_var, event_time_var, event_var, covars),
                              names(data))),
             drop = FALSE]
  names(ev)[names(ev) == subject_var] <- "id"
  names(ev)[names(ev) == event_time_var] <- "time_to_event"
  names(ev)[names(ev) == event_var] <- "event"
  if (!all(c("id", "time_to_event", "event") %in% names(ev))) {
    stop("Could not locate the subject / event-time / event columns in the ",
         "data. Set 'subject_var', 'event_time_var' and 'event_var' to match ",
         "your column names.")
  }
  ev$event <- as.integer(ev$event)

  m <- match(as.character(ev$id), subj_ids)
  # The longitudinal submodel may fall back to a random intercept only, in
  # which case the empirical-Bayes random slope is zero.
  ev$.b0 <- if ("(Intercept)" %in% colnames(eb)) eb[m, "(Intercept)"] else 0
  ev$.b1 <- if ("time" %in% colnames(eb)) eb[m, "time"] else 0
  fe_b0 <- if ("(Intercept)" %in% names(fe)) fe[["(Intercept)"]] else 0
  fe_b1 <- if ("time" %in% names(fe)) fe[["time"]] else 0
  ev$.y_current <- fe_b0 + ev$.b0 + (fe_b1 + ev$.b1) * ev$time_to_event

  # ---- cause-specific hazards ----
  cause_fits <- lapply(seq_len(n_causes), function(k) {
    ek <- as.integer(ev$event == k)
    f_k <- as.formula(paste("Surv(time_to_event,", paste0("ev_", k), ") ~",
                            rhs_txts[[k]], "+",
                            paste(assoc_names, collapse = " + ")))
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
