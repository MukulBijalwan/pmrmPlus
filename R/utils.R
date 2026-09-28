# Shared utilities

#' Extract subject-level summary (EB intercept and slope)
#'
#' @param fit A \code{pmrm_nlme} object.
#' @return Data frame with \code{id}, \code{b0} and \code{b1}.
#' @export
subject_summaries <- function(fit) {
  eb <- nlme::random.effects(fit)
  data.frame(id = rownames(eb), b0 = eb[, "(Intercept)"], b1 = eb[, "time"])
}

#' Reshape long trial data to per-subject event data
#'
#' @param long_data Long-format data frame.
#' @param id_var,time_to_event_var,event_var Column names.
#' @return One-row-per-subject data frame.
#' @export
to_event_data <- function(long_data,
                          id_var = "id",
                          time_to_event_var = "time_to_dropout",
                          event_var = "dropout") {
  uq <- !duplicated(long_data[[id_var]])
  keep <- intersect(c(id_var, time_to_event_var, event_var), names(long_data))
  long_data[uq, keep, drop = FALSE]
}

#' Extract the right-hand side of a survival formula as text
#'
#' Accepts one-sided (\code{~ x}) or two-sided (\code{Surv(...) ~ x})
#' formulas and returns the right-hand side, so callers can rebuild a
#' survival formula with additional association terms.
#'
#' @param f A formula.
#' @return A character string with the right-hand side (defaults to \code{"1"}).
#' @noRd
surv_rhs_text <- function(f) {
  if (length(f) == 3L) return(paste(deparse(f[[3L]]), collapse = ""))
  if (length(f) == 2L) return(paste(deparse(f[[2L]]), collapse = ""))
  "1"
}

`%||%_u` <- NULL  # reserved namespace placeholder

# Names referenced through ggplot2 non-standard evaluation.
utils::globalVariables(c(".data", ".pred", "delay", "loglik", "time"))
