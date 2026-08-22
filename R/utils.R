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

`%||%_u` <- NULL  # reserved namespace placeholder
