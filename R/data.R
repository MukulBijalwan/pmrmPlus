# Package data sets

#' Simulated ADNI-style Alzheimer's trial data
#'
#' Longitudinal data from a simulated ADNI-style Alzheimer's disease trial
#' with linear progression, a treatment effect with delayed onset, and
#' subject-level dropout information.
#'
#' @format A data frame with 2160 rows and 6 variables:
#' \describe{
#'   \item{rid}{Subject identifier.}
#'   \item{time}{Visit time in months.}
#'   \item{trt}{Treatment indicator: 0 = control, 1 = active.}
#'   \item{adas_cog}{ADAS-Cog score (higher values indicate worse cognition).}
#'   \item{time_to_dropout}{Subject-level time to dropout.}
#'   \item{dropout}{Dropout indicator: 1 = dropped out, 0 = censored.}
#' }
#' @source Simulated with \code{simulate_progression_trial()}; see
#'   \code{data-raw/create_data.R}.
#' @examples
#' head(adni_simulated)
#' table(adni_simulated$trt)
#' summary(adni_simulated$adas_cog)
"adni_simulated"

#' Simulated Parkinson's trial data with informative dropout
#'
#' Longitudinal UPDRS data from a simulated Parkinson's disease trial in
#' which progression drives an informative dropout process, so that faster
#' progressors leave the trial earlier. Intended for use with
#' \code{pmrm_joint()} and \code{pmrm_competing()}.
#'
#' @format A data frame with 1800 rows and 6 variables:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{time}{Visit time in months.}
#'   \item{trt}{Treatment indicator: 0 = control, 1 = active.}
#'   \item{updrs}{UPDRS score (higher values indicate worse symptoms).}
#'   \item{time_to_dropout}{Subject-level time to dropout.}
#'   \item{dropout}{Dropout indicator: 1 = dropped out, 0 = censored.}
#' }
#' @source Simulated with \code{simulate_joint_trial()}; see
#'   \code{data-raw/create_data.R}.
#' @examples
#' head(parkinsons_simulated)
#' with(parkinsons_simulated, table(trt, dropout))
#' summary(parkinsons_simulated$updrs)
"parkinsons_simulated"
