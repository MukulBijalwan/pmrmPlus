# Biomarker-informed adaptive dosing

#' PK/PD helpers for biomarker dosing
#'
#' One-compartment oral PK concentration and an Emax biomarker response.
#'
#' @param time Time(s) since dose.
#' @param dose Administered dose.
#' @param ka Absorption rate constant.
#' @param ke Elimination rate constant.
#' @param V Volume of distribution.
#' @return Concentration at \code{time}.
#' @export
pkpd_concentration <- function(time, dose, ka = 0.5, ke = 0.1, V = 10) {
  (dose * ka) / (V * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
}

#' Emax biomarker model
#'
#' @param conc Drug concentration.
#' @param b_baseline Baseline biomarker level.
#' @param emax Maximum effect.
#' @param ec50 Concentration giving half-maximal effect.
#' @export
emax_biomarker <- function(conc, b_baseline = 100, emax = 30, ec50 = 5) {
  b_baseline + emax * conc / (ec50 + conc)
}

#' Biomarker-informed adaptive dosing rule
#'
#' Adjusts the dose one step up/down when the measured biomarker falls
#' outside the target range, otherwise maintains the current dose.
#'
#' @param biomarker_value Observed biomarker value.
#' @param target_range Two numeric values giving the acceptable range.
#' @param dose_levels Ordered vector of allowed doses.
#' @param current_dose_idx Index of the subject's current dose.
#' @return List with \code{new_dose}, \code{new_dose_idx}, and
#'   \code{adjustment} ("up", "down", or "maintain").
#' @export
biomarker_dosing_rule <- function(biomarker_value,
                                  target_range = c(80, 120),
                                  dose_levels = c(0, 50, 100, 150, 200),
                                  current_dose_idx = 3) {
  if (!is.numeric(current_dose_idx)) {
    current_dose_idx <- which(dose_levels == current_dose_idx)
    if (length(current_dose_idx) == 0) stop("current_dose_idx not found in dose_levels")
  }
  lo <- min(target_range); hi <- max(target_range)
  new_idx <- if (biomarker_value > hi) {
    min(current_dose_idx + 1, length(dose_levels))
  } else if (biomarker_value < lo) {
    max(current_dose_idx - 1, 1)
  } else {
    current_dose_idx
  }
  list(new_dose = dose_levels[new_idx],
       new_dose_idx = new_idx,
       adjustment = c("down", "maintain", "up")[sign(new_idx - current_dose_idx) + 2])
}

#' Initialize subjects for an adaptive trial simulation
#'
#' @param n_subjects Number of subjects.
#' @param n_visits Number of planned visits.
#' @param dose_levels Allowed dose levels.
#' @param start_dose_idx Starting dose index.
#' @param trt_prob Probability of assignment to active treatment.
#' @param seed Random seed.
#' @return Data frame of subject-level attributes.
#' @export
initialize_subjects <- function(n_subjects, n_visits = 12,
                                dose_levels = c(0, 50, 100, 150, 200),
                                start_dose_idx = 3, trt_prob = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    id = seq_len(n_subjects),
    trt = sample(0:1, n_subjects, replace = TRUE, prob = c(1 - trt_prob, trt_prob)),
    dose_idx = rep(start_dose_idx, n_subjects),
    discontinued = rep(FALSE, n_subjects),
    stringsAsFactors = FALSE
  )
}

#' Simulate biomarker measurements under the PK/PD model
#'
#' @param subjects Subject data frame from \code{initialize_subjects()}.
#' @param visit Current visit index.
#' @param visit_interval Interval between visits.
#' @param dose_levels Dose levels vector.
#' @param noise_sd Residual measurement SD.
#' @param emax,ec50,b_baseline Emax-model parameters.
#' @return Named numeric vector of observed biomarkers.
#' @export
simulate_biomarker <- function(subjects, visit, visit_interval = 4,
                               dose_levels = c(0, 50, 100, 150, 200),
                               noise_sd = 3,
                               emax = 30, ec50 = 5, b_baseline = 100) {
  doses <- dose_levels[pmin(subjects$dose_idx, length(dose_levels))]
  doses <- ifelse(subjects$trt == 1, doses, 0)
  conc <- pkpd_concentration(visit * visit_interval, doses)
  obs <- emax_biomarker(conc, b_baseline = b_baseline, emax = emax, ec50 = ec50)
  setNames(rnorm(length(obs), obs, noise_sd), subjects$id)
}

#' Run a biomarker-informed adaptive trial simulation
#'
#' Simulates visits sequentially: measure biomarkers with noise, apply the
#' dosing rule, simulate progression and terminal events.
#'
#' @param n_subjects Number of subjects.
#' @param progression_model List describing progression:
#'   intercept, slope, treatment effect on slope, residual sd.
#' @param pkpd_model List with elements \code{emax}, \code{ec50},
#'   \code{b_baseline}, \code{noise_sd}.
#' @param dosing_rule Function like \code{biomarker_dosing_rule()}.
#' @param n_visits Number of visits.
#' @param visit_interval Weeks between visits.
#' @param dropout_rate Annual dropout hazard.
#' @param target_range Target biomarker range for dosing.
#' @param dose_levels Allowed dose levels.
#' @param seed Random seed.
#' @return A list with long-format outcomes \code{$long_data} and final
#'   subject states \code{$subjects}.
#' @export
pmrm_adaptive_trial <- function(n_subjects,
                                progression_model = list(intercept = 20, slope = 1.5,
                                                         trt_effect = -0.6, sigma = 2),
                                pkpd_model = list(emax = 30, ec50 = 5,
                                                  b_baseline = 100, noise_sd = 3),
                                dosing_rule = biomarker_dosing_rule,
                                n_visits = 12,
                                visit_interval = 4,
                                dropout_rate = 0.15,
                                target_range = c(pkpd_model$b_baseline - 10,
                                                 pkpd_model$b_baseline + 10),
                                dose_levels = c(0, 50, 100, 150, 200),
                                seed = NULL) {
  subjects <- initialize_subjects(n_subjects, n_visits = n_visits, seed = seed)
  rows <- vector("list", n_visits * n_subjects)
  k <- 0
  for (v in seq_len(n_visits)) {
    t_v <- v * visit_interval
    bio <- simulate_biomarker(subjects, v, visit_interval = visit_interval,
                              dose_levels = dose_levels,
                              emax = pkpd_model$emax, ec50 = pkpd_model$ec50,
                              b_baseline = pkpd_model$b_baseline,
                              noise_sd = pkpd_model$noise_sd)
    for (i in seq_len(n_subjects)) {
      if (!subjects$discontinued[i] && subjects$trt[i] == 1) {
        dec <- dosing_rule(bio[i], target_range = target_range,
                           current_dose_idx = subjects$dose_idx[i])
        subjects$dose_idx[i] <- dec$new_dose_idx
      }
      if (!subjects$discontinued[i]) {
        haz <- dropout_rate * visit_interval / 52
        if (runif(1) < haz) subjects$discontinued[i] <- TRUE
      }
    }
    dose_now <- dose_levels[pmin(subjects$dose_idx, length(dose_levels))]
    eff_frac <- delayed_effect(t_v, shape = "exponential", lambda = 0.05)
    y_mean <- progression_model$intercept + progression_model$slope * t_v +
      subjects$trt * progression_model$trt_effect * t_v / 24 * eff_frac
    for (i in seq_len(n_subjects)) {
      y_obs <- if (!subjects$discontinued[i])
        rnorm(1, y_mean[i], progression_model$sigma) else NA_real_
      k <- k + 1
      rows[[k]] <- data.frame(id = i, visit = v, time = t_v,
                              trt = subjects$trt[i],
                              dose = dose_now[i],
                              biomarker = unname(bio[i]),
                              y = y_obs,
                              dropped = subjects$discontinued[i])
    }
  }
  list(long_data = do.call(rbind, rows), subjects = subjects)
}

