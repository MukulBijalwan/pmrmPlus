# Simulation-based power and assurance

#' Simulation-based power for progression trials
#'
#' @param n_per_arm Subjects per arm.
#' @param effect_size Treatment effect on slope (per year).
#' @param delay Delayed onset in months.
#' @param dropout_rate Annual dropout rate.
#' @param measurement_times Visit schedule in months.
#' @param n_sims Number of simulations.
#' @param alpha One-sided... actually two-sided type-I error level.
#' @param seed Random seed.
#' @param ... Passed to \code{simulate_progression_trial()}.
#' @return List with power, mean estimate/SE, bias, MSE, and the
#'   non-convergence rate.
#' @export
power_pmrm <- function(n_per_arm,
                       effect_size = -0.5,
                       delay = 0,
                       dropout_rate = 0.15,
                       measurement_times = seq(0, 24, by = 6),
                       n_sims = 200,
                       alpha = 0.05,
                       seed = NULL,
                       ...) {
  rejections <- 0
  estimates <- rep(NA_real_, n_sims)
  se_estimates <- rep(NA_real_, n_sims)

  for (s in seq_len(n_sims)) {
    sim_data <- simulate_progression_trial(
      n_control = n_per_arm, n_treatment = n_per_arm,
      effect_size = effect_size, delay = delay,
      dropout_rate = dropout_rate, times = measurement_times,
      seed = if (is.null(seed)) NULL else seed + s, ...
    )

    fit <- tryCatch(
      pmrm_nlme(y ~ time * trt, data = sim_data, progression = "linear"),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      tt <- tryCatch(summary(fit)$tTable, error = function(e) NULL)
      row_name <- grep("time:trt|trt:time", rownames(tt), value = TRUE)[1]
      if (!is.na(row_name)) {
        estimates[s] <- tt[row_name, "Value"]
        se_estimates[s] <- tt[row_name, "Std.Error"]
        p_val <- tt[row_name, "p-value"]
        if (p_val < alpha) rejections <- rejections + 1
      }
    }
  }

  list(
    power = rejections / n_sims,
    mean_estimate = mean(estimates, na.rm = TRUE),
    mean_se = mean(se_estimates, na.rm = TRUE),
    bias = mean(estimates - effect_size, na.rm = TRUE),
    mse = mean((estimates - effect_size)^2, na.rm = TRUE),
    non_convergence_rate = mean(is.na(estimates)),
    n_per_arm = n_per_arm
  )
}

#' Find sample size by simulation
#'
#' @param target_power Desired power.
#' @param n_grid Grid of per-arm sample sizes to evaluate.
#' @param ... Arguments passed to \code{power_pmrm()}.
#' @return List with \code{optimal_n} and a \code{power_curve} data frame.
#' @export
find_sample_size <- function(target_power = 0.8,
                             n_grid = seq(50, 300, by = 25),
                             ...) {
  results <- lapply(n_grid, function(n) {
    message("Simulating n = ", n)
    suppressWarnings(power_pmrm(n_per_arm = n, ...))
  })
  powers <- vapply(results, `[[`, numeric(1), "power")
  best_n <- n_grid[which.min(abs(powers - target_power))]
  list(optimal_n = best_n,
       power_curve = data.frame(n = n_grid, power = powers),
       results = results)
}

#' Scenario builders for trial simulation studies
#'
#' @description Named scenario lists consumed by \code{power_pmrm()} via
#'   \code{do.call}.
#' @return A list of simulation arguments.
#' @rdname scenarios
#' @export
scenario_base <- function() {
  list(effect_size = -0.5, delay = 6, dropout_rate = 0.15,
       measurement_times = seq(0, 24, by = 3))
}

#' @rdname scenarios
#' @export
scenario_optimistic <- function() {
  list(effect_size = -1.0, delay = 3, dropout_rate = 0.10,
       measurement_times = seq(0, 24, by = 3))
}

#' @rdname scenarios
#' @export
scenario_pessimistic <- function() {
  list(effect_size = -0.25, delay = 12, dropout_rate = 0.25,
       measurement_times = seq(0, 24, by = 3))
}
