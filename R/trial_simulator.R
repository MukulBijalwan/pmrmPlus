# End-to-end clinical trial simulator

#' Simulate a progression trial
#'
#' @param n_control,n_treatment Subjects per arm.
#' @param effect_size Treatment effect on the slope per unit of time
#'   (negative = slower progression).
#' @param delay Delayed onset (same units as \code{times}).
#' @param dropout_rate Annual dropout probability.
#' @param times Vector of measurement times.
#' @param intercept Population intercept at time 0.
#' @param slope Control-arm slope per unit time.
#' @param random_sd SD of subject random slope.
#' @param residual_sd Residual SD.
#' @param seed Random seed.
#' @return Long-format data frame with columns \code{id, time, trt, y}.
#' @examples
#' sim <- simulate_progression_trial(n_control = 20, n_treatment = 20,
#'                                   effect_size = -0.5, delay = 6,
#'                                   times = seq(0, 24, by = 6), seed = 1)
#' head(sim)
#' aggregate(y ~ trt, data = sim, FUN = mean, na.rm = TRUE)
#' @export
simulate_progression_trial <- function(n_control = 50, n_treatment = 50,
                                       effect_size = -0.5,
                                       delay = 0,
                                       dropout_rate = 0.15,
                                       times = seq(0, 24, by = 6),
                                       intercept = 20,
                                       slope = 1,
                                       random_sd = 0.3,
                                       residual_sd = 2,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- n_control + n_treatment
  trt <- c(rep(0, n_control), rep(1, n_treatment))
  b0i <- rnorm(n, 0, random_sd * sqrt(times[2] - times[1]))
  b1i <- rnorm(n, 0, random_sd)

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    # exponential onset fraction; full effect realized asymptotically
    frac <- if (delay > 0) 1 - exp(-0.2 * pmax(0, times - delay)) else rep(1, length(times))
    mu <- intercept + b0i[i] +
      (slope + b1i[i]) * times +
      trt[i] * effect_size * times * frac
    keep <- runif(length(times)) >= dropout_rate / length(times)
    y <- ifelse(keep, rnorm(length(times), mu, residual_sd), NA)
    rows[[i]] <- data.frame(id = i, time = times, trt = trt[i], y = y)
  }
  do.call(rbind, rows)
}

#' Simulate a trial with informative dropout via a joint model mechanism
#'
#' Progression drives a Weibull hazard for dropout through the current
#' latent value, so faster progressors leave earlier.
#'
#' @inheritParams simulate_progression_trial
#' @param alpha Association of current value with the log-hazard.
#' @param weibull_scale,weibull_shape Weibull parameters for baseline hazard.
#' @examples
#' sim <- simulate_joint_trial(n_control = 20, n_treatment = 20,
#'                             times = seq(0, 24, by = 6), seed = 1)
#' head(sim)
#' # share of subjects with an observed dropout event
#' mean(tapply(sim$dropout, sim$id, function(z) z[1]))
#' @export
simulate_joint_trial <- function(n_control = 50, n_treatment = 50,
                                 effect_size = -0.5, delay = 0,
                                 times = seq(0, 24, by = 3),
                                 intercept = 20, slope = 1,
                                 random_sd = 0.5, residual_sd = 2,
                                 alpha = 0.05,
                                 weibull_scale = 100, weibull_shape = 1.4,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- n_control + n_treatment
  trt <- c(rep(0, n_control), rep(1, n_treatment))
  b0i <- rnorm(n, 0, random_sd)
  b1i <- rnorm(n, 0, random_sd / 10)

  surv_time <- numeric(n); event <- integer(n)
  for (i in seq_len(n)) {
    # approximate hazard integral over a fine grid using latent value
    tt <- seq(0.25, max(times) + 12, by = 0.25)
    lat <- intercept + b0i[i] + (slope + b1i[i]) * tt +
      trt[i] * effect_size * tt
    hz <- weibull_shape / weibull_scale *
      (tt / weibull_scale)^(weibull_shape - 1) * exp(alpha * lat)
    H <- cumsum(hz) * 0.25
    u <- runif(1)
    idx <- which(H >= -log(u))[1]
    if (is.na(idx)) { surv_time[i] <- max(times); event[i] <- 0 }
    else { surv_time[i] <- tt[idx]; event[i] <- 1 }
  }

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    keep <- times <= surv_time[i]
    frac <- if (delay > 0) 1 - exp(-0.2 * pmax(0, times - delay)) else rep(1, length(times))
    mu <- intercept + b0i[i] + (slope + b1i[i]) * times +
      trt[i] * effect_size * times * frac
    y <- ifelse(keep, rnorm(length(times), mu, residual_sd), NA)
    rows[[i]] <- data.frame(id = i, time = times, trt = trt[i], y = y,
                            time_to_dropout = surv_time[i], dropout = event[i])
  }
  do.call(rbind, rows)
}
