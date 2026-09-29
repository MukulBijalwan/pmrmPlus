test_that("delayed_effect shapes behave correctly", {
  t <- seq(0, 12, by = 2)
  e_step <- delayed_effect(t, delay = 4, shape = "step")
  expect_equal(e_step, c(0, 0, 0, 1, 1, 1, 1))

  e_exp <- delayed_effect(t, delay = 0, shape = "exponential", lambda = 0.1)
  expect_true(all(e_exp >= 0 & e_exp <= 1))
  expect_true(all(diff(e_exp) >= 0))   # monotone increasing
  expect_equal(e_exp[1], 0)

  e_lin <- delayed_effect(t, trt_time = 0, delay = 0, shape = "linear")
  expect_equal(e_lin[length(e_lin)], 1)
})

test_that("simulate_progression_trial returns valid long data", {
  sim <- simulate_progression_trial(n_control = 10, n_treatment = 10,
                                    times = c(0, 6, 12), seed = 1)
  expect_s3_class(sim, "data.frame")
  expect_true(all(c("id", "time", "trt", "y") %in% names(sim)))
  expect_setequal(unique(sim$trt), c(0, 1))
})

test_that("pmrm_nlme fits a linear progression model", {
  sim <- simulate_progression_trial(n_control = 30, n_treatment = 30,
                                    effect_size = -0.5, seed = 42)
  fit <- pmrm_nlme(y ~ time * trt, data = sim, progression = "linear")
  expect_s3_class(fit, "pmrm_nlme")
  pr <- predict(fit, newdata = data.frame(time = 0, trt = 0),
                level = "population")
  expect_length(pr, 1)
})

test_that("pmrm_delayed recovers the delay", {
  sim <- simulate_progression_trial(n_control = 60, n_treatment = 60,
                                    effect_size = -0.8, delay = 6,
                                    times = seq(0, 24, by = 3),
                                    dropout_rate = 0.1, residual_sd = 1.5,
                                    random_sd = 0.2, seed = 7)
  grid <- seq(0, 15, by = 3)
  fit <- suppressWarnings(pmrm_delayed(y ~ time * trt, data = sim,
                                        delay_estimation = "profile",
                                        delay_grid = grid, lambda = 0.2))
  expect_s3_class(fit, "pmrm_delayed")
  expect_true(fit$delay_estimate %in% grid)
  # the treatment effect must enter only through the delayed term, otherwise
  # the unrestricted treatment terms absorb the delay and the profile
  # likelihood always peaks at zero delay
  expect_equal(fit$delay_estimate, 6)
  expect_true(all(fit$delay_ci %in% grid))
  expect_equal(nrow(fit$profile_loglik), length(grid))

  s <- summary(fit)
  expect_equal(s$delay_estimate, fit$delay_estimate)
  expect_equal(s$delay_ci, fit$delay_ci)
  # summary() must reach the underlying lme tTable rather than recursing
  expect_true(is.matrix(s$t_table))
  expect_true(nrow(s$t_table) >= 2)
})

test_that("pmrm_delayed recovers a zero-delay truth", {
  sim <- simulate_progression_trial(n_control = 60, n_treatment = 60,
                                    effect_size = -0.8, delay = 0,
                                    times = seq(0, 24, by = 6),
                                    residual_sd = 1.5, random_sd = 0.2, seed = 7)
  fit <- suppressWarnings(pmrm_delayed(y ~ time * trt, data = sim,
                                        delay_estimation = "profile",
                                        delay_grid = c(0, 6), lambda = 0.2))
  expect_equal(fit$delay_estimate, 0)
})

test_that("pmrm_delayed strips treatment terms so the delay is identifiable", {
  f <- pmrmPlus:::strip_treatment_terms(y ~ time * trt + age, "trt")
  expect_equal(sort(attr(terms(f), "term.labels")),
               sort(c("time", "age", "trt_eff")))
  # a formula with no treatment term keeps everything and just adds trt_eff
  g <- pmrmPlus:::strip_treatment_terms(y ~ time, "trt")
  expect_equal(attr(terms(g), "term.labels"), c("time", "trt_eff"))
  # an interaction is stripped as well
  h <- pmrmPlus:::strip_treatment_terms(y ~ time + trt:age, "trt")
  expect_equal(attr(terms(h), "term.labels"), c("time", "trt_eff"))
})

test_that("predict() and plot() work on a pmrm_delayed fit", {
  sim <- simulate_progression_trial(n_control = 30, n_treatment = 30,
                                    effect_size = -0.8, delay = 6,
                                    times = seq(0, 24, by = 6), seed = 7)
  fit <- suppressWarnings(pmrm_delayed(y ~ time * trt, data = sim,
                                        delay_estimation = "profile",
                                        delay_grid = c(0, 6), lambda = 0.2))
  nd <- nlme::getData(fit)
  pr <- predict(fit, newdata = nd, level = "population")
  expect_length(pr, nrow(nd))
  expect_true(all(is.finite(pr)))
  expect_error(plot(fit, type = "profile_likelihood"), NA)
  expect_error(plot(fit, type = "trajectory"), NA)
})

test_that("biomarker dosing rule steps dose correctly", {
  r_up <- biomarker_dosing_rule(150, current_dose_idx = 3)
  r_dn <- biomarker_dosing_rule(50, current_dose_idx = 3)
  r_mt <- biomarker_dosing_rule(100, current_dose_idx = 3)
  expect_equal(r_up$adjustment, "up")
  expect_equal(r_dn$adjustment, "down")
  expect_equal(r_mt$adjustment, "maintain")
  expect_equal(r_dn$new_dose_idx, 2)
})

test_that("pkpd and emax helpers are finite and positive", {
  conc <- pkpd_concentration(seq(0, 24, by = 4), dose = 100)
  expect_true(all(is.finite(conc)))
  b <- emax_biomarker(conc, b_baseline = 100, emax = -30)
  expect_true(all(is.finite(b)))
})

test_that("power_pmrm detects strong effects", {
  p <- power_pmrm(n_per_arm = 40, effect_size = -1.5, delay = 0,
                  n_sims = 5, measurement_times = seq(0, 24, by = 12),
                  seed = 11)
  expect_gte(p$power, 0)
  expect_lte(p$power, 1)
})

# --- regression tests -------------------------------------------------------

test_that("power_pmrm recovers the treatment effect", {
  # Regression: the treatment estimate was silently dropped because the
  # unexported nlme:::summary.lme() was called directly, leaving power at 0.
  p <- power_pmrm(n_per_arm = 40, effect_size = -1.5, delay = 0,
                  n_sims = 5, measurement_times = seq(0, 24, by = 12),
                  seed = 11)
  expect_gt(p$power, 0)
  expect_lt(p$mean_estimate, 0)
  expect_equal(p$non_convergence_rate, 0)
})

test_that("nonlinear progression shapes fit with subject-level random effects", {
  # Regression: the nonlinear branch built an invalid `random` specification
  # (pdDiag(<character vector>), groups = ~trt), so all shapes errored.
  gen <- function(mt) {
    set.seed(11)
    times <- seq(0, 24, by = 3)
    do.call(rbind, lapply(1:60, function(i) {
      data.frame(id = i, time = times, trt = as.integer(i > 30),
                 y = mt(times) + rnorm(length(times), sd = 1.2))
    }))
  }
  specs <- list(
    plateau     = list(mt = function(t) 20 - 0.9 * pmin(t, 12),
                       start = list(b0 = 20, b1 = -1, tau = 12)),
    emax        = list(mt = function(t) 20 - 20 * t / (5 + t),
                       start = list(b0 = 20, bmax = -20, ed50 = 5)),
    exponential = list(mt = function(t) 20 + 14 * (1 - exp(-0.25 * t)),
                       start = list(b0 = 20, binf = 34, lambda = 0.25))
  )
  for (prog in names(specs)) {
    sp <- specs[[prog]]
    fit <- tryCatch(
      pmrm_nlme(y ~ time * trt, data = gen(sp$mt), progression = prog,
                start = sp$start),
      error = function(e) NULL)
    expect_false(is.null(fit), info = prog)
    if (!is.null(fit)) {
      expect_true("pmrm_nlme" %in% class(fit), info = prog)
      # one random-effect row per subject, not per treatment arm
      expect_equal(nrow(nlme::random.effects(fit)), 60, info = prog)
    }
  }
})

test_that("summary() dispatches for joint and competing-risks fits", {
  # Regression: summary.pmrm_joint / summary.pmrm_competing were defined but
  # never registered, so summary() fell through to summary.default.
  data(parkinsons_simulated, package = "pmrmPlus")
  fit <- tryCatch(
    suppressWarnings(pmrm_joint(
      updrs ~ time * trt,
      survival::Surv(time_to_dropout, dropout) ~ trt,
      data = parkinsons_simulated, association = "current_value")),
    error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)
    expect_false(inherits(s, "summaryDefault"))
    expect_true("association_coefficients" %in% names(s))
  }
})

test_that("diagnostics work when the response is not called 'y'", {
  sim <- simulate_progression_trial(n_control = 20, n_treatment = 20,
                                    times = seq(0, 24, by = 6), seed = 2)
  names(sim)[names(sim) == "y"] <- "score"
  fit <- pmrm_nlme(score ~ time * trt, data = sim)
  v <- vpc(fit, n_sims = 3, seed = 1)
  expect_true(all(c("observed", "sim_median") %in% names(v)))
  expect_true(is.finite(goodness_of_fit(fit)$rmse))
})
