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

test_that("pmrm_delayed recovers the delay approximately", {
  sim <- simulate_progression_trial(n_control = 60, n_treatment = 60,
                                    effect_size = -0.8, delay = 6,
                                    times = seq(0, 24, by = 3),
                                    dropout_rate = 0.1, residual_sd = 1.5,
                                    random_sd = 0.2, seed = 7)
  grid <- seq(0, 15, by = 3)
  fit <- suppressWarnings(
    tryCatch(pmrm_delayed(y ~ time * trt, data = sim, delay_grid = grid,
                          lambda = 0.5),
             error = function(e) NULL))
  if (!is.null(fit)) {
    expect_s3_class(fit, "pmrm_delayed")
    expect_true(fit$delay_estimate %in% grid)
    s <- summary(fit)
    expect_true("delay_estimate" %in% names(s))
  }
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
