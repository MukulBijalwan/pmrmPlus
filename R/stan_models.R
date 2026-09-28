# Expert-elicited default priors

#' Default priors for Bayesian progression models
#'
#' @param scale Outcome scale factor (e.g., SD of the response) used to
#'   make weakly-informative priors scale-aware.
#' @return Named list of prior family and parameter values for fixed
#'   effects, random-effect SDs, residual SD, association parameters, and
#'   survival shape.
#' @examples
#' default_priors(scale = 5)
#' @export
default_priors <- function(scale = 1) {
  list(
    fixed_effects = list(family = "normal", mean = 0, sd = 10 * scale),
    random_effect_sd = list(family = "half-normal", scale = 2 * scale),
    residual_sd = list(family = "exponential", rate = 1 / scale),
    association_alpha = list(family = "normal", mean = 0, sd = 1),
    treatment_on_hazard = list(family = "normal", mean = 0, sd = 1),
    weibull_shape = list(family = "gamma", shape = 2, rate = 1),
    delay_tau = list(family = "uniform", lower = 0, upper = 24),
    onset_rate_lambda = list(family = "exponential", rate = 0.1)
  )
}

# Stan compilation & fitting wrappers (require rstan; guarded at runtime)

stan_file <- function(name) {
  system.file("stan", name, package = "pmrmPlus", mustWork = FALSE)
}

has_rstan <- function() requireNamespace("rstan", quietly = TRUE)

fit_joint_stan <- function(stan_data, chains = 4, iter = 2000, seed = NULL, ...) {
  if (!has_rstan()) {
    stop("rstan is required for Bayesian joint modeling. ",
         "Install it with install.packages('rstan').")
  }
  code <- paste(readLines(stan_file("joint_longitudinal_survival.stan")),
                collapse = "\n")
  rstan::stan(model_code = code, data = stan_data,
              chains = chains, iter = iter, seed = seed, ...)
}
