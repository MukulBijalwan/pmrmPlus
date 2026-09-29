# progme

Extended Progression Models for Repeated Measures.

Tools for progressive-disease clinical trials: nonlinear mixed-effects
progression models, delayed treatment onset estimation via profile
likelihood, joint longitudinal-survival models for informative dropout,
competing risks, biomarker-informed adaptive dosing, and simulation-based
power/sample-size calculation.

## Install

```r
# install.packages("remotes")
remotes::install_github("MukulBijalwan/progme")
```

## Quick start

```r
library(progme)

# Simulate a trial with a 6-month delayed treatment effect
sim <- simulate_progression_trial(n_control = 60, n_treatment = 60,
                                  effect_size = -0.5, delay = 6,
                                  times = seq(0, 24, by = 3))

# Estimate the delay by profile likelihood
fit <- pmrm_delayed(y ~ time * trt, data = sim,
                    delay_estimation = "profile",
                    delay_grid = seq(0, 18, by = 3))
print(fit)
plot(fit, type = "profile_likelihood")

# Simulation-based power
p <- power_pmrm(n_per_arm = 100, effect_size = -0.5, delay = 6, n_sims = 100)
p$power
```

## Authors

* Mukul Bijalwan (mukulbijalwan555@gmail.com) --
  ORCID [0009-0001-3040-6912](https://orcid.org/0009-0001-3040-6912)
* Gunjan Aggrwal
* Mukul Jain
