# Data creation script for pmrmPlus simulated datasets
# Run once from the package root:
#   Rscript data-raw/create_data.R

library(pmrmPlus)

set.seed(2024)

# ---- ADNI-style Alzheimer's simulated dataset (linear progression + delay) ----
adni_simulated <- simulate_progression_trial(
  n_control = 120, n_treatment = 120,
  effect_size = -0.6, delay = 6,
  dropout_rate = 0.18,
  times = seq(0, 24, by = 3),
  intercept = 22, slope = 1.2,
  random_sd = 0.4, residual_sd = 2.5
)
# add subject-level event info (dropout time)
ids <- unique(adni_simulated$id)
ev <- data.frame(id = ids,
                 time_to_dropout = runif(length(ids), 12, 30),
                 dropout = sample(0:1, length(ids), TRUE, prob = c(0.7, 0.3)))
adni_simulated <- merge(adni_simulated, ev, by = "id")
names(adni_simulated)[names(adni_simulated) == "y"] <- "adas_cog"
names(adni_simulated)[names(adni_simulated) == "id"] <- "rid"

# ---- Parkinson's simulated dataset with informative dropout ----
parkinsons_simulated <- simulate_joint_trial(
  n_control = 100, n_treatment = 100,
  effect_size = -0.5, alpha = 0.04,
  times = seq(0, 24, by = 3),
  intercept = 25, slope = 1.5
)
names(parkinsons_simulated)[names(parkinsons_simulated) == "y"] <- "updrs"
names(parkinsons_simulated)[names(parkinsons_simulated) == "time_to_dropout"] <-
  "time_to_dropout"

save(adni_simulated, file = "data/adni_simulated.rda", compress = "bzip2")
save(parkinsons_simulated, file = "data/parkinsons_simulated.rda",
     compress = "bzip2")
