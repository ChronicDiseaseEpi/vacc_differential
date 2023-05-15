#04_review_models
library(brms)
library(tidyverse)


if(R.version$os == "linux-gnu") {
  setwd("~/Documents/Vaccine_Efficacy/Vaccine_Efficacy_LMICs")
  cmdstanr::set_cmdstan_path("/opt/cmdstanr/cmdstan-2.29.2")
  mybak = "cmdstanr"
} else {
  mybak = "rstan"
}

## read in models
a <- list.files("Scratch_data/", patt = "Rds")
b <- map(a, ~ readRDS(file = paste0("Scratch_data/", .x)))


b_smry <- map_int(b, ~rstan::get_num_divergent(.x$fit))

for(i in seq_along(b)){
  print(i)
  print(summary(b))
}
