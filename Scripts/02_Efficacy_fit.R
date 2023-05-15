# Efficacy_fit
## Packages ----
library(tidyverse)
library(brms)
library(tidybayes)
library(writexl)

if(R.version$os == "linux-gnu") {
  setwd("~/Documents/Vaccine_Efficacy/Vaccine_Efficacy_LMICs")
  cmdstanr::set_cmdstan_path("/opt/cmdstanr/cmdstan-2.29.2")
  mybak = "cmdstanr"
} else {
  mybak = "rstan"
}
## Function to create summary statistics from posteriors
MakeSummaries <- function(drws){
  ## Get list of drug classes
  dcs <- vacc1$who_atc_lbl %>% unique() %>% sort() %>% str_replace_all("\\s", ".")
  # Get effect estimates for low income at drug class level
  dc_lvl_lm <- map(dcs, ~ {
    # take effect estimate for high income (ie the intercept) and low income and drug classes with these estimates in order to add them together to get the
    # drug class level effect estimate for low income trials
    lm <- drws %>% select(b_Intercept, b_income_cat_finlm, contains(.x))
    lm <- lm [ , str_detect(names(lm), "Intercept|income_cat_finlm")]
    lm %>% 
      as.matrix() %>%
      rowSums() %>%
      mean_qi
  })
  # Get effect estimates for high income at drug class level
  dc_lvl_hi <- map(dcs, ~ {
    # take effect estimate for high income (ie the intercept) and drug classes with these estimates
    hi <- drws %>% select(b_Intercept, contains(.x))
    hi <- hi [ , str_detect(names(hi), "Intercept")]
    hi %>% 
      as.matrix() %>%
      rowSums() %>%
      mean_qi
  }
  )
  # Get interaction effect estimates for low income income at drug class level
  dc_lvl_nter <- map(dcs, ~ {
    # take effect estimate for low income (ie the interaction) and drug classes with these estimates
    hi <- drws %>% select(b_income_cat_finlm, contains(.x))
    hi <- hi [ , str_detect(names(hi), "income_cat_finlm")]
    hi %>% 
      as.matrix() %>%
      rowSums() %>%
      mean_qi
  }
  )
  
  # Get interaction effect estimates for low income income at whole analysis level
  overall_nter <- drws %>% select(b_income_cat_finlm) %>% 
    as.matrix() %>%
    rowSums() %>%
    mean_qi()
  overall_lm <- drws %>% select(b_income_cat_finlm, b_Intercept) %>% 
    as.matrix() %>%
    rowSums() %>%
    mean_qi()
  overall_hi <- drws %>% select(b_Intercept) %>% 
    as.matrix() %>%
    rowSums() %>%
    mean_qi()
  overall <- bind_rows(inter = overall_nter,
                       lm = overall_lm,
                       hi = overall_hi,
                       .id = "income") %>% 
    mutate(drug_class = "overall")
  
  ## Name drug classes
  names(dc_lvl_lm) <- dcs
  names(dc_lvl_hi) <- dcs
  names(dc_lvl_nter) <- dcs
  ## Bind into single dataframes
  dc_lvl_lm <- bind_rows(dc_lvl_lm, .id = "drug_class")
  dc_lvl_hi <- bind_rows(dc_lvl_hi, .id = "drug_class")
  dc_lvl_nter  <- bind_rows(dc_lvl_nter, .id = "drug_class")
  dc_lvl <- bind_rows(hi = dc_lvl_hi,
                      lm = dc_lvl_lm,
                      inter = dc_lvl_nter,
                      .id = "income")
  ## Add overall estimates onto drug class to compare
  dc_lvl <- bind_rows(dc_lvl,
                      overall)
  
  # Rounding function which avoids truncating decimals
  MyRound <- function(x) {
    formatC(round(x,2), digits = 2, format = "f")
  }
  # Produce odds ratio table
  or <- dc_lvl %>% 
    mutate(across(c(y, ymin, ymax), ~ MyRound(exp(.x))),
           res = paste0(y, " (", ymin, " to ", ymax, ")")) %>% 
    select(income, drug_class, res) %>% 
    spread(income, res)
  # produce log-odds ratio (original scale model fit on) table
  lor <- dc_lvl %>% 
    mutate(across(c(y, ymin, ymax), ~ MyRound(.x)),
           res = paste0(y, " (", ymin, " to ", ymax, ")")) %>% 
    select(income, drug_class, res) %>% 
    spread(income, res)
  # produce log-odds ratio (original scale model fit on) variation summaries
  smries <- list(
    sd_class_hi   = c(m = mean(drws$sd_who_atc_lbl__Intercept), s = sd(drws$sd_who_atc_lbl__Intercept)) %>% MyRound(),
    sd_class_nter = c(m = mean(drws$sd_who_atc_lbl__income_cat_finlm), s = sd(drws$sd_who_atc_lbl__income_cat_finlm)) %>% MyRound(),
    sd_trial_hi   = c(m = mean(drws$sd_study_id__Intercept), s = sd(drws$sd_study_id__Intercept)) %>% MyRound(),
    sd_trial_nter = c(m = mean(drws$sd_study_id__income_cat_finlm), s = sd(drws$sd_study_id__income_cat_finlm)) %>% MyRound()
  ) %>% bind_rows(.id = "smry")
  
  list(or = or, lor = lor, lor_smries = smries, data_for_plots = dc_lvl)
}

## read in data ----
vacc1 <- readRDS("Outputs/vacc_cleaned_data.Rds")

## Standardise age within each vaccine category ----
vacc1 <- vacc1 %>% 
  group_by(who_atc_lbl) %>% 
  mutate(age_m = mean(age_imputed_years),
         age_s = sd(age_imputed_years)) %>% 
  ungroup()
## two trials in Hep E which have the same age, so for thse set both ages to zero
vacc1 <- vacc1 %>% 
  mutate(age_std_m = if_else(who_atc_lbl == "J07BC01 Hepatitis vaccines", 0, (age_imputed_years - age_m)/age_s),
         age_std_s = if_else(who_atc_lbl == "J07BC01 Hepatitis vaccines", 0.001, sd_imputed_years/age_s))

## Run Bayesian hierarchical model random effect for trial ----
## need to run for each additional covariate
my_priors <- c(
  prior(student_t(3, 0, 2.5), class = Intercept), # population estimate treatmetn effect
  prior(student_t(3, 0, 2.5), class = b), # population estimate treatment  interactions
  prior(student_t(3, 0, 2.5), class = sd))    # variation of between trial, class and condition prior

  ## delete the model if want to re-run it
  mod_re_trial <- brm(ce|se(se) ~ income_cat_fin + (income_cat_fin|who_atc_lbl) 
                      + (income_cat_fin|study_id) 
                      +  me(age_std_m, age_std_s),
                      data = vacc1, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15), 
                      prior = my_priors, 
                      cores = 4,
                      iter = 3000,
                      backend = mybak)
  saveRDS(mod_re_trial, "Outputs/re_model1.Rds")
mod_re_trial <- readRDS("Outputs/re_model1.Rds")
summary(mod_re_trial)
main_results <- MakeSummaries(as_draws_df(mod_re_trial))
rm(mod_re_trial)
gc()

#controltype_harmonised
controltype_harmonised <- brm( ce|se(se) ~ income_cat_fin + controltype_harmonised 
                               + (income_cat_fin+ controltype_harmonised|who_atc_lbl) 
                               + (income_cat_fin+ controltype_harmonised|study_id)
                               + me(age_std_m, age_std_s),
                               data = vacc1, 
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               prior = my_priors, 
                               iter = 2500,
                               cores = 4,
                               backend = mybak)
saveRDS(controltype_harmonised, "Outputs/re_model_control.Rds")
controltype_harmonised_smry <- MakeSummaries(as_draws_df(controltype_harmonised))
controltype_harmonised_smry$or 
write.csv(controltype_harmonised_smry$or,"Outputs/controltype_harmonised_smry_odds.csv", row.names=F)

## compare. The most relevant comparisons are the interaction term. 
controltype_harmonised_smry$lor 
main_results$lor
main_results$or
write.csv(main_results$or,"Outputs/main_results_odds.csv", row.names=F)

#export csv data file
write.csv(controltype_harmonised_smry$lor,"Outputs/controltype_harmonised_smry.csv", row.names=F)
write.csv(main_results$lor,"Outputs/main_results_logodds.csv", row.names=F)

#year_harmonised
## standardise year
vacc1 <- vacc1 %>% 
  mutate(year_hrm_stnd = (year_harmonised - 2003)/8)
year_harmonised <- brm( ce|se(se) ~ income_cat_fin + year_hrm_stnd 
                        + (income_cat_fin+ year_hrm_stnd|who_atc_lbl)
                        + (income_cat_fin+ year_hrm_stnd|study_id)
                        + me(age_std_m, age_std_s),
                        data = vacc1, 
                        control = list(adapt_delta = 0.99, max_treedepth = 15), 
                        prior = my_priors, 
                        iter = 2500,
                        cores = 4,
                        backend = mybak)
saveRDS(year_harmonised, "Scratch_data/year_harmonised.Rds")
year_harmonised_smry <- MakeSummaries(as_draws_df(year_harmonised))
year_harmonised_smry$lor
main_results$lor
write.csv(year_harmonised_smry$lor,"Outputs/year_harmonised_smry.csv", row.names=F)
year_harmonised_smry$or
write.csv(year_harmonised_smry$or,"Outputs/year_harmonised_smry_odds.csv", row.names=F)

#phase_harmonised
phase_harmonised <- brm( ce|se(se) ~ income_cat_fin + phase_harmonised 
                         + (income_cat_fin+ phase_harmonised|who_atc_lbl) 
                         + (income_cat_fin+ phase_harmonised|study_id)
                         + me(age_std_m, age_std_s),
                         data = vacc1, 
                         control = list(adapt_delta = 0.999, max_treedepth = 15), 
                         prior = my_priors, 
                         iter = 2500,
                         cores = 4,
                         backend = mybak)
saveRDS(phase_harmonised, "Scratch_data/phase_harmonised.Rds")
phase_harmonised_smry <- MakeSummaries(as_draws_df(phase_harmonised))
phase_harmonised_smry$lor
main_results$lor
write.csv(phase_harmonised_smry$lor,"Outputs/phase_harmonised_smry.csv", row.names=F)
phase_harmonised_smry$or
write.csv(phase_harmonised_smry$or,"Outputs/phase_harmonised_smry_odds.csv", row.names=F)

#blinding_harmonised
table(vacc1$blinding_harmonised)
vacc1$blinding_harmonised[vacc1$blinding_harmonised=="double-blind"]="double-blind"
vacc1$blinding_harmonised[vacc1$blinding_harmonised %in% c("open-label" , "single-blind") ]="others"
table(vacc1$blinding_harmonised)

blinding_harmonised <- brm( ce|se(se) ~ income_cat_fin + blinding_harmonised 
                            + (income_cat_fin+ blinding_harmonised|who_atc_lbl) 
                            + (income_cat_fin+ blinding_harmonised|study_id)
                            + me(age_std_m, age_std_s),
                            data = vacc1, 
                            control = list(adapt_delta = 0.99, max_treedepth = 15), 
                            prior = my_priors, 
                            iter = 2500,
                            cores = 4,
                            backend = mybak)
saveRDS(blinding_harmonised, "Scratch_data/blinding_harmonised.Rds")
blinding_harmonised_smry <- MakeSummaries(as_draws_df(blinding_harmonised))
blinding_harmonised_smry$lor
main_results$lor
write.csv(blinding_harmonised_smry$lor,"Outputs/blinding_harmonised_smry.csv", row.names=F)
blinding_harmonised_smry$or
write.csv(blinding_harmonised_smry$or,"Outputs/blinding_harmonised_smry_odds.csv", row.names=F)

#administration_harmonised
administration_harmonised <- brm( ce|se(se) ~ income_cat_fin + administration_harmonised 
                                  + (income_cat_fin+ administration_harmonised|who_atc_lbl) 
                                  + (income_cat_fin+ administration_harmonised|study_id)
                                  + me(age_std_m, age_std_s),
                                  data = vacc1, 
                                  control = list(adapt_delta = 0.99, max_treedepth = 15), 
                                  prior = my_priors, 
                                  iter = 2500,
                                  cores = 4,
                                  backend = mybak)
saveRDS(administration_harmonised, "Scratch_data/administration_harmonised.Rds")
administration_harmonised_smry <- MakeSummaries(as_draws_df(administration_harmonised))
administration_harmonised_smry$lor
main_results$lor
write.csv(administration_harmonised_smry$lor,"Outputs/administration_harmonised_smry.csv", row.names=F)
administration_harmonised_smry$or
write.csv(administration_harmonised_smry$or,"Outputs/administration_harmonised_smry_odds.csv", row.names=F)

#analysis_harmonised
analysis_harmonised <- brm( ce|se(se) ~ income_cat_fin + analysis_harmonised 
                            + (income_cat_fin+ analysis_harmonised|who_atc_lbl) 
                            + (income_cat_fin+ analysis_harmonised|study_id)
                            + me(age_std_m, age_std_s),
                            data = vacc1, 
                            control = list(adapt_delta = 0.99, max_treedepth = 15), 
                            prior = my_priors, 
                            iter = 2500,
                            cores = 4,
                            backend = mybak)
saveRDS(analysis_harmonised, "Scratch_data/analysis_harmonised.Rds")
analysis_harmonised_smry <- MakeSummaries(as_draws_df(analysis_harmonised))
analysis_harmonised_smry$lor
main_results$lor
write.csv(analysis_harmonised_smry$lor,"Outputs/analysis_harmonised_smry.csv", row.names=F)
analysis_harmonised_smry$or
write.csv(analysis_harmonised_smry$or,"Outputs/analysis_harmonised_smry_odds.csv", row.names=F)

#vactype_harmonised
vactype_harmonised <- brm( ce|se(se) ~ income_cat_fin + vactype_harmonised 
                           + (income_cat_fin+ vactype_harmonised|who_atc_lbl) 
                           + (income_cat_fin+ vactype_harmonised|study_id)
                           + me(age_std_m, age_std_s),
                           data = vacc1, 
                           control = list(adapt_delta = 0.99, max_treedepth = 15), 
                           prior = my_priors, 
                           iter = 2500,
                           cores = 4,
                           backend = mybak)
saveRDS(vactype_harmonised, "Scratch_data/vactype_harmonised.Rds")
vactype_harmonised_smry <- MakeSummaries(as_draws_df(vactype_harmonised))
vactype_harmonised_smry$lor
main_results$lor
write.csv(vactype_harmonised_smry$lor,"Outputs/vactype_harmonised_smry.csv", row.names=F)
vactype_harmonised_smry$or
write.csv(vactype_harmonised_smry$or,"Outputs/vactype_harmonised_smry_odds.csv", row.names=F)

# #age_imputed
age_imputed_years <- brm( ce|se(se) ~ income_cat_fin + age_std_m
                        + (income_cat_fin+ age_std_m|who_atc_lbl)
                        + (income_cat_fin+ age_std_m|study_id),
                        data = vacc1,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        prior = my_priors,
                        iter = 2500,
                        cores = 4,
                        backend = mybak)
saveRDS(age_imputed_years, "Scratch_data/age_imputed_years.Rds")
age_imputed_years_smry <- MakeSummaries(as_draws_df(age_imputed_years))
age_imputed_years_smry$lor
main_results$lor
write.csv(age_imputed_years_smry$lor,"Outputs/age_imputed_years_smry.csv", row.names=F)
age_imputed_years_smry$or
write.csv(age_imputed_years_smry$or,"Outputs/age_imputed_years_smry_odds.csv", row.names=F)

write_xlsx(list( mainresult=main_results$or,
                 controltype=controltype_harmonised_smry$or,
                 year=year_harmonised_smry$or,
                 phase=phase_harmonised_smry$or,
                 blinding=blinding_harmonised_smry$or,
                 administration=administration_harmonised_smry$or,
                 analysis=analysis_harmonised_smry$or,
                 vactype=vactype_harmonised_smry$or,
                 age=age_imputed_years_smry$or),
           "Outputs/covariate_effect_estimate_odds.xlsx")

write_xlsx(list( mainresult=main_results$lor,
                 controltype=controltype_harmonised_smry$lor,
                 year=year_harmonised_smry$lor,
                 phase=phase_harmonised_smry$lor,
                 blinding=blinding_harmonised_smry$lor,
                 administration=administration_harmonised_smry$lor,
                 analysis=analysis_harmonised_smry$lor,
                 vactype=vactype_harmonised_smry$lor,
                 age=age_imputed_years_smry$lor),
           "Outputs/covariate_effect_estimate_logodds.xlsx")
