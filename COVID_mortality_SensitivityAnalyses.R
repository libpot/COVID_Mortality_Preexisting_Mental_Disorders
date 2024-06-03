#############################
### SENSITIVITY ANALYSES ###
############################

# loading libraries
library(tidyverse)
library(broom)
library(survival)
library(data.table)
library(EValue)
library(patchwork)
library(survminer)


# importing data
data_main_analysis <- read_csv("XXXData/data_main_analysis.csv")

# absolute frequencies of negative control outcomes
sensitivity_abs_freq <- data_main_analysis %>%
  filter(!(cohort == "psychotic disorders" & time_period ==  "period 1"),
         !(cohort == "substance use disorders" & time_period ==  "period 1")) %>% 
  pivot_longer(cols = ends_with("_5letpred"),
               names_to = "condition_type", 
               values_to = "condition_binary") %>%  
  mutate(condition_type = str_remove(condition_type, "_5letpred")) %>% 
  group_by(time_period, condition_type, cohort, exposure_definition) %>% 
  nest() %>% 
  mutate(counts_covid_28 = map_dbl(.x = data,
                                   ~ sum(.x$death_covid_28 == 1 & .x$condition_binary == 1)),
         counts_covid_60 = map_dbl(.x = data,
                                   ~ sum(.x$death_covid_60 == 1 & .x$condition_binary == 1)),
         counts_allcause_28 = map_dbl(.x = data,
                                      ~ sum(.x$death_allcause_28 == 1 & .x$condition_binary == 1)),
         counts_allcause_60 = map_dbl(.x = data,
                                      ~ sum(.x$death_allcause_60 == 1 & .x$condition_binary == 1))) %>% 
  select(cohort, exposure_definition, time_period, condition_type, starts_with("counts")) %>%
  pivot_wider(names_from = "condition_type", 
              values_from = c("counts_covid_28", "counts_covid_60", "counts_allcause_28", "counts_allcause_60"))

write_csv(sensitivity_abs_freq, "XXXResults/Tables/Supplement/sensitivity_abs_freq.csv")



### NEGATIVE CONTROL ANALYSIS ###

# pre-vaccination period
sensitivity_pre_vacc <- data_main_analysis %>%
  filter(time_period %in% c("period 1", "period 2")) %>% # pre-vaccination period
  filter(!(cohort == "psychotic disorders" & time_period ==  "period 1"),
         !(cohort == "substance use disorders" & time_period ==  "period 1")) %>% 
  pivot_longer(cols = ends_with("_5letpred"),
               names_to = "condition_type", 
               values_to = "condition_binary") %>%  
  mutate(condition_type = str_remove(condition_type, "_5letpred")) %>% 
  group_by(time_period, condition_type, cohort, exposure_definition) %>% 
  nest() %>%   
  mutate(covid_28 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ condition_binary + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                data = .x)),
         covid_60 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ condition_binary + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                data = .x)),
         all_cause_28 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ condition_binary + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         all_cause_60 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ condition_binary + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)))

# vaccination_period
sensitivity_vacc <- data_main_analysis %>%
  filter(time_period %in% c("period 3", "period 4", "period 5")) %>% # vaccination period
  pivot_longer(cols = ends_with("_5letpred"),
               names_to = "condition_type",
               values_to = "condition_binary") %>%  
  mutate(condition_type = str_remove(condition_type, "_5letpred")) %>% 
  group_by(time_period, condition_type, cohort, exposure_definition) %>% 
  nest() %>%   
  mutate(covid_28 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ condition_binary + vaccination + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                data = .x)),
         covid_60 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ condition_binary + vaccination + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                data = .x)),
         all_cause_28 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ condition_binary + vaccination + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         all_cause_60 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ condition_binary + vaccination + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)))

sensitivity_results <- bind_rows(sensitivity_pre_vacc, sensitivity_vacc) %>% 
  pivot_longer(cols = covid_28:all_cause_60,
               names_to = "mortality") %>%
  mutate(model_estimates = map(.x = value,
                               ~ .x %>% tidy(conf.int = TRUE, 
                                             exponentiate = TRUE) %>%
                                 filter(term == "condition_binary") %>% 
                                 select(estimate, conf.low, conf.high))) %>% 
  select(condition_type, cohort, time_period, mortality, model_estimates) %>%
  unnest() %>% 
  ungroup() %>% 
  na.omit() %>% 
  filter(!is.infinite(conf.high)) %>% 
  group_by(cohort, time_period) %>%
  summarise(n_test = n(),
            non_null_test = n_test - sum(conf.low < 1 & conf.high > 1),
            prop_non_null_test = formatC(round(non_null_test/n_test*100, 2), format = "f", digits = 2),
            mean_estimate = formatC(round(mean(estimate),2), format = "f", digits = 2))

write_csv(sensitivity_results, "XXXResults/Tables/Supplement/negative_control.csv")


# E-values 
E_values <- models_results %>%
  select(-results) %>% 
  ungroup() %>% 
  filter(!is.na(estimate),
         adjusted == "adjusted") %>% 
  mutate(E_value = pmap_dbl(across(c(estimate, conf.low, conf.high)), 
                            ~ as_tibble(evalues.HR(..1, ..2, ..3, rare = TRUE))$point[[2]]),
         E_value = ifelse(conf.low < 1 & conf.high > 1, NA, E_value)) %>% 
  select(cohort, time_period, mortality, follow_up, exposure_definition, E_value) %>%
  mutate(E_value = formatC(round(E_value, 2), format = "f", digits = 2)) %>% 
  pivot_wider(names_from = c("follow_up", "exposure_definition"),
              values_from = "E_value") %>% 
  arrange(mortality) %>% 
  select(cohort, time_period, mortality, 
        '28_diagnosed', '28_diagnosed and treated', '60_diagnosed', '60_diagnosed and treated')

write_csv(E_values, "XXXResults/Tables/Supplement/E_values.csv")

