###########################
### DATA PRE-PROCESSING ###
###########################

# loading libraries
library(tidyverse)
library(broom)
library(data.table)
library(fastDummies)

# importing data
data_main_analysis <- read_csv("XXXData/data_all_cohorts_all_outcomes.csv")

# re-coding the data
data_main_analysis <- data_main_analysis %>% 
  mutate(exposure = ifelse(str_detect(cohort, "KS"), "control", "exposed"),
         kohorta = cohort,
         cohort = str_remove(cohort, " KS"),
         time_to_death_allcause_28 = time_to_umrti_28,
         time_to_death_allcause_60 = time_to_umrti_60,
         time_to_death_covid_28 = time_to_umrti_covid_28,
         time_to_death_covid_60 = time_to_umrti_covid_60,
         death_covid_28 = umrti_covid_28,
         death_covid_60 = umrti_covid_60,
         death_allcause_28 = umrti_28,
         death_allcause_60 = umrti_60, 
         age_group = vek_skupina,
         sex = pohlavi,
         vaccination = ockovani,
         matching = propojeni_cohort_KS,
         number_matched = pocet_napojenych_KS,
         infection_year = pozitivita_rok,
         infection_month = pozitivita_mesic,
         ID_patient = ID_pacienta,
         region_residence = Bydliste_Kraj,
         DCCI_cat = case_when(DCCI == 0 ~ "0",
                              DCCI == 1 ~ "1",
                              DCCI == 2 ~ "2",
                              DCCI == 3 ~ "3",
                              DCCI >= 4 ~ "4+"),
         number_outpatient = case_when(cohort %in% c("F1", "F1 leky") ~ pocet_kontaktu_bezF1,
                                       cohort %in% c("F2", "F2 leky") ~ pocet_kontaktu_bezF2,
                                       cohort %in% c("F3", "F3 leky") ~ pocet_kontaktu_bezF3,
                                       cohort %in% c("F4", "F4 leky") ~ pocet_kontaktu_bezF4,
                                       cohort %in% c("F1F2F3F4", "F1F2F3F4 leky") ~ pocet_kontaktu_bezF1F2F3F4),
         number_inpatient = case_when(cohort %in% c("F1", "F1 leky") ~ pocet_hosp_bezF1,
                                      cohort %in% c("F2", "F2 leky") ~ pocet_hosp_bezF2,
                                      cohort %in% c("F3", "F3 leky") ~ pocet_hosp_bezF3,
                                      cohort %in% c("F4", "F4 leky") ~ pocet_hosp_bezF4,
                                      cohort %in% c("F1F2F3F4", "F1F2F3F4 leky") ~ pocet_hosp_bezF1F2F3F4)) %>% 
  mutate(exposure_definition = ifelse(str_detect(cohort, "leky"), "diagnosed and treated", "diagnosed"),
         exp_def = ifelse(exposure_definition == "diagnosed", "dg", "dg & tr"),
         cohort_ICD = factor(str_remove(cohort, " leky"),
                             levels = c("F1F2F3F4", "F1", "F2", "F3", "F4"),
                             labels = c("Fx", "F1", "F2", "F3", "F4"),
                             ordered = TRUE),
         cohort = factor(str_remove(cohort, " leky"),
                         levels = c("F1F2F3F4", "F1", "F2", "F3", "F4"),
                         labels = c("any mental disorder", "substance use disorders", "psychotic disorders", "affective disorders", "anxiety disorders"),
                         ordered = TRUE),
         time_period = factor(as.character(paste0("period ", obdobi)),
                              levels = c("period 1", "period 2", "period 3", "period 4", "period 5"),
                              labels = c("period 1", "period 2", "period 3", "period 4", "period 5"),
                              ordered = TRUE)) %>% 
  filter(age_group != "0-4",
         age_group != "5-9") # removing children 9 years old and younger

# N (%) of unmatched individuals
N_unmatched <- data_main_analysis %>% 
  filter(exposure == "exposed") %>% 
  group_by(cohort, time_period, exposure_definition) %>% 
  summarise(N_overall = n(),
            N_unmatched = paste0(sum(number_matched == 0), " (", formatC(round(sum(number_matched == 0)/n()*100, 2), format = "f", digits = 2), ")")) %>% 
  pivot_wider(names_from = c("exposure_definition"), 
              values_from = c("N_overall", "N_unmatched")) %>% 
  select(cohort, time_period, N_overall_diagnosed, N_unmatched_diagnosed, `N_overall_diagnosed and treated`, `N_unmatched_diagnosed and treated`)

write_csv(N_unmatched, "XXXResults/Tables/Supplement/N_unmatched.csv")

# sample characteristics: matched vs. unmatched individuals
descriptives_unmatched <- data_main_analysis %>% 
  filter(exposure == "exposed") %>% 
  mutate(matching_success = ifelse(number_matched == 0, "unmatched", "matched")) %>% 
  group_by(cohort_ICD, matching_success, exposure_definition, time_period) %>% 
  nest() %>% 
  mutate(N = map_chr(.x = data,
                     ~ as.character(length(.x$ID_patient))),
         female = map_chr(.x = data,
                          ~ paste0(sum(.x$sex == "Z"), " (", formatC(round(sum(.x$sex == "Z")/length(.x$sex)*100,2), format = "f", digits = 2), ")")),
         not_vaccinated = map_chr(.x = data,
                                  ~ paste0(sum(.x$vaccination == "neockovan"), " (", formatC(round(sum(.x$vaccination == "neockovan")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         first_only = map_chr(.x = data,
                              ~ paste0(sum(.x$vaccination == "pouze prvni davka"), " (", formatC(round(sum(.x$vaccination == "pouze prvni davka")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         first_extra = map_chr(.x = data,
                               ~ paste0(sum(.x$vaccination == "prvni extra davka"), " (", formatC(round(sum(.x$vaccination == "prvni extra davka")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         complete = map_chr(.x = data,
                            ~ paste0(sum(.x$vaccination == "ukoncene ockovani"), " (", formatC(round(sum(.x$vaccination == "ukoncene ockovani")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         across(.cols = c(first_only, first_extra, complete, not_vaccinated),
                ~ ifelse(time_period %in% c("1", "2"), NA, .x)),
         month_infection = map_chr(.x = data, 
                                   ~ paste0(median(.x$infection_month), " (", IQR(.x$infection_month), ")")),
         year_infection = map_chr(.x = data,
                                  ~ paste0(median(.x$infection_year), " (", IQR(.x$infection_year), ")")),
         DCCI = map_chr(.x = data,
                        ~ paste0(round(mean(.x$DCCI), 2), " (", round(sd(.x$DCCI), 2),  ")"))) %>% 
  select(-data) %>% 
  ungroup() %>% 
  pivot_longer(cols = N:DCCI,
               names_to = "variable") %>% 
  pivot_wider(names_from = c("cohort_ICD", "matching_success"), 
              values_from = c("value")) %>%
  mutate(time_period = str_remove(time_period, "period ")) %>% 
  select(time_period, variable, exposure_definition, 
         Fx_unmatched, Fx_matched, F1_unmatched, F1_matched, F2_unmatched, F2_matched, F3_unmatched, F3_matched, F4_unmatched, F4_matched) %>% 
  arrange(exposure_definition, time_period) %>% 
  filter(!((time_period %in% c("1", "2")) & (variable %in% c("not_vaccinated", "first_only", "first_extra", "complete"))))

write_csv(descriptives_unmatched, "XXXResults/Tables/descriptives_unmatched.csv")

# number of macthed controls (0, 1, 2.. 5)
number_matched_controls <- data_main_analysis %>% 
  group_by(cohort, exposure_definition, time_period) %>% 
  count(number_matched) %>% 
  filter(number_matched != -1) %>% 
  pivot_wider(names_from = "number_matched",
              values_from = "n") %>% 
  mutate(total = `1` + `2`+ `3`+ `4` + `5`) %>% 
  mutate(across(c(`1`:`5`),
                ~ paste0(.x, " (", formatC(round(.x/total*100, 2), format = "f", digits = 2), ")"))) 

number_matched_controls_diag <- number_matched_controls %>% 
  filter(exposure_definition == "diagnosed") %>% 
  select(-c(total, exposure_definition))

number_matched_controls_diag_treat <- number_matched_controls %>% 
  filter(exposure_definition == "diagnosed and treated") %>% 
  select(-c(total, exposure_definition))

write_csv(number_matched_controls_diag, "XXXResults/Tables/number_matched_controls_diag.csv")
write_csv(number_matched_controls_diag_treat, "XXXResults/Tables/number_matched_controls_diag_treat.csv")

# excluding unmatched individuals
data_main_analysis <- data_main_analysis %>% 
  filter(number_matched != 0) 

# saving the data file for main analysis
write_csv(data_main_analysis, "XXXData/data_main_analysis.csv")



### QUALITY CONTROL ###
# checking the matching success per matching variable 
map_dfr(c("age_group", "sex", "vaccination", "DCCI_cat", "infection_year", "infection_month"), 
        ~ data_main_analysis %>% 
        group_by(cohort, time_period, exposure_definition, matching) %>% 
          summarise(cond = all(.[[.x]] == first(.[[.x]]))) %>% 
          ungroup() %>% 
          count(cond))

# checking the number of matched individuals within each stratum
data_main_analysis %>% 
  group_by(cohort, time_period, exposure_definition, matching) %>% 
  count() %>% 
  ungroup() %>% 
  summarise(n_matched = n >= 2) %>% # a minimal number of individuals per strata should equal two (one exposed and one unexposed individual)
  count(n_matched)





### DESCRIPTIVES ###

# descriptive statistics for matching variables
cohort_characteristics <- data_main_analysis %>% 
  group_by(cohort_ICD, exposure, exposure_definition, time_period) %>% 
  nest() %>% 
  mutate(N = map_chr(.x = data,
                     ~ as.character(length(.x$ID_patient))),
         female = map_chr(.x = data,
                          ~ paste0(sum(.x$sex == "Z"), " (", formatC(round(sum(.x$sex == "Z")/length(.x$sex)*100,2), format = "f", digits = 2), ")")),
         not_vaccinated = map_chr(.x = data,
                                  ~ paste0(sum(.x$vaccination == "neockovan"), " (", formatC(round(sum(.x$vaccination == "neockovan")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         first_only = map_chr(.x = data,
                              ~ paste0(sum(.x$vaccination == "pouze prvni davka"), " (", formatC(round(sum(.x$vaccination == "pouze prvni davka")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         first_extra = map_chr(.x = data,
                               ~ paste0(sum(.x$vaccination == "prvni extra davka"), " (", formatC(round(sum(.x$vaccination == "prvni extra davka")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         complete = map_chr(.x = data,
                            ~ paste0(sum(.x$vaccination == "ukoncene ockovani"), " (", formatC(round(sum(.x$vaccination == "ukoncene ockovani")/length(.x$vaccination)*100, 2), format = "f", digits = 2), ")")),
         across(.cols = c(first_only, first_extra, complete, not_vaccinated),
                ~ ifelse(time_period %in% c("1", "2"), NA, .x)),
         month_infection = map_chr(.x = data, 
                                   ~ paste0(median(.x$infection_month), " (", IQR(.x$infection_month), ")")),
         year_infection = map_chr(.x = data,
                                  ~ paste0(median(.x$infection_year), " (", IQR(.x$infection_year), ")")),
         DCCI = map_chr(.x = data,
                        ~ paste0(round(mean(.x$DCCI), 2), " (", round(sd(.x$DCCI), 2),  ")"))) %>% 
  select(-data) %>% 
  ungroup() %>% 
  pivot_longer(cols = N:DCCI,
               names_to = "variable") %>% 
  pivot_wider(names_from = c("cohort_ICD", "exposure"), 
              values_from = c("value")) %>%
  mutate(time_period = str_remove(time_period, "period ")) %>% 
  select(time_period, variable, exposure_definition,
         Fx_control, Fx_exposed, F1_control, F1_exposed, F2_control, F2_exposed, F3_control, F3_exposed, F4_control, F4_exposed) %>% 
  arrange(exposure_definition, time_period) %>% 
  filter(!((time_period %in% c("1", "2")) & (variable %in% c("not_vaccinated", "first_only", "first_extra", "complete"))))

write_csv(cohort_characteristics, "XXXResults/Tables/descriptives.csv")


# descriptive statistics for additional covariates
descriptives_contacts <- data_main_analysis %>% 
  select(cohort_ICD, exposure, exposure_definition, time_period,
         number_inpatient, number_outpatient) %>% 
  group_by(cohort_ICD, exposure, exposure_definition, time_period) %>% 
  summarise(across(c("number_inpatient", "number_outpatient"),
                   ~ paste0(formatC(round(mean(.x), 2), format = "f", digits = 2), " (", formatC(round(sd(.x), 2), format = "f", digits = 2),  ")"))) %>% 
  ungroup()

descriptives_covariates <- data_main_analysis %>% 
  select(cohort_ICD, exposure, exposure_definition, time_period,
         region_residence, ends_with("1rokpred")) %>% 
  dummy_cols(select_columns = "region_residence") %>%
  mutate(across(.cols = c(starts_with("region_"), ends_with("1rokpred")),
                ~ as.double(.x))) %>% 
  select(-c("region_residence")) %>% 
  group_by(cohort_ICD, exposure, exposure_definition, time_period) %>% 
  summarise(across(.cols = c(starts_with("region_"), ends_with("1rokpred")),
                   ~ paste0(sum(.x == 1), " (", formatC(round(sum(.x == 1)/length(.x)*100,2), format = "f", digits = 2), ")"))) %>% 
  ungroup() %>%
  merge(descriptives_contacts, by = c("cohort_ICD", "exposure", "exposure_definition", "time_period")) %>% 
  as_tibble() %>% 
  pivot_longer(cols = region_residence_CZ010:number_outpatient,
               names_to = "variable") %>% 
  pivot_wider(names_from = c("cohort_ICD", "exposure"), 
              values_from = c("value")) %>%
  mutate(time_period = str_remove(time_period, "period ")) %>% 
  select(time_period, variable, exposure_definition,
         Fx_control, Fx_exposed, F1_control, F1_exposed, F2_control, F2_exposed, F3_control, F3_exposed, F4_control, F4_exposed) %>% 
  arrange(exposure_definition, time_period) %>% 
  mutate(variable = str_remove(variable, "region_residence_"),
         variable = str_remove(variable, "_1rokpred"),
         variable = case_when(variable == "CZ031" ~	"South-Bohemian region",
                              variable == "CZ064" ~	"South-Moravian region",
                              variable == "CZ041"	~ "Karlovy Vary region",
                              variable == "CZ052"	~ "Hradec Kralove region",
                              variable == "CZ051"	~ "Liberec region",
                              variable == "CZ080"	~ "Moravian-Silesian region",
                              variable == "CZ071"	~ "Olomouc region",
                              variable == "CZ053"	~ "Pardubice region",
                              variable == "CZ032"	~ "Plzen region",
                              variable == "CZ010"	~ "Prague region",
                              variable == "CZ020"	~ "Central-Bohemian region",
                              variable == "CZ042"	~ "Usti region",
                              variable == "CZ063"	~ "Vysocina region",
                              variable == "CZ072" ~	"Zlin region",
                              variable == "CZ099" ~ "living abroad",
                              variable == "C02" ~ "antihypertensives", 
                              variable == "B01AC06_N02BA01_N02BA51" ~ "aspirin", 
                              variable == "C10AA" ~ "statins",
                              variable == "B01" ~ "antithrombotic agents",
                              variable == "M01A" ~ "non-steroidal anti-inflammatory medications",
                              variable == "M05BA_M05BB" ~ "bisphosphonates",
                              variable == "G03A" ~ "oral contraceptives", 
                              variable == "G03C_G03D_G03F" ~ "hormone replacement therapy",
                              variable == "N03" ~ "anticonvulsants",
                              variable == "L01A_L01B_L01C_L01D_L01E_L01X" ~ "cytostatic chemotherapy",
                              variable == "V10" ~ "radiotherapy",
                              variable == "L04" ~ "immunosuppressant medication",
                              variable == "R03AC13_R03AC12" ~ "long-acting beta-agonist" ,
                              variable == "R03DC" ~ "leukotriene receptor antagonists",
                              variable == "R03BA" ~ "inhaled glucocorticoids",
                              variable == "number_outpatient" ~ "number_outpatient",
                              variable == "number_inpatient" ~ "number_inpatient"))

write_csv(descriptives_covariates, "XXXResults/Tables/descriptives_covariates.csv")

