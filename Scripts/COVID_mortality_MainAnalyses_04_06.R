#####################
### MAIN ANALYSES ###
#####################

# loading libraries
library(tidyverse)
library(broom)
library(survival)
library(data.table)
library(patchwork)
library(survminer)
library(ggpubr)
library(scales)
library(Cairo)
library(extrafont)
library(ggh4x)


# importing data
data_main_analysis <- read_csv("XXXData/data_main_analysis.csv")

# absolute risk
absolute_risk <- data_main_analysis %>% 
  group_by(cohort, time_period, exposure_definition, exposure) %>% 
  summarise(N = length(ID_patient),
            across(c("death_allcause_28", "death_allcause_60", "death_covid_28", "death_covid_60"), 
                   ~ paste0(sum(.x == 1), " (", formatC(round(sum(.x == 1)/length(.x)*100,2), format = "f", digits = 2), ")"))) %>% 
  pivot_wider(names_from = c("exposure"), 
              values_from = c("N", starts_with("death"))) %>% 
  arrange(exposure_definition)

write_csv(absolute_risk, "XXXResults/Tables/Supplement/absolute_risk.csv")



### MODEL ASSUMPTIONS ###

# Schoenfeld residuals
schoenfeld_residuals_pre_vacc <- data_main_analysis %>% 
  filter(time_period %in% c("period 1", "period 2")) %>% 
  group_by(cohort, time_period) %>% 
  nest() %>%   
  mutate(covid_28 = map(.x = data,
                        ~ cox.zph(coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x))),
         covid_60 = map(.x = data,
                        ~ cox.zph(coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x))),
         all_cause_28 = map(.x = data,
                            ~ cox.zph(coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                            + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                            data = .x))),
         all_cause_60 = map(.x = data,
                            ~ cox.zph(coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                            + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                            data = .x))),
         across(ends_with(c("_60", "_28")),
                ~ as.data.frame(.x$table)$p[[1]] < 0.05))
schoenfeld_residuals_pre_vacc

schoenfeld_residuals_vacc <- data_main_analysis %>% 
  filter(time_period %in% c("period 3", "period 4", "period 5")) %>% 
  group_by(cohort, time_period) %>% 
  nest() %>%   
  mutate(covid_28 = map(.x = data,
                        ~ cox.zph(coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x))),
         covid_60 = map(.x = data,
                        ~ cox.zph(coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x))),
         all_cause_28 = map(.x = data,
                            ~ cox.zph(coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                            + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                            data = .x))),
         all_cause_60 = map(.x = data,
                            ~ cox.zph(coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                            + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                            data = .x))),
         across(ends_with(c("_60", "_28")),
                ~ as.data.frame(.x$table)$p[[1]] < 0.05))
schoenfeld_residuals_vacc



### SURVIVAL MODELS ###

# pre-vaccination period
models_pre_vacc <- data_main_analysis %>% 
  filter(time_period %in% c("period 1", "period 2")) %>% 
  group_by(cohort, exposure_definition, time_period) %>% 
  nest() %>%  
  mutate(covid_28 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                data = .x)),
         covid_60 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                data = .x)),
         all_cause_28 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                    data = .x)),
         all_cause_60 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching),
                                    data = .x)),
         covid_28_adj = map(.x = data,
                            ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         covid_60_adj = map(.x = data,
                            ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching)
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         all_cause_28_adj = map(.x = data,
                                ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching)
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x)),
         all_cause_60_adj = map(.x = data,
                                ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + sex + infection_year + infection_month + DCCI + strata(matching)
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x)))

# vaccination period
models_vacc <- data_main_analysis %>% 
  filter(time_period %in% c("period 3", "period 4", "period 5")) %>% 
  group_by(cohort, exposure_definition, time_period) %>% 
  nest() %>%  
  mutate(covid_28 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                data = .x)),
         covid_60 = map(.x = data,
                        ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                data = .x)),
         all_cause_28 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                    data = .x)),
         all_cause_60 = map(.x = data,
                            ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching),
                                    data = .x)),
         covid_28_adj = map(.x = data,
                            ~ coxph(Surv(time_to_death_covid_28, death_covid_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching) 
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         covid_60_adj = map(.x = data,
                            ~ coxph(Surv(time_to_death_covid_60, death_covid_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching)
                                    + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                    data = .x)),
         all_cause_28_adj = map(.x = data,
                                ~ coxph(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching)
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x)),
         all_cause_60_adj = map(.x = data,
                                ~ coxph(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure + age_group + vaccination + sex + infection_year + infection_month + DCCI + strata(matching)
                                        + region_residence + number_outpatient + number_inpatient + C02_1rokpred + B01AC06_N02BA01_N02BA51_1rokpred + C10AA_1rokpred + B01_1rokpred + M01A_1rokpred + M05BA_M05BB_1rokpred + G03A_1rokpred + G03C_G03D_G03F_1rokpred + N03_1rokpred + L01A_L01B_L01C_L01D_L01E_L01X_1rokpred + V10_1rokpred + L04_1rokpred + R03AC13_R03AC12_1rokpred + R03DC_1rokpred + R03BA_1rokpred,
                                        data = .x)))


models_results <- bind_rows(models_pre_vacc, models_vacc) %>% 
  pivot_longer(cols = covid_28:all_cause_60_adj,
               names_to = "mortality") %>%
  mutate(model_estimates = map(.x = value,
                               ~ .x %>% tidy(conf.int = TRUE, 
                                             exponentiate = TRUE) %>%
                                 filter(term == "exposureexposed") %>% 
                                 select(estimate, conf.low, conf.high, p.value))) %>% 
  select(cohort, time_period, mortality, model_estimates) %>% 
  unnest() %>% 
  mutate(follow_up = str_extract(mortality, "28|60"),
         adjusted = factor(ifelse(is.na(str_extract(mortality, "adj")), "raw", "adjusted"), ordered = TRUE),
         mortality = str_extract(mortality, "covid|all_cause"),
         mortality = ifelse(mortality == "covid", "COVID-19 mortality", "All-cause mortality"),
         cohort_time = as.factor(paste0(cohort, ": ", time_period)),
         results = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                          " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                          formatC(round(conf.high, 2), format = "f", digits = 2), ")"),
         results = ifelse(str_detect(results, "Inf"), "NA", results),
         p.value = as.character(ifelse(p.value < 0.001, "<0.001",  formatC(round(p.value, 3), format = "f", digits = 3))),
         across(.cols = c(estimate, conf.low, conf.high, p.value),
                ~ ifelse(results == "NA", NA, .x)),
         order_cohort = as.numeric(cohort),
         order_time_period = as.numeric(time_period)) %>% 
  arrange(order_cohort, order_time_period)

write_csv(model_results, "XXXResults/Tables/model_results.csv")


main_results_table <- models_results %>% 
  select(cohort, time_period, mortality, follow_up, exposure_definition, adjusted, results, p.value) %>% 
  pivot_wider(names_from = c("follow_up", "adjusted", "exposure_definition"),
              values_from = c("results", "p.value")) %>% 
  mutate(time_period = str_remove(time_period, "period ")) %>% 
  arrange(mortality) %>% 
  select(cohort, time_period, mortality, 
         'results_28_raw_diagnosed', 'p.value_28_raw_diagnosed', 'results_28_adjusted_diagnosed', 'p.value_28_adjusted_diagnosed',  
         'results_28_raw_diagnosed and treated', 'p.value_28_raw_diagnosed and treated', 'results_28_adjusted_diagnosed and treated', 'p.value_28_adjusted_diagnosed and treated',
         'results_60_raw_diagnosed', 'p.value_60_raw_diagnosed', 'results_60_adjusted_diagnosed', 'p.value_60_adjusted_diagnosed', 
         'results_60_raw_diagnosed and treated', 'p.value_60_raw_diagnosed and treated', 'results_60_adjusted_diagnosed and treated',  'p.value_60_adjusted_diagnosed and treated')

write_csv(main_results_table, "XXXResults/Tables/main_results.csv")




### PLOTTING ###

my_colors2   <- c("up to 28 days" ="#0072B2", 
                  "up to 60 days" = "#D55E00")   

font_import()
y
loadfonts(device="win")
windowsFonts(Times = windowsFont("Times New Roman"))

# All-cause mortality
plot_allcause <- models_results %>% 
  mutate(follow_up = ifelse(follow_up == "28", "up to 28 days", "up to 60 days"),
         time_period = str_replace(time_period, "period", "epoch")) %>% 
  filter(adjusted == "adjusted",
         mortality == "All-cause mortality",
         cohort_time != "psychotic disorders: period 1",
         cohort_time != "substance use disorders: period 1") %>% 
  arrange(order_cohort, order_time_period) %>% 
  ggplot(aes(x = time_period,
             y = estimate,
             color = follow_up)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                linewidth = 1,
                width = 0.25,
                position = position_dodge(width = 0.4)) +
  scale_y_continuous(trans = "log10", 
                     limit = c(0.4, 3.4), 
                     breaks = c(0.5, 1, 2)) +
  geom_hline(yintercept = 1, 
             linetype = 2,
             linewidth = 0.5) +
  labs(x = "",
       y = "adjusted hazard ratio (95% confidence interval)",
       color = "Follow-up period",
       title = "Relative risk of all-cause mortality following SARS-CoV-2 infection in people with pre-existing mental disorders") +
  theme_minimal() +
  facet_nested_wrap(~ cohort + exposure_definition,
                           scales = "free",
                           nrow = 5,
                           strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
                                                by_layer_x = TRUE)) +
  theme(panel.spacing = unit(3, "lines"),
        plot.title = element_text(size = 22,
                                  face = "bold",
                                  hjust = 0.5),
        plot.subtitle = element_text(size = 20,
                                     hjust = 0.5),
        legend.position="bottom",
        legend.text=element_text(size=18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(hjust=0),
        text = element_text(size=22,
                            family = "Times"),
        axis.line.x = element_line(color="black", linewidth = 1),
        axis.line.y = element_line(color="black", linewidth = 1)) +
  guides(col = guide_legend(title.position = "top", title.hjust =0.5, 
                            title.theme = element_text(size = 18)))

ggsave(filename = "plot_allcause.eps",
       path = "XXXResults/Graphs", 
       width = 24, 
       height = 16, 
       device= cairo_ps, 
       bg="white",
       dpi=700)



# COVID-19 mortality
plot_covid19 <- models_results %>% 
  mutate(follow_up = ifelse(follow_up == "28", "up to 28 days", "up to 60 days"),
         time_period = str_replace(time_period, "period", "epoch")) %>% 
  filter(adjusted == "adjusted",
         mortality == "COVID-19 mortality",
         cohort_time != "psychotic disorders: period 1",
         cohort_time != "substance use disorders: period 1") %>% 
  arrange(order_cohort, order_time_period) %>% 
  ggplot(aes(x = time_period,
             y = estimate,
             color = follow_up)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                linewidth = 1,
                width = 0.25,
                position = position_dodge(width = 0.4)) +
  scale_y_continuous(trans = "log10", 
                     limit = c(0.45, 4.4), 
                     breaks = c(0.5, 1, 2)) +
  geom_hline(yintercept = 1, 
             linetype = 2,
             linewidth = 0.5) +
  labs(x = "",
       y = "adjusted hazard ratio (95% confidence interval)",
       color = "Follow-up period",
       title = "Relative risk of COVID-19 mortality following SARS-CoV-2 infection in people with pre-existing mental disorders") +
  theme_minimal() +
  facet_nested_wrap(~ cohort + exposure_definition,
                           scales = "free",
                           nrow = 5,
                           strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
                                                by_layer_x = TRUE)) +
  theme(panel.spacing = unit(3, "lines"),
        plot.title = element_text(size = 22,
                                  face = "bold",
                                  hjust = 0.5),
        plot.subtitle = element_text(size = 20,
                                     hjust = 0.5),
        legend.position="bottom",
        legend.text=element_text(size=18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(hjust=0),
        text = element_text(size=22,
                            family = "Times"),
        axis.line.x = element_line(color="black", linewidth = 1),
        axis.line.y = element_line(color="black", linewidth = 1)) +
  guides(col = guide_legend(title.position = "top", title.hjust =0.5, 
                            title.theme = element_text(size = 18)))

ggsave(filename = "plot_covid19.eps",
       path = "XXXResults/Graphs", 
       width = 24, 
       height = 16, 
       device= cairo_ps, 
       bg="white",
       dpi=700)


# Kaplan-Meier plots
KM_models <- data_main_analysis %>% 
  group_by(cohort, time_period, exposure_definition) %>% 
  nest() %>%   
  mutate(covid_28 = map(.x = data,
                        ~ survfit(Surv(time_to_death_covid_28, death_covid_28) ~ exposure, data = .x)),
         covid_60 = map(.x = data,
                        ~ survfit(Surv(time_to_death_covid_60, death_covid_60) ~ exposure, data = .x)),
         all_cause_28 = map(.x = data,
                            ~ survfit(Surv(time_to_death_allcause_28, death_allcause_28) ~ exposure, data = .x)),
         all_cause_60 = map(.x = data,
                            ~ survfit(Surv(time_to_death_allcause_60, death_allcause_60) ~ exposure, data = .x)))

KM_plots <- KM_models %>%   
  mutate(title_name = list(cohort),
         exp_def = list(case_when(exposure_definition == "diagnosed" ~ "(diagnosed)",
                                  exposure_definition == "diagnosed and treated" ~ "(diagnosed and treated)")),
         time_name = list(str_remove(as.character(time_period), "period")),
         plot_name = paste(cohort, exposure_definition, time_period, ".eps", sep = "_")) %>% 
  mutate(KM_covid_28 = pmap(across(c(covid_28, data, title_name, exp_def, time_name)),
                            ~ ggsurvplot(fit = ..1,
                                         data = ..2,
                                         title = paste0("COVID-19 mortality up to 28 days in people with ", ..3, " ", ..4, " in epoch", ..5),
                                         fun = "event",
                                         legend.title	= "",
                                         legend.labs = c(paste0("Matched counterparts without ", ..3), paste0("Individuals with ", ..3)),
                                         ylab = "Cumulative event (95% CI)",
                                         tables.y.text.col = FALSE,
                                         risk.table = TRUE,
                                         cumevents = TRUE,
                                         cumcensor = TRUE,
                                         censor = FALSE,
                                         tables.height = 0.15,
                                         conf.int = TRUE)),
         KM_covid_60 = pmap(across(c(covid_60, data, title_name, exp_def, time_name)),
                            ~ ggsurvplot(fit = ..1,
                                         data = ..2,
                                         title = paste0("COVID-19 mortality up to 60 days in people with ", ..3, " ", ..4, " in epoch", ..5),
                                         fun = "event",
                                         legend.title	= "",
                                         legend.labs = c(paste0("Matched counterparts without ", ..3), paste0("Individuals with ", ..3)),
                                         ylab = "Cumulative event (95% CI)",
                                         tables.y.text.col = FALSE,
                                         risk.table = TRUE,
                                         cumevents = TRUE,
                                         cumcensor = TRUE,
                                         censor = FALSE,
                                         tables.height = 0.15,
                                         conf.int = TRUE)),
         KM_allcause_28 = pmap(across(c(all_cause_28, data, title_name, exp_def, time_name)),
                               ~ ggsurvplot(fit = ..1,
                                            data = ..2,
                                            title = paste0("All-cause mortality up to 28 days in people with ", ..3, " ", ..4, " in epoch", ..5),
                                            fun = "event",  
                                            legend.title	= "",
                                            legend.labs = c(paste0("Matched counterparts without ", ..3), paste0("Individuals with ", ..3)),
                                            ylab = "Cumulative event (95% CI)",
                                            tables.y.text.col = FALSE,
                                            risk.table = TRUE,
                                            cumevents = TRUE,
                                            cumcensor = TRUE,
                                            censor = FALSE,
                                            tables.height = 0.15,
                                            conf.int = TRUE)),
         KM_allcause_60 = pmap(across(c(all_cause_60, data, title_name, exp_def, time_name)),
                               ~ ggsurvplot(fit = ..1,
                                            data = ..2,
                                            title = paste0("All-cause mortality up to 60 days in people with ", ..3, " ", ..4, " in epoch", ..5),
                                            fun = "event",
                                            legend.title	= "",
                                            legend.labs = c(paste0("Matched counterparts without ", ..3), paste0("Individuals with ", ..3)),
                                            ylab = "Cumulative event (95% CI)",
                                            tables.y.text.col = FALSE,
                                            risk.table = TRUE,
                                            cumevents = TRUE,
                                            cumcensor = TRUE,
                                            censor = FALSE,
                                            tables.height = 0.15,
                                            conf.int = TRUE)))


KMplots_covid_28 <- lapply(KM_plots$KM_covid_28, ggsave_workaround)
KMplots_covid_60 <- lapply(KM_plots$KM_covid_60, ggsave_workaround)
KMplots_allcause_28 <- lapply(KM_plots$KM_allcause_28, ggsave_workaround)
KMplots_allcause_60 <- lapply(KM_plots$KM_allcause_60, ggsave_workaround)

pwalk(list(KM_plots$plot_name, KMplots_covid_28),
      ggsave, 
      path = "XXXResults/Graphs/Supplement/covid_28", 
      width = 24, 
      height = 16, 
      device= cairo_ps,
      dpi=700)

pwalk(list(KM_plots$plot_name, KMplots_covid_60),
      ggsave, 
      path = "XXXResults/Graphs/Supplement/covid_60", 
      width = 24, 
      height = 16, 
      device= cairo_ps,
      dpi=700)

pwalk(list(KM_plots$plot_name, KMplots_allcause_28),
      ggsave, 
      path = "XXXResults/Graphs/Supplement/allcause_28", 
      width = 24, 
      height = 16, 
      device= cairo_ps,
      dpi=700)

pwalk(list(KM_plots$plot_name, KMplots_allcause_60),
      ggsave, 
      path = "XXXResults/Graphs/Supplement/allcause_60", 
      width = 24, 
      height = 16, 
      device= cairo_ps,
      dpi=700)


