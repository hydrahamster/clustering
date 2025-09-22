## ---------------------------
##
## Script name: age-onset-expl.R
##
## Purpose of script: As an upshot of the ML classification work, differential complications trajectories by age of onset became apparent so this is to explore this story further
##
## Author: Dr. Liza Kretzschmar
##
## Date Created: 2025-09-17
##
##
## ---------------------------
##
## Notes:
##   
##   
##
## ---------------------------

## back end
pacman::p_load(
  tidyverse,
  janitor,
  DataExplorer,
  bodlR,
  gtsummary,
  ggstatsplot
)

data_base <- readRDS("../metadata-processing/data_base_mar25.rds")

data_phenotypes <- readRDS("../metadata-processing/data_phenotypes.rds")

## ---------------------------

# initial clue of things going on
data_expl <- data_phenotypes %>% 
  dplyr::select(prophecyid,
                t2dm_macrovascular_microvascular) %>%
  left_join((data_base %>% dplyr::select(prophecyid,
                                         t2dm_diag_age_onset,
                                         t2dm_diag_yrs_since_adj,
                                         sex)), by = "prophecyid") %>%
  dplyr::filter(!is.na(t2dm_macrovascular_microvascular)) %>%
  dplyr::filter(t2dm_macrovascular_microvascular!= "-T2DM -MACRO -MICRO")

data_expl %>%
  ggplot(aes(x = t2dm_diag_age_onset, #, #t2dm_diag_yrs_since_adj, #, #clinicage_rnd, #, 
             y = t2dm_diag_yrs_since_adj)) +
  geom_point(aes(shape = t2dm_macrovascular_microvascular, color = t2dm_macrovascular_microvascular), size = 2.5)+
  #, position=position_jitter(h=0.05,w=0.15)) +
  scale_shape_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c(1, 19, 0, 15)) +
  theme_bw() +
  scale_colour_bodl(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Participant age at T2DM diagnosis",
    y = "Time since T2DM diagnosis",
    colour = "Participant's complication(s)",
    shape = "Participant's complication(s)",
    title = "Participant's complications by age of T2DM onset and duration of T2DM"
  )

## ---------------------------

# Age of onset strtification and characterisation

data_onset <- data_phenotypes %>% 
  dplyr::select(-c(dr_early,
                   ckd_early,
                   drpresence, 
                   t2dm)) %>%
  left_join(data_base, by = "prophecyid") %>%
  mutate(onset_cat = case_when(
    t2dm_diag_age_onset < 36 ~ "early",
    t2dm_diag_age_onset > 59 ~ "late",
    is.na(t2dm_diag_age_onset) ~ NA_character_,
    TRUE ~ "mid"
  )) %>%
  mutate(onset_cat = as_factor(onset_cat)) %>%
  mutate(onset_cat = fct_relevel(onset_cat, "late", after = Inf))  %>%
  mutate(onset_cat = fct_relevel(onset_cat, "early")) %>%
  dplyr::filter(!is.na(onset_cat)) %>%
  mutate(duration_cat = case_when(
    t2dm_diag_yrs_since_adj < 11 ~ "0-10 yrs",
    t2dm_diag_yrs_since_adj < 21 ~ "11-20 yrs",
    t2dm_diag_yrs_since_adj < 31 ~ "21-30 yrs",
    t2dm_diag_yrs_since_adj < 41 ~ "31-40 yrs",
    t2dm_diag_yrs_since_adj < 51 ~ "41-50 yrs",
    is.na(t2dm_diag_yrs_since_adj) ~ NA_character_,
    TRUE ~ "WTF"
  )) %>%
  mutate(duration_cat = as_factor(duration_cat)) %>%
  mutate(duration_cat = fct_relevel(duration_cat, "11-20 yrs", after = 1))  %>%
  mutate(duration_cat2 = case_when(
    t2dm_diag_yrs_since_adj < 5 ~ "0-4 yrs",
    t2dm_diag_yrs_since_adj < 10 ~ "5-9 yrs",
    t2dm_diag_yrs_since_adj < 15 ~ "10-14 yrs",
    t2dm_diag_yrs_since_adj < 20 ~ "15-19 yrs",
    t2dm_diag_yrs_since_adj < 25 ~ "20-24 yrs",
    t2dm_diag_yrs_since_adj < 30 ~ "25-29 yrs",
    t2dm_diag_yrs_since_adj > 29 ~ "30+ yrs",
    is.na(t2dm_diag_yrs_since_adj) ~ NA_character_,
    TRUE ~ "WTF"
  )) %>%
  mutate(duration_cat2 = as_factor(duration_cat2)) %>%
  mutate(duration_cat2 = fct_relevel(duration_cat2, c("0-4 yrs", "5-9 yrs", "10-14 yrs", "15-19 yrs", "20-24 yrs", "25-29 yrs", "30+ yrs")))  %>%
  mutate(across(where(is.factor), fct_drop))

# table for comparing 
data_onset %>%
  dplyr::select(onset_cat,
                t2dm_macrovascular_microvascular,
                sex,
                clinicage_rnd,
                remoteness2,
                diabhxfam,
                maritalstatus2,
                educhighest_adj,
                unemployed,,
                incomedayslast2grp,
                #residents,
                incomeperscat,
                educfurther1,
                smoking,
                alcoauditc1,
                #phq9score,
                #chronicstress_score,
                #racism_score,
                bmi, #anthropomorphic
                waist_mean,
                bpdia, #BP
                bpsys,
                hba1c, #glycemic
                glucose,
                chol_total,#lipids
                hdl_comb,
                ldl_comb,
                trig_comb,
                lpa,
                apoa1,
                apob,
                potassium_final, #blood chem
                sodium_final,
                chloride_final,
                blood_urea_final,
                alt, #liver
                ast,
                ggt,
                crphs, #inflammation
                #haematocrit_final,#blood content
                haemoglob_final,
                egfr_combined, #renal
                acr) %>%
  dplyr::filter(diabhxfam != "Refused/Not stated")%>%
  mutate(across(where(is.factor), fct_drop)) %>%
  tbl_summary(by = onset_cat,
              missing = "no",
              label = list(clinicage_rnd ~ "Participant age",
                           remoteness2 ~ "Remoteness",
                           diabhxfam ~ "Family history T2DM",
                           maritalstatus2 ~ "Martial status",
                           educhighest_adj ~ "Education level",
                           incomedayslast2grp ~ "How long does income last",
                           incomeperscat ~ "Fortnightly income",
                           t2dm_macrovascular_microvascular ~ "Cardiometabolic complications")
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>%
  add_n() %>%
  bold_labels()

# Overall comparison of complications by onset category
data_onset %>% 
  ggplot(aes(onset_cat)) + 
  geom_bar(aes(fill = t2dm_macrovascular_microvascular)) +
  theme_bw() +
  scale_fill_bodl(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"))

# Complications status by age of onset and time since
data_onset %>% 
  ggplot(aes(duration_cat2)) + 
  geom_bar(aes(fill = t2dm_macrovascular_microvascular)) +
  theme_bw() +
  scale_fill_bodl(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"))+
  facet_grid( ~ onset_cat)

# extratc which participants do not develop complications after 20 years
check <- data_onset %>%
  dplyr::filter(t2dm_macrovascular_microvascular == "+T2DM -MACRO -MICRO") %>% 
  dplyr::filter( duration_cat2 != "0-10 yrs") %>% 
  dplyr::select(t2dm_diag_yrs_since_adj, clinicage_rnd, duration_cat, onset_cat)

no_trouble <- data_onset %>%
  dplyr::filter(t2dm_macrovascular_microvascular == "+T2DM -MACRO -MICRO") %>%
  dplyr::filter(duration_cat2 != "0-4 yrs" & duration_cat2 != "5-9 yrs" & duration_cat2 != "10-14 yrs")  

no_trouble %>%
  dplyr::filter(onset_cat == "mid") %>%
  dplyr::select(prophecyid)

# baseline participants
data_zerodx <- data_onset %>%
  dplyr::filter(t2dm_diag_yrs_since_adj < 1)

data_zerodx %>%
  ggplot(aes(t2dm_macrovascular_microvascular)) + 
  geom_bar(aes(fill = t2dm_macrovascular_microvascular)) +
  theme_bw() +
  scale_fill_bodl(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO")) +
  labs(x = "Complications") +
  scale_x_discrete(label = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"))

data_zerodx %>%
  ggplot(aes(onset_cat)) + 
  geom_bar(aes(fill = t2dm_macrovascular_microvascular)) +
  theme_bw() +
  scale_fill_bodl(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"))
