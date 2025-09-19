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
  gtsummary
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


