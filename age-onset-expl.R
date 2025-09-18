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
  bodlR
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
  scale_colour_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c("#009E73", "#E69F00", "#CC79A7", "#D55E00")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Participant age at T2DM diagnosis",
    y = "Time since T2DM diagnosis",
    colour = "Participant's complication(s)",
    shape = "Participant's complication(s)",
    title = "Participant's complications by age of T2DM onset and duration of T2DM"
  )

## ---------------------------
# AUSDRISK comparison

data_ausdr <- data_base %>%
  dplyr::select(prophecyid,
                clinicage_rnd,
                sex,
                diabhxfam,
                diabhxmoth, # not 0
                diabhxfath,
                diabhxsis,
                diabhxbro,
                t2dm,
                prevpregdiabhx_adj,
                currentpregdiab,
                pbs_bp_med,
                smoking,
                nutfruit_comb, # 4 or 5
                nutveg,# 4 or 5
                physmodtimemins_comb_avg_cut, # 0
                waist_mean ) %>%
  mutate(q1_age = case_when(
    clinicage_rnd < 35 ~ 0,
    clinicage_rnd < 45 ~ 2,
    clinicage_rnd < 65 ~ 4,
    clinicage_rnd > 64 ~ 8,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q2_gender = case_when(
    sex == "Female" ~ 0,
    sex == "Male" ~ 3,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q3_ethn = 2)%>%
  mutate(q4_fam = case_when(
    diabhxfam == "No" ~ 0, # no fam history of T2DM
    ((diabhxmoth != 0 & !is.na(diabhxmoth)) | (diabhxfath != 0 & !is.na(diabhxfath)) | (diabhxsis != 0 & !is.na(diabhxsis))  | (diabhxbro != 0 & !is.na(diabhxbro))  ) ~ 3, # yes if mother/father/sister/brother count of T2DM isn't zero or missing
    diabhxfam == "Yes" ~ 0, # remaining entries for family history would capture grandparents which isn't in the risk score
    TRUE ~ NA_integer_
  ))%>%
  mutate(q5_gluc = case_when(
    prevpregdiabhx_adj == "Yes" ~ 6, 
    t2dm == "Yes" ~ 6,
    currentpregdiab == "Yes" ~ 6,
    TRUE ~ 0# All participants have T2DM status, so no NA capture needed
  )) %>%
  mutate(q6_bpmed = case_when(
    pbs_bp_med == "Yes - BP treatment" ~ 2,
    pbs_bp_med == "No - BP treatment" ~ 0,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q7_smoke = case_when(
    smoking == "Yes" ~ 2,
    is.na(smoking) ~ NA_integer_,
    TRUE ~ 0
  ))%>%
  mutate(q8_food = case_when(
    (nutfruit_comb == "Once a day" | nutfruit_comb == "More than once a day" | nutveg == "Once a day" | nutveg == "More than once a day") ~ 0,
    (is.na(nutfruit_comb) & is.na(nutveg)) ~ NA_integer_,
    TRUE ~ 1
  ))%>%
  mutate(q9_sed = case_when(
    physmodtimemins_comb_avg_cut == "Yes - Physical Activity >=150 minutes" ~ 0,
    is.na(physmodtimemins_comb_avg_cut) ~ NA_integer_,
    TRUE ~ 2
  ))%>%
  mutate(q10_waist = case_when(
    (sex == "Female" & waist_mean < 80) ~ 0,
    (sex == "Female" & waist_mean < 91) ~ 4,
    (sex == "Female" & waist_mean > 90) ~ 7,
    (sex == "Male" & waist_mean < 90) ~ 0,
    (sex == "Male" & waist_mean < 101) ~ 4,
    (sex == "Male" & waist_mean > 100) ~ 7,
    TRUE ~ NA_integer_
  )) %>%
  mutate(ausdr_na = rowSums(is.na(across(q1_age:q10_waist)))) %>%
  mutate(ausdr_score = rowSums(across(q1_age:q10_waist), na.rm = T))

data_ausdr %>%
  dplyr::filter(ausdr_na == 0) %>%
  ggplot(aes(x = ausdr_score, #, #t2dm_diag_yrs_since_adj, #, #clinicage_rnd, #, 
             y = clinicage_rnd)) +
  geom_rect(xmin = -Inf, xmax = 6, ymin = -Inf, ymax = Inf, fill = "#75B24D", alpha = 0.002) +
  geom_rect(xmin = 6, xmax = 11, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.002) +
  geom_rect(xmin = 11, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#9E2900", alpha = 0.002) +
  geom_point(aes(color = t2dm), size = 2.5)+
  #, position=position_jitter(h=0.05,w=0.15)) +
  #scale_shape_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c(1, 19, 0, 15)) +
  theme_bw() +
  scale_colour_bodl(palette = "bodl_full", values = c("DesertFlameYellow", "BushTomatoRed")) +
  #scale_colour_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c("#009E73", "#E69F00", "#CC79A7", "#D55E00")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Participant AUSDRisk score",
    y = "Participant age",
    colour = "T2DM status",
    title = "Participant's AUSDRisk score"
  )

data_ausdr9 <- data_base %>%
  dplyr::select(prophecyid,
                clinicage_rnd,
                sex,
                diabhxfam,
                diabhxmoth, # not 0
                diabhxfath,
                diabhxsis,
                diabhxbro,
                t2dm,
                prevpregdiabhx_adj,
                currentpregdiab,
                pbs_bp_med,
                smoking,
                nutfruit_comb, # 4 or 5
                nutveg,# 4 or 5
                physmodtimemins_comb_avg_cut, # 0
                waist_mean ) %>%
  mutate(q1_age = case_when(
    clinicage_rnd < 35 ~ 0,
    clinicage_rnd < 45 ~ 2,
    clinicage_rnd < 65 ~ 4,
    clinicage_rnd > 64 ~ 8,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q2_gender = case_when(
    sex == "Female" ~ 0,
    sex == "Male" ~ 3,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q3_ethn = 2)%>%
  mutate(q4_fam = case_when(
    diabhxfam == "No" ~ 0, # no fam history of T2DM
    ((diabhxmoth != 0 & !is.na(diabhxmoth)) | (diabhxfath != 0 & !is.na(diabhxfath)) | (diabhxsis != 0 & !is.na(diabhxsis))  | (diabhxbro != 0 & !is.na(diabhxbro))  ) ~ 3, # yes if mother/father/sister/brother count of T2DM isn't zero or missing
    diabhxfam == "Yes" ~ 0, # remaining entries for family history would capture grandparents which isn't in the risk score
    TRUE ~ NA_integer_
  ))%>%
  mutate(q5_gluc = case_when(
    prevpregdiabhx_adj == "Yes" ~ 6, 
    t2dm == "Yes" ~ 6,
    currentpregdiab == "Yes" ~ 6,
    TRUE ~ 0# All participants have T2DM status, so no NA capture needed
  )) %>%
  mutate(q6_bpmed = case_when(
    pbs_bp_med == "Yes - BP treatment" ~ 2,
    pbs_bp_med == "No - BP treatment" ~ 0,
    TRUE ~ NA_integer_
  ))%>%
  mutate(q7_smoke = case_when(
    smoking == "Yes" ~ 2,
    is.na(smoking) ~ NA_integer_,
    TRUE ~ 0
  ))%>%
  mutate(q8_food = case_when(
    (nutfruit_comb == "Once a day" | nutfruit_comb == "More than once a day" | nutveg == "Once a day" | nutveg == "More than once a day") ~ 0,
    (is.na(nutfruit_comb) & is.na(nutveg)) ~ NA_integer_,
    TRUE ~ 1
  ))%>%
  # mutate(q9_sed = case_when(
  #   physmodtimemins_comb_avg_cut == "Yes - Physical Activity >=150 minutes" ~ 0,
  #   is.na(physmodtimemins_comb_avg_cut) ~ NA_integer_,
  #   TRUE ~ 2
  # ))%>%
  mutate(q10_waist = case_when(
    (sex == "Female" & waist_mean < 80) ~ 0,
    (sex == "Female" & waist_mean < 91) ~ 4,
    (sex == "Female" & waist_mean > 90) ~ 7,
    (sex == "Male" & waist_mean < 90) ~ 0,
    (sex == "Male" & waist_mean < 101) ~ 4,
    (sex == "Male" & waist_mean > 100) ~ 7,
    TRUE ~ NA_integer_
  )) %>%
  mutate(ausdr_na = rowSums(is.na(across(q1_age:q10_waist)))) %>%
  mutate(ausdr_score = rowSums(across(q1_age:q10_waist), na.rm = T))

data_ausdr9 %>%
  dplyr::filter(ausdr_na == 0) %>%
  ggplot(aes(x = ausdr_score, #, #t2dm_diag_yrs_since_adj, #, #clinicage_rnd, #, 
             y = clinicage_rnd)) +
  geom_rect(xmin = -Inf, xmax = 4, ymin = -Inf, ymax = Inf, fill = "#75B24D", alpha = 0.002) +
  geom_rect(xmin = 4, xmax = 9, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.002) +
  geom_rect(xmin = 9, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#9E2900", alpha = 0.002) +
  geom_point(aes(color = t2dm), size = 2.5)+
  #, position=position_jitter(h=0.05,w=0.15)) +
  #scale_shape_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c(1, 19, 0, 15)) +
  theme_bw() +
  scale_colour_bodl(palette = "bodl_full", values = c("DesertFlameYellow", "BushTomatoRed")) +
  #scale_colour_manual(labels = c("-MACRO -MICRO", "-MACRO +MICRO", "+MACRO -MICRO", "+MACRO +MICRO"), values = c("#009E73", "#E69F00", "#CC79A7", "#D55E00")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Participant AUSDRisk score",
    y = "Participant age",
    colour = "T2DM status",
    title = "Participant's AUSDRisk score"
  )

## ---------------------------
