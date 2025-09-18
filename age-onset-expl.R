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
