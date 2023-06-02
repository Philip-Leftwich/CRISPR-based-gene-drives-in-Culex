rm(list=ls())
# Packages====

library(tidyverse)
library(janitor)
library(glmmTMB)
library(DHARMa)
library(emmeans)
#______________====

# Data ====

load("data/white_cross.Rdata")

#______________====

# Cross 1 =====
## Assessing gene drive at the white locus using an inheritance-bias approach

## Data formatting====

gfp_female <- gfp_female %>% mutate(Parent = "female")
gfp_male <- gfp_male %>% mutate(Parent = "male")

gfp_total <- rbind(gfp_male,gfp_female)
gfp_total <- janitor::clean_names(gfp_total)

gfp_total <- 
  gfp_total %>% 
  drop_na(cross_number) %>% 
  mutate(id = row_number()) %>% 
  mutate(gfp_w_indel = as.numeric(str_replace(gfp_w_indel, "-", "0"))) %>% 
  mutate(gfp_w_wt = as.numeric(str_replace(gfp_w_wt, "-", "0"))) %>% 
  mutate(win = gfp_donor_copy,
       loss = gfp_w_indel + gfp_w_wt)

## Mixed-effects model ====

### Single models for each sex====

single_sex_model <- function(sex){
  model <- gfp_total %>% 
    filter(parent == {{sex}}) %>% 
    glmmTMB((cbind(win,loss))~ 1 + (1|parent/id), family=binomial, data=.)
  summary(model)
}

single_sex_model("male")
single_sex_model("female")


## Full model with Cas9 marker ====


rfp_female <- rfp_female %>% mutate(Parent = "female")
rfp_male <- rfp_male %>% mutate(Parent = "male")

rfp_total <- rbind(rfp_female, rfp_male)

rfp_total <- rfp_total %>% 
  mutate(DsRed. = as.numeric(str_replace(DsRed., "-", "0"))) %>% 
  drop_na() %>% 
  rename(Win = DsRed.,
         Loss = DsRed..1,
         id = Cross.number) %>% 
  mutate(Fluoro = "Red")

rfp_total <- janitor::clean_names(rfp_total)

gfp_total <- gfp_total %>% 
  mutate(fluoro = "Green") %>% 
  select(win, loss, id, parent, fluoro)

cross_1 <- rbind(rfp_total, gfp_total)

## Mixed effects model ====

lmer1 <- glmmTMB((cbind(win,loss))~ parent*fluoro + (1|parent/id), family=binomial, data=cross_1)

simulateResiduals(fittedModel = lmer1, plot = T)

means_summary <- emmeans::emmeans(lmer1, specs = ~ parent*fluoro, type="response")


#______________====

# Cross 2====
## Assessing gene drive at the white locus using a marked chromosome approach

## Data formatting====

cross_2 <- lapply(list(cross_2_gfp, cross_2_rfp), janitor::clean_names)

cross_2 <- list(cross_2_gfp, cross_2_rfp) %>% 
  map(~ janitor::clean_names(.) %>% 
        drop_na(f2_raft))

cross_2[[1]] <- 
  cross_2[[1]] %>% 
  mutate(gfp_w_copy = as.numeric(str_replace(gfp_w_copy, "-", "0"))) %>% 
  mutate(Win = gfp_w_donor+gfp_w_copy,
         Loss = non_gfp_w_nhej_wt) %>% 
  mutate(fluoro = "Green") 

cross_2[[2]] <- 
  cross_2[[2]] %>% 
    mutate(Win = ds_red,
         Loss = ds_red_1) %>% 
  mutate(fluoro = "Red")

cross_2_total <- cross_2 %>% 
  map_df(~ select(., f2_raft, f1_cage, Win, Loss, fluoro))




## Mixed-effects model====

lmer2<- glmmTMB((cbind(Win,Loss))~ fluoro + (1|f2_raft), family=binomial, data=cross_2_total)

simulateResiduals(fittedModel = lmer2, plot = T)

means_summary2 <- emmeans::emmeans(lmer2, specs= "fluoro", type="response")

## w- conversion analysis ====

lmer3<- glmmTMB(cbind(gfp_w_copy, non_gfp_w_nhej_wt)~ 1, family=binomial, data=cross_2[[1]])

simulateResiduals(fittedModel = lmer3, plot = T)

means_summary3 <- emmeans::emmeans(lmer3, specs= ~ 1, type="response")


## W- vs biased w+====

cross_2_long <- cross_2[[1]] %>% pivot_longer(cols = gfp_w_donor:gfp_w_copy, names_to="chromosome")

lmer4<- glmmTMB(cbind(value, (total-value))~ chromosome + (1|f2_raft), family=binomial, data=cross_2_long)

emmeans::emmeans(lmer4, specs= "chromosome", type="response")
