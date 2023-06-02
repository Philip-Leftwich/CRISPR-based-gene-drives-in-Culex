rm(list=ls())
# Packages====

library(tidyverse)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(glmmTMB)
#______________====

# Data ====
## Assessing gene drive at the kmo locus using an inheritance bias approach
load("data/kmo_cross.RData")

## Data formatting====

kmo_cross <- janitor::clean_names(kmo_cross)

kmo_cross <- kmo_cross %>% 
  mutate(grandparent=case_when(grandparent=="F"~"Female",
                                        grandparent=="M"~"Male")) %>% 
  mutate(parent=case_when(parent=="F"~"Female",
                                     parent=="M"~"Male")) %>% 
  mutate(win = vasa_cas9_2069_g_rna + x2069_g_rna_only,
         loss = vasa_cas9_only + wt)

# Fixed effects model====

lm.model2 <- glm((cbind(win,loss))~ parent*grandparent, family=binomial, data=kmo_cross, na.action=na.fail)

dredged <- dredge(lm.model2)
best.model <- get.models(dredged, subset=1)[[1]]

# Mixed model====

### total homing ignoring sex
lmm.modela <- glmmTMB((cbind(win,loss))~ 1+(1|id), family=binomial, data=kmo_cross)

### best fit model
lmm.model <- glmmTMB((cbind(win,loss))~ parent+(1|id), family=binomial, data=kmo_cross)

### mean estimates
emmeans::emmeans(lmm.model, specs=pairwise  ~parent, type="response")

emmeans::emmeans(lmm.model, specs=~1, type="response")

DHARMa::simulateResiduals(lmm.model, plot = T)





