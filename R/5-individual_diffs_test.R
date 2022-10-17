rm(list=ls())
source("dmc/dmc.R")
source("dmc/dmc_extras.R")
source("R/0-analysis_functions.R")
load_model("LBA", "lba_B.R")


fitpath <- file.path(set_fit_path(), "Reliability_Analysis")
loadpath <- create_loadpath(fitpath)
savepath <- create_savepath(fitpath)

loadpath("CA_top_samples_A_lb_final.RData")

CA_top_samples <- CA_top_samples_A_lb

tmp <- 1:24
tmp <- tmp[!tmp==4]
full_balance <- c(tmp, 28) 

CA_top_samples <- CA_top_samples[names(CA_top_samples) %in% full_balance]

theme_set(theme_simple())




# Next look at correlations across participants between inhibition/excitation
# and automation accuracy effects (benefits on automation correct, 
# costs on automation incorrect)

#Define a function that obtains this measure

automation_effects <- function (currentsim) {

  benefit=NA;cost=NA
  
  nonfail_accuracy_manual <- 
    mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="nonf"],2,2)==
                                    tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="nonf"]) 
  )
  
  nonfail_accuracy_H <- 
    mean(substr(currentsim$S[currentsim$cond=="H" & currentsim$failtrial=="nonf"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="H" & currentsim$failtrial=="nonf"])
  )
  
  nonfail_accuracy_L <- 
    mean(substr(currentsim$S[currentsim$cond=="L" & currentsim$failtrial=="nonf"],2,2)==
                               tolower(currentsim$R[currentsim$cond=="L" & currentsim$failtrial=="nonf"])
  )
  
  fail_accuracy_manual <- 
    mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="fail"])) 
  

  fail_accuracy_H <- 
    mean(substr(currentsim$S[currentsim$cond=="H" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="H" & currentsim$failtrial=="fail"]))
  
  
  fail_accuracy_L <- 
    mean(substr(currentsim$S[currentsim$cond=="L" & currentsim$failtrial=="fail"],2,2)==
                            tolower(currentsim$R[currentsim$cond=="L" & currentsim$failtrial=="fail"]))
  
  
  out <- c(nonfail_accuracy_H -nonfail_accuracy_manual, fail_accuracy_manual- fail_accuracy_H,
           nonfail_accuracy_L -nonfail_accuracy_manual, fail_accuracy_manual- fail_accuracy_L)
  names(out) <- c("H_benefit", "H_cost", "L_benefit","L_cost")
  out
  
}

#
library(dplyr)
library(readxl)
library(tidyverse)

Auto_Trust_Scores <- read_csv("Auto_Trust_Scores.csv")

Auto_Trust_Long <- Auto_Trust_Scores %>%  
  pivot_longer(cols=starts_with("D"), 
               names_to = c("Session", "Condition", "Q"),
             names_pattern= "D(.)_(.)_Q(.)")

p_trust_scores <- Auto_Trust_Long %>% group_by(Participant, Condition)%>% 
  summarise(score=mean(value, na.rm=T))

names(p_trust_scores)[1]<- "s"
p_trust_scores$s <- as.character(p_trust_scores$s)

#get list of participant data

data<- get.hdata.dmc(CA_top_samples)

#make data frame containing relevant effects

for (i in unique(data$s)) {
  effects<- automation_effects(data[data$s==i,])
  # browser()
  effects["s"] <- i
  if (i ==1) out <- effects else out <- rbind(out,effects)
}

out <- as.data.frame(
  out) %>% 
  left_join(p_trust_scores %>% filter(Condition=="H"), by="s")

# for(i in unique(cleandats$s)) {
#   print(i)
#   effects2 <- automation_effects(cleandats %>% filter(s==i))
#   effects2$s <- i
#   
#   if (i ==1) out2 <- effects2 else out2 <- rbind(out2,effects2)
#   
# }


out <- as.data.frame(
  out) %>% 
  left_join(p_trust_scores %>% filter(Condition=="L"), by="s")

#same functions as up above but cor.plausible function collapses all dimensions
# before calculation

inhibition_corplausible_H <- function(x) (
                     (
                        #non-conflict, non-fail, false reduction   
                        (x["mean_v.nn.M.nonf.false"] - x["mean_v.nn.H.nonf.false"]) +
                        #non-conflict, fail, true reduction
                        (x["mean_v.nn.M.fail.true"] - x["mean_v.nn.H.fail.true"]) +
                        #conflict, non-fail, false reduction  
                        (x["mean_v.cc.M.nonf.false"] - x["mean_v.cc.H.nonf.false"]) +
                        #conflict, fail, true reduction
                        (x["mean_v.cc.M.fail.true"] - x["mean_v.cc.H.fail.true"]) 
                     )/4
)


excitation_corplausible_H <- function(x) (
                     (
                        #non-conflict non-fail true increase   
                        (x["mean_v.nn.H.nonf.true"] - x["mean_v.nn.M.nonf.true"]) +
                        #non-conflict fail false increase
                        (x["mean_v.nn.H.fail.false"] - x["mean_v.nn.M.fail.false"]) +
                        #conflict non-fail true increase
                        (x["mean_v.cc.H.nonf.true"] - x["mean_v.cc.M.nonf.true"]) +
                        #conflict fail false increase
                        (x["mean_v.cc.H.fail.false"] - x["mean_v.cc.M.fail.false"]) 
                     )/4
)

inhibition_corplausible_L <- function(x) (
  (
    #non-conflict, non-fail, false reduction   
    (x["mean_v.nn.M.nonf.false"] - x["mean_v.nn.L.nonf.false"]) +
      #non-conflict, fail, true reduction
      (x["mean_v.nn.M.fail.true"] - x["mean_v.nn.L.fail.true"]) +
      #conflict, non-fail, false reduction  
      (x["mean_v.cc.M.nonf.false"] - x["mean_v.cc.L.nonf.false"]) +
      #conflict, fail, true reduction
      (x["mean_v.cc.M.fail.true"] - x["mean_v.cc.L.fail.true"]) 
  )/4
)


excitation_corplausible_L <- function(x) (
  (
    #non-conflict non-fail true increase   
    (x["mean_v.nn.L.nonf.true"] - x["mean_v.nn.M.nonf.true"]) +
      #non-conflict fail false increase
      (x["mean_v.nn.L.fail.false"] - x["mean_v.nn.M.fail.false"]) +
      #conflict non-fail true increase
      (x["mean_v.cc.L.nonf.true"] - x["mean_v.cc.M.nonf.true"]) +
      #conflict fail false increase
      (x["mean_v.cc.L.fail.false"] - x["mean_v.cc.M.fail.false"]) 
  )/4
)

#See dmc tutorial "plausible" https://osf.io/pbwx8/


#inhibition, costs and benefits

get_corplausible_MCI <- function(samples, cv, p.name, n, fun, kappa=1){
  
  cor.r <- cor.plausible(samples,fun=fun,
                         cv=cv, p.name=p.name)
  
  dens.r <- postRav(r=cor.r,n=n,spacing=.01,kappa=kappa) 
  
  MCIs <- c(postRav.mean(dens.r), postRav.ci(dens.r,interval=c(.025,.975)))
  names(MCIs) <- c("M", "LCI", "HCI")
  MCIs
}

out$H_benefit <- as.numeric(out$H_benefit)
out$H_cost <- as.numeric(out$H_cost)

out$L_benefit <- as.numeric(out$L_benefit)
out$L_cost <- as.numeric(out$L_cost)

inh_H_benefit <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)


ex_H_benefit <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)




inh_L_benefit <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)

ex_L_benefit <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)





inh_H_cost <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)

ex_H_cost <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)



inh_L_cost <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)

ex_L_cost <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)






inh_H_trust <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="score.x", n=24)

ex_H_trust <-get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="score.x", n=24)


inh_L_trust <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="score.y", n=24)

ex_L_trust <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="score.y", n=24)


tmp <- rbind(inh_H_benefit, ex_H_benefit, inh_L_benefit, ex_L_benefit,
      inh_H_cost, ex_H_cost, inh_L_cost, ex_L_cost,
      inh_H_trust, ex_H_trust, inh_L_trust, ex_L_trust
      )


#same thing but with RT effects




automation_effects <- function (currentsim) {
  
  H_nonf_decrease <- NA; L_nonf_decrease <- NA;
  H_fail_increase <- NA; L_fail_increase <- NA;
  
  RT_tibble <- currentsim %>% filter(toupper(substr(S,1,1))==R) %>%
    group_by(failtrial, cond) %>% 
    summarise(RT= mean(RT))
  
  H_nonf_decrease <- RT_tibble$RT[RT_tibble$failtrial=="nonf" & RT_tibble$cond=="M"] -
    RT_tibble$RT[RT_tibble$failtrial=="nonf" & RT_tibble$cond=="H"]
  
  L_nonf_decrease <- RT_tibble$RT[RT_tibble$failtrial=="nonf" & RT_tibble$cond=="M"] -
    RT_tibble$RT[RT_tibble$failtrial=="nonf" & RT_tibble$cond=="L"]
  
  H_fail_increase <- RT_tibble$RT[RT_tibble$failtrial=="fail" & RT_tibble$cond=="H"] -
    RT_tibble$RT[RT_tibble$failtrial=="fail" & RT_tibble$cond=="M"]
  
  L_fail_increase <- RT_tibble$RT[RT_tibble$failtrial=="fail" & RT_tibble$cond=="L"] -
    RT_tibble$RT[RT_tibble$failtrial=="fail" & RT_tibble$cond=="M"]
  
  
  out <- c(H_nonf_decrease, L_nonf_decrease, H_fail_increase, L_fail_increase)
  
  names(out) <- c(
    "H_benefit", "L_benefit", "H_cost", "L_cost"
  )
  
  out
  
}


rm(out)
#
library(dplyr)
library(readxl)
library(tidyverse)

Auto_Trust_Scores <- read_csv("Auto_Trust_Scores.csv")

Auto_Trust_Long <- Auto_Trust_Scores %>%  
  pivot_longer(cols=starts_with("D"), 
               names_to = c("Session", "Condition", "Q"),
               names_pattern= "D(.)_(.)_Q(.)")

p_trust_scores <- Auto_Trust_Long %>% group_by(Participant, Condition)%>% 
  summarise(score=mean(value, na.rm=T))

names(p_trust_scores)[1]<- "s"
p_trust_scores$s <- as.character(p_trust_scores$s)

#get list of participant data

data<- get.hdata.dmc(CA_top_samples)

#make data frame containing relevant effects

for (i in unique(data$s)) {
  effects<- automation_effects(data[data$s==i,])
  effects["s"] <- i
  if (i ==1) out <- effects else out <- rbind(out,effects)
}

out <- as.data.frame(
  out) %>% 
  left_join(p_trust_scores %>% filter(Condition=="H"), by="s")

# for(i in unique(cleandats$s)) {
#   print(i)
#   effects2 <- automation_effects(cleandats %>% filter(s==i))
#   effects2$s <- i
#   
#   if (i ==1) out2 <- effects2 else out2 <- rbind(out2,effects2)
#   
# }


out <- as.data.frame(
  out) %>% 
  left_join(p_trust_scores %>% filter(Condition=="L"), by="s")

#same functions as up above but cor.plausible function collapses all dimensions
# before calculation

inhibition_corplausible_H <- function(x) (
  (
    #non-conflict, non-fail, false reduction   
    (x["mean_v.nn.M.nonf.false"] - x["mean_v.nn.H.nonf.false"]) +
      #non-conflict, fail, true reduction
      (x["mean_v.nn.M.fail.true"] - x["mean_v.nn.H.fail.true"]) +
      #conflict, non-fail, false reduction  
      (x["mean_v.cc.M.nonf.false"] - x["mean_v.cc.H.nonf.false"]) +
      #conflict, fail, true reduction
      (x["mean_v.cc.M.fail.true"] - x["mean_v.cc.H.fail.true"]) 
  )/4
)


excitation_corplausible_H <- function(x) (
  (
    #non-conflict non-fail true increase   
    (x["mean_v.nn.H.nonf.true"] - x["mean_v.nn.M.nonf.true"]) +
      #non-conflict fail false increase
      (x["mean_v.nn.H.fail.false"] - x["mean_v.nn.M.fail.false"]) +
      #conflict non-fail true increase
      (x["mean_v.cc.H.nonf.true"] - x["mean_v.cc.M.nonf.true"]) +
      #conflict fail false increase
      (x["mean_v.cc.H.fail.false"] - x["mean_v.cc.M.fail.false"]) 
  )/4
)

inhibition_corplausible_L <- function(x) (
  (
    #non-conflict, non-fail, false reduction   
    (x["mean_v.nn.M.nonf.false"] - x["mean_v.nn.L.nonf.false"]) +
      #non-conflict, fail, true reduction
      (x["mean_v.nn.M.fail.true"] - x["mean_v.nn.L.fail.true"]) +
      #conflict, non-fail, false reduction  
      (x["mean_v.cc.M.nonf.false"] - x["mean_v.cc.L.nonf.false"]) +
      #conflict, fail, true reduction
      (x["mean_v.cc.M.fail.true"] - x["mean_v.cc.L.fail.true"]) 
  )/4
)


excitation_corplausible_L <- function(x) (
  (
    #non-conflict non-fail true increase   
    (x["mean_v.nn.L.nonf.true"] - x["mean_v.nn.M.nonf.true"]) +
      #non-conflict fail false increase
      (x["mean_v.nn.L.fail.false"] - x["mean_v.nn.M.fail.false"]) +
      #conflict non-fail true increase
      (x["mean_v.cc.L.nonf.true"] - x["mean_v.cc.M.nonf.true"]) +
      #conflict fail false increase
      (x["mean_v.cc.L.fail.false"] - x["mean_v.cc.M.fail.false"]) 
  )/4
)

#See dmc tutorial "plausible" https://osf.io/pbwx8/


#inhibition, costs and benefits

get_corplausible_MCI <- function(samples, cv, p.name, n, fun, kappa=1){
  
  cor.r <- cor.plausible(samples,fun=fun,
                         cv=cv, p.name=p.name)
  
  dens.r <- postRav(r=cor.r,n=n,spacing=.01,kappa=kappa) 
  
  MCIs <- c(postRav.mean(dens.r), postRav.ci(dens.r,interval=c(.025,.975)))
  names(MCIs) <- c("M", "LCI", "HCI")
  MCIs
}

out$H_benefit <- as.numeric(out$H_benefit)
out$H_cost <- as.numeric(out$H_cost)

out$L_benefit <- as.numeric(out$L_benefit)
out$L_cost <- as.numeric(out$L_cost)

inh_H_benefit_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)

ex_H_benefit_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)




inh_L_benefit_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)

ex_L_benefit_RT <-get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)





inh_H_cost_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)

ex_H_cost_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)



inh_L_cost_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)

ex_L_cost_RT <- get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)









tmp <- rbind(inh_H_benefit, ex_H_benefit, inh_L_benefit, ex_L_benefit,
             inh_H_cost, ex_H_cost, inh_L_cost, ex_L_cost,
             inh_H_benefit_RT, ex_H_benefit_RT, inh_L_benefit_RT, ex_L_benefit_RT,
             inh_H_cost_RT, ex_H_cost_RT, inh_L_cost_RT, ex_L_cost_RT,
             inh_H_trust, ex_H_trust, inh_L_trust, ex_L_trust
)

tmp <- as.data.frame(tmp)


tmp$effect <- rownames(tmp)






save(ex_L_benefit, inh_L_benefit, ex_H_benefit, inh_H_benefit,
  ex_L_cost, inh_L_cost, ex_H_cost, inh_H_cost,
  
  ex_L_benefit_RT, inh_L_benefit_RT, ex_H_benefit_RT, inh_H_benefit_RT,
  ex_L_cost_RT, inh_L_cost_RT, ex_H_cost_RT, inh_H_cost_RT,
  
  ex_L_trust, inh_L_trust, ex_H_trust, inh_H_trust,
  file= "model_correlations_A_lb_TRANSFORMED.RData")



MCIs <- function(MCI){
  MCI <- round(MCI, 2)
  paste0(MCI["M"], " (", MCI["LCI"], " - ", MCI["HCI"], ")")
  
}

cortab <- rbind(
  c(MCIs(ex_L_benefit), MCIs(inh_L_benefit), MCIs(ex_H_benefit), MCIs(inh_H_benefit)),
  c(MCIs(ex_L_cost), MCIs(inh_L_cost), MCIs(ex_H_cost), MCIs(inh_H_cost)),
  
  c(MCIs(ex_L_benefit_RT), MCIs(inh_L_benefit_RT), MCIs(ex_H_benefit_RT), MCIs(inh_H_benefit_RT)),
  c(MCIs(ex_L_cost_RT), MCIs(inh_L_cost_RT), MCIs(ex_H_cost_RT), MCIs(inh_H_cost_RT)),
  
  c(MCIs(ex_L_trust), MCIs(inh_L_trust), MCIs(ex_H_trust), MCIs(inh_H_trust))

)

rownames(cortab) <- c("Accuracy Benefit", "Accuracy Cost", "RT Benefit",
                        "RT Cost", "Trust")
colnames(cortab) <- c("Excitation (L)", "Inhibition (L)", "Excitation (H)", "Inhibition (H)")

pandoc.table(cortab)



tmp$cond <- "High Reliability"
tmp$cond[grepl("_L_", tmp$effect)] <- "Low Reliability"

tmp$mech <- "Inhibition"
tmp$mech[grepl("ex", tmp$effect)] <- "Excitation"

tmp$var <- "Accuracy Cost"
tmp$var[grepl("benefit$", tmp$effect)] <- "Accuracy Benefit"

tmp$var[grepl("cost_RT", tmp$effect)] <- "RT Cost"
tmp$var[grepl("benefit_RT", tmp$effect)] <- "RT Benefit"

tmp$var[grepl("trust", tmp$effect)]  <- "Trust"

inhextab <- rbind(
  c("", "", "", ""),
  c(L_conf_ex_success, L_conf_inh_success, H_conf_ex_success, H_conf_inh_success),
  c(L_conf_ex_fail, L_conf_inh_fail, H_conf_ex_fail, H_conf_inh_fail),
  c("", "", "", ""),
  c(L_nonconf_ex_success, L_nonconf_inh_success, H_conf_ex_success, H_nonconf_inh_success),
  c(L_nonconf_ex_fail, L_nonconf_inh_fail, H_conf_ex_fail, H_nonconf_inh_fail)
  
)

rownames(inhextab) <- c("Conflict Trials", "Automation Correct", "Automation Incorrect",
                        "Non-conflict Trials", "Automation Correct", "Automation Incorrect")
colnames(inhextab) <- c("Excitation (L)", "Inhibition (L)", "Excitation (H)", "Inhibition (H)")

pandoc.table(inhextab)




















loadpath("CA_top_samples_pp.RData")






options(dplyr.summarise.inform = FALSE)

loadpath("CA_top_samples_pp.RData")
loadpath("no_thres.RData")
loadpath("no_ex.RData")
loadpath("no_inh.RData")





subject_postexp<- get.subj.effects.m(list(pp, no_inh, no_ex, no_thres), automation_effects , 
                          c("Full Model", "Inhibition Removed", "Excitation Removed", "Fixed Thresholds"))

subject_postexp$effect <- factor(subject_postexp$effect, levels=c("H_benefit", "H_cost",
                                                                  "L_benefit", "L_cost"),
                      labels = c("Automation Benefit (H)", "Automation Cost (H)",
                                 "Automation Benefit (L)", "Automation Cost (L)"))

subject_postexp$model <- factor(subject_postexp$model, levels=c("Full Model","Fixed Thresholds",
                                                                "Inhibition Removed",
                                          "Excitation Removed"))

ggplot(subject_postexp, aes(data, mean)) + geom_point(size=1) + geom_abline(slope=1, intercept=0) +
facet_grid(effect~model, scales = "free") + geom_errorbar(aes(ymax = upper, ymin = lower), alpha=0.3) +
  ylab("Model") + xlab("Observed")+ theme(text = element_text(size = 14))


