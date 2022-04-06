source("dmc/dmc.R")
source("dmc/dmc_extras.R")
load("samples_data/CA_top_samples.RData")
load_model("LBA", "lba_B.R")


path <- "C:/Users/282952E/Dropbox/fits/Reliability_Analysis/"

loadpath <- function(filenam) {
  path <- "C:/Users/282952E/Dropbox/fits/Reliability_Analysis/"
  load(paste0(path, filenam), env=parent.frame())
}


tmp <- 1:24
tmp <- tmp[!tmp==4]
full_balance <- c(tmp, 28) 

CA_top_samples <- CA_top_samples[names(CA_top_samples) %in% full_balance]

theme_set(theme_simple())

#First get plot of individual participant inhibition/excitation levels

#apply a function to an individual participant's thetas
# and get mean + 95% CIs

individual_summary <- function(theta_FUN, samples){
  individual_thetas <- samples$theta
  effect_summary <- theta_FUN(individual_thetas)
  c(mean= mean(effect_summary), LCI= quantile(effect_summary, 0.025),
    HCI= quantile(effect_summary, 0.975))
}

# Deploy individual_summary to get mean + 95% CIs of a function of thetas
# for participant list

get_individual_summaries <- function(theta_FUN, hsamples, effect_name){
  
  effect_df<- as.data.frame(
    t(
      sapply(FUN=function(x)  individual_summary(x, theta_FUN=theta_FUN), 
        hsamples)
      )
  )
  effect_df$effect <- effect_name
  effect_df$participant <- rownames(effect_df)
  colnames(effect_df) <- c("M", "LCI", "HCI", "effect", "participant")
  effect_df
}


#Inhibition and excitation functions- averaged over the effects

H_inhibition <- function(x) (
                     (
                        #non-conflict, non-fail, false reduction   
                        (x[,"mean_v.nn.M.nonf.false",,drop=FALSE] - x[,"mean_v.nn.H.nonf.false",,drop=FALSE]) +
                        #non-conflict, fail, true reduction
                        (x[,"mean_v.nn.M.fail.true",,drop=FALSE] - x[,"mean_v.nn.H.fail.true",,drop=FALSE]) +
                        #conflict, non-fail, false reduction  
                        (x[,"mean_v.cc.M.nonf.false",,drop=FALSE] - x[,"mean_v.cc.H.nonf.false",,drop=FALSE]) +
                        #conflict, fail, true reduction
                        (x[,"mean_v.cc.M.fail.true",,drop=FALSE] - x[,"mean_v.cc.H.fail.true",,drop=FALSE]) 
                     )/4
)


H_excitation <- function(x) (
                     (
                        #non-conflict, non-fail, true increase   
                        (x[,"mean_v.nn.H.nonf.true",,drop=FALSE] - x[,"mean_v.nn.M.nonf.true",,drop=FALSE]) +
                        #non-conflict, fail, false increase
                        (x[,"mean_v.nn.H.fail.false",,drop=FALSE] - x[,"mean_v.nn.M.fail.false",,drop=FALSE]) +
                        #conflict, non-fail, true increase
                        (x[,"mean_v.cc.H.nonf.true",,drop=FALSE] - x[,"mean_v.cc.M.nonf.true",,drop=FALSE]) +
                        #conflict, fail, false increase
                        (x[,"mean_v.cc.H.fail.false",,drop=FALSE] - x[,"mean_v.cc.M.fail.false",,drop=FALSE]) 
                     )/4
)

L_inhibition <- function(x) (
  (
    #non-conflict, non-fail, false reduction   
    (x[,"mean_v.nn.M.nonf.false",,drop=FALSE] - x[,"mean_v.nn.L.nonf.false",,drop=FALSE]) +
      #non-conflict, fail, true reduction
      (x[,"mean_v.nn.M.fail.true",,drop=FALSE] - x[,"mean_v.nn.L.fail.true",,drop=FALSE]) +
      #conflict, non-fail, false reduction  
      (x[,"mean_v.cc.M.nonf.false",,drop=FALSE] - x[,"mean_v.cc.L.nonf.false",,drop=FALSE]) +
      #conflict, fail, true reduction
      (x[,"mean_v.cc.M.fail.true",,drop=FALSE] - x[,"mean_v.cc.L.fail.true",,drop=FALSE]) 
  )/4
)


L_excitation <- function(x) (
  (
    #non-conflict, non-fail, true increase   
    (x[,"mean_v.nn.L.nonf.true",,drop=FALSE] - x[,"mean_v.nn.M.nonf.true",,drop=FALSE]) +
      #non-conflict, fail, false increase
      (x[,"mean_v.nn.L.fail.false",,drop=FALSE] - x[,"mean_v.nn.M.fail.false",,drop=FALSE]) +
      #conflict, non-fail, true increase
      (x[,"mean_v.cc.L.nonf.true",,drop=FALSE] - x[,"mean_v.cc.M.nonf.true",,drop=FALSE]) +
      #conflict, fail, false increase
      (x[,"mean_v.cc.L.fail.false",,drop=FALSE] - x[,"mean_v.cc.M.fail.false",,drop=FALSE]) 
  )/4
)

               




effects_H_inhibition <- get_individual_summaries(
  theta_FUN = H_inhibition, hsamples=CA_top_samples,
                         effect_name= "H_Inhibition")

effects_H_excitation <- get_individual_summaries(
  theta_FUN = H_excitation, hsamples=CA_top_samples,
                         effect_name= "H_Excitation")


effects_L_inhibition <- get_individual_summaries(
  theta_FUN = L_inhibition, hsamples=CA_top_samples,
  effect_name= "L_Inhibition")

effects_L_excitation <- get_individual_summaries(
  theta_FUN = L_excitation, hsamples=CA_top_samples,
  effect_name= "L_Excitation")


all_effects <- rbind(effects_H_excitation,
                            effects_H_inhibition,
                     effects_L_inhibition,
                     effects_L_excitation)

all_effects$participant <- factor(as.numeric(all_effects$participant))

all_effects$Condition <- "L"
all_effects$Condition[grep("H", all_effects$effect)] <- "H"

all_effects$DV <- "Excitation"
all_effects$DV[grep("Inhibition", all_effects$effect)] <- "Inhibition"


ggplot(all_effects, aes(participant, M)) +
    geom_point(stat = "identity", size=2.5) +
    geom_errorbar(aes(ymax = HCI , ymin = LCI, width = 0.3)) +
    geom_hline(aes(yintercept=0), linetype=2)+
    facet_grid(Condition~DV) +xlab("Participant") +ylab ("Effect")


# Next look at correlations across participants between inhibition/excitation
# and automation accuracy effects (benefits on automation correct, 
# costs on automation incorrect)

#Define a function that obtains this measure

automation_effects <- function (currentsim) {

  benefit=NA;cost=NA
  
  nonfail_accuracy_manual <- mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="nonf"],2,2)==
                                    tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="nonf"]))  
  
  nonfail_accuracy_H <- mean(substr(currentsim$S[currentsim$cond=="H" & currentsim$failtrial=="nonf"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="H" & currentsim$failtrial=="nonf"]))  
  
  nonfail_accuracy_L <- mean(substr(currentsim$S[currentsim$cond=="L" & currentsim$failtrial=="nonf"],2,2)==
                               tolower(currentsim$R[currentsim$cond=="L" & currentsim$failtrial=="nonf"])) 
  
  fail_accuracy_manual <- mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="fail"]))  

  fail_accuracy_H <- mean(substr(currentsim$S[currentsim$cond=="H" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="H" & currentsim$failtrial=="fail"])) 
  
  fail_accuracy_L <- mean(substr(currentsim$S[currentsim$cond=="L" & currentsim$failtrial=="fail"],2,2)==
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
  effects["s"] <- i
  if (i ==1) out <- effects else out <- rbind(out,effects)
}

out <- as.data.frame(
  out) %>% 
  left_join(p_trust_scores %>% filter(Condition=="H"), by="s")

for(i in unique(cleandats$s)) {
  print(i)
  effects2 <- automation_effects(cleandats %>% filter(s==i))
  effects2$s <- i
  
  if (i ==1) out2 <- effects2 else out2 <- rbind(out2,effects2)
  
}


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

get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_benefit", n=24)




get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_benefit", n=24)





get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="H_cost", n=24)



get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="L_cost", n=24)






get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="score.x", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_H,
                     cv=as.data.frame(out), 
                     p.name="score.x", n=24)


get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="score.y", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible_L,
                     cv=as.data.frame(out), 
                     p.name="score.y", n=24)






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


