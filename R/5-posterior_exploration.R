rm(list=ls())
source("dmc/dmc.R")
source("R/0-analysis_functions.R")
source("dmc/dmc_extras.R")
load_model ("LBA","lba_B.R")
library(stringr)


fitpath <- file.path(set_fit_path(), "Reliability_Analysis")
loadpath <- create_loadpath(fitpath)
savepath <- create_savepath(fitpath)

# loadpath("CA_top_samples_A_lb_final.RData")
# 
# CA_top_samples <- CA_top_samples_A_lb

theme_set(theme_simple())

loadpath("CA_top_samples_pp_A_lb.RData")

options(dplyr.summarise.inform = FALSE)

pnames <- colnames(CA_top_samples[[1]]$theta)

# thres_pickps_set <- rep(
#   pnames[
#    grep("B\\.M", pnames)
#   ], 2
# )
# 
# thres_pickps_other <- c(
#   pnames[grep("B\\.L", pnames)],
#   pnames[grep("B\\.H", pnames)]
#   )
# 
# 
# no_thres <- pickps.h.post.predict.dmc(CA_top_samples,
#                                       save.simulation = TRUE,
#                                       pickps_set=thres_pickps_set,
#                                       pickps_other=thres_pickps_other,
#                                       n.post=200)
# #  
# savepath(no_thres, file= "no_thres_A_lb.RData")

loadpath("no_thres_A_lb.RData")

# 
# ex <- c("mean_v.nn.M.nonf.true", "mean_v.cc.M.nonf.true",
#                    "mean_v.nn.M.fail.false", "mean_v.cc.M.fail.false")
# 
# ex_pickps_set <- rep(ex, 2)
# 
# ex_pickps_other <- c(str_replace(ex, "M", "L"), str_replace(ex, "M", "H"))
# 
# no_ex <- pickps.h.post.predict.dmc(CA_top_samples,
#                                       save.simulation = TRUE,
#                                       pickps_set=ex_pickps_set,
#                                       pickps_other=ex_pickps_other,
#                                    n.post=200)
# 
# savepath(no_ex, file= "no_ex_A_lb.RData")

loadpath("no_ex_A_lb.RData")


# # 
# # 
# # 
# inh <- c("mean_v.nn.M.nonf.false", "mean_v.cc.M.nonf.false",
#         "mean_v.nn.M.fail.true", "mean_v.cc.M.fail.true")
# # 
# inh_pickps_set <- rep(inh, 2)
# 
# inh_pickps_other <- c(str_replace(inh, "M", "L"), str_replace(inh, "M", "H"))
# 
# no_inh <- pickps.h.post.predict.dmc(CA_top_samples,
#                                    save.simulation = TRUE,
#                                    pickps_set=inh_pickps_set,
#                                    pickps_other=inh_pickps_other,
#                                    n.post=200)
# # #
# savepath(no_inh, file= "no_inh_A_lb.RData")

loadpath("no_inh_A_lb.RData")


tmp <- 1:24
tmp <- tmp[!tmp==4]
full_balance <- c(tmp, 28) 

pp <- pp[names(pp) %in% full_balance]
no_inh <- no_inh[names(no_inh) %in% full_balance]
no_ex <- no_ex[names(no_ex) %in% full_balance]
no_thres <- no_thres[names(no_thres) %in% full_balance]


acc_effects <- function (currentsim) {
  
  H_nonf_gain <- NA; L_nonf_gain <- NA;
  H_fail_loss <- NA; L_fail_loss <- NA;
  
  acc_tibble <- currentsim %>% group_by(failtrial, cond) %>% 
    summarise(acc= mean(toupper(substr(S,1,1))==R))
  
  H_nonf_gain <- acc_tibble$acc[acc_tibble$failtrial=="nonf" & acc_tibble$cond=="H"] -
                    acc_tibble$acc[acc_tibble$failtrial=="nonf" & acc_tibble$cond=="M"] 
  
  L_nonf_gain <- acc_tibble$acc[acc_tibble$failtrial=="nonf" & acc_tibble$cond=="L"] -
                     acc_tibble$acc[acc_tibble$failtrial=="nonf" & acc_tibble$cond=="M"] 
  
  H_fail_loss <- acc_tibble$acc[acc_tibble$failtrial=="fail" & acc_tibble$cond=="M"] -
                    acc_tibble$acc[acc_tibble$failtrial=="fail" & acc_tibble$cond=="H"]
  
  L_fail_loss <- acc_tibble$acc[acc_tibble$failtrial=="fail" & acc_tibble$cond=="M"] -
                    acc_tibble$acc[acc_tibble$failtrial=="fail" & acc_tibble$cond=="L"]
     

  
  out <- c(H_nonf_gain, L_nonf_gain, H_fail_loss, L_fail_loss)
  
  names(out) <- c(
    "H_nonf_gain", "L_nonf_gain", "H_fail_loss", "L_fail_loss"
    )
  
  out
  
}

full <- get.effects.dmc(pp, fun = acc_effects)
noinh <- get.effects.dmc(no_inh, fun = acc_effects)
noex <- get.effects.dmc(no_ex, fun = acc_effects)
nothres <- get.effects.dmc(no_thres, fun = acc_effects)

full$model <- "Full"; noinh$model <- "No Inhibition"
noex$model <- "No Excitation"; nothres$model <-"Fixed Thresholds"
full$eff <- rownames(full); noinh$eff <- rownames(noinh); noex$eff <- rownames(noex); nothres$eff <- rownames(nothres)

combined <- rbind(full,noinh, noex, nothres)

combined$model <- factor(combined$model,
                         levels = c("Full", "Fixed Thresholds",
                                    "No Excitation", "No Inhibition"))

combined$Auto <- "Automation Correct (Accuracy Benefit)"
combined$Auto[grep("fail", combined$eff)] <- "Automation Incorrect (Accuracy Cost)"

combined$eff <- factor(combined$eff, levels=c("H_fail_loss",
                                              "H_nonf_gain",
                                              "L_fail_loss",
                                              "L_nonf_gain"),
                       labels=c("High Reliability", "High Reliability",
                                "Low Reliability", "Low Reliability"))

combined %>% ggplot(aes(y=mean, x=eff)) + geom_point(size=2)+ 
  geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes_string(x = 'eff', y= 'data'), pch=21, size=3, colour="black")+
 geom_line(aes(group=1,y=data), linetype=2) +facet_grid(model~Auto) +xlab("")+
  ylab("")




RT_effects <- function (currentsim) {
  
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
    "H_nonf_decrease", "L_nonf_decrease", "H_fail_increase", "L_fail_increase"
  )
  
  out
  
}


full_RT <- get.effects.dmc(pp, fun = RT_effects)
noinh_RT <- get.effects.dmc(no_inh, fun = RT_effects)
noex_RT <- get.effects.dmc(no_ex, fun = RT_effects)
nothres_RT <- get.effects.dmc(no_thres, fun = RT_effects)

full_RT$model <- "Full"; noinh_RT$model <- "No Inhibition"
noex_RT$model <- "No Excitation"; nothres_RT$model <-"Fixed Thresholds"
full_RT$eff <- rownames(full); noinh_RT$eff <- rownames(noinh)
noex_RT$eff <- rownames(noex); nothres_RT$eff <- rownames(nothres)


combined_RT <- rbind(full_RT,noinh_RT, noex_RT, nothres_RT)

combined_RT$model <- factor(combined_RT$model,
                         levels = c("Full", "Fixed Thresholds",
                                    "No Excitation", "No Inhibition"))

combined_RT$Auto <- "Automation Correct (RT Benefit)"
combined_RT$Auto[grep("fail", combined_RT$eff)] <- "Automation Incorrect (RT Cost)"

combined_RT$eff <- factor(combined_RT$eff, levels=c("H_fail_loss",
                                              "H_nonf_gain",
                                              "L_fail_loss",
                                              "L_nonf_gain"),
                       labels=c("High Reliability", "High Reliability",
                                "Low Reliability", "Low Reliability"))


combined_RT %>% ggplot(aes(y=mean, x=eff)) + geom_point(size=2)+ 
  geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
  geom_point(aes_string(x = 'eff', y= 'data'), pch=21, size=3, colour="black")+
  geom_line(aes(group=1,y=data), linetype=2) +facet_grid(model~Auto) +xlab("")+
  ylab("")

save(combined, combined_RT, file="postexp_summaries_A_lb.RData")

