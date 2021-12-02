source("R/functions.R")
source("dmc/dmc.R")
source("dmc/dmc_extras.R")
load_model ("LBA","lba_B.R")
load("samples_data/CA_top_samples.RData")
theme_set(theme_simple())
load("samples_data/CA_top_samples_pp.RData")

options(dplyr.summarise.inform = FALSE)

pnames <- colnames(CA_top_samples[[1]]$theta)

thres_pickps_set <- rep(
  pnames[
   grep("B\\.M", pnames)
  ], 2
)

thres_pickps_other <- c(
  pnames[grep("B\\.L", pnames)],
  pnames[grep("B\\.H", pnames)]
  )


no_thres <- pickps.h.post.predict.dmc(CA_top_samples, 
                                      save.simulation = TRUE,
                                      pickps_set=thres_pickps_set, 
                                      pickps_other=thres_pickps_other)

save(no_thres, file= "samples_data/no_thres.RData")


ex <- c("mean_v.nn.M.nonf.true", "mean_v.cc.M.nonf.true",
                   "mean_v.nn.M.fail.false", "mean_v.cc.M.fail.false")

ex_pickps_set <- rep(ex, 2)

ex_pickps_other <- c(str_replace(ex, "M", "L"), str_replace(ex, "M", "H"))

no_ex <- pickps.h.post.predict.dmc(CA_top_samples, 
                                      save.simulation = TRUE,
                                      pickps_set=ex_pickps_set, 
                                      pickps_other=ex_pickps_other)

save(no_ex, file= "samples_data/no_ex.RData")


inh <- c("mean_v.nn.M.nonf.false", "mean_v.cc.M.nonf.false",
        "mean_v.nn.M.fail.true", "mean_v.cc.M.fail.true")

inh_pickps_set <- rep(inh, 2)

inh_pickps_other <- c(str_replace(inh, "M", "L"), str_replace(inh, "M", "H"))

no_inh <- pickps.h.post.predict.dmc(CA_top_samples, 
                                   save.simulation = TRUE,
                                   pickps_set=inh_pickps_set, 
                                   pickps_other=inh_pickps_other)

save(no_inh, file= "samples_data/no_inh.RData")





acc_effects <- function (currentsim) {
  
  H_nonf_gain <- NA; L_nonf_gain <- NA;
  H_fail_loss <- NA; L_fail_loss <- NA;
  
  acc_tibble <- currentsim %>% group_by(failtrial, cond) %>% 
    summarise(acc= mean(toupper(substr(S,1,1))==R), war)
  
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

full$model <- "full"; noinh$model <- "noinh"; noex$model <- "noex"
full$eff <- rownames(full); noinh$eff <- rownames(noinh); noex$eff <- rownames(noex)

combined <- rbind(full,noinh, noex)

theme_set(theme_classic())

combined %>% ggplot(aes(y=mean, x=eff)) + geom_point(aes(shape=model), size=2)
