source("dmc/dmc.R")
source("dmc/dmc_extras.R")
source("R/0-analysis_functions.R")

fitpath <- file.path(set_fit_path(), "Reliability_Analysis")
loadpath <- create_loadpath(fitpath)
savepath <- create_savepath(fitpath)

loadpath("CA_top_thresholdsmult_simple_samples_A_lb.RData")

autos<- sapply(
  CA_top_thresholdsmult_simple_samples_A_lb,
  function(x) attr(x, "auto"))

auto_names <- names(autos[!is.na(autos)])

tmp <- 1:24
tmp <- tmp[!tmp==4]
full_balance <- c(tmp, 28) 

converged <- auto_names[auto_names %in% full_balance]

feasible_samples <- CA_top_thresholdsmult_samples[converged]

save(feasible_samples, file="feasible_samples.RData")
