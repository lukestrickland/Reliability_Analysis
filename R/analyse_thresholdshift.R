load("samples_data/CA_top_thresholdsmult_samples.RData")

autos<- lapply(
  CA_top_thresholdsmult_samples,
  function(x) attr(x, "auto"))

auto_names <- names(autos[!is.na(autos)])

tmp <- 1:24
tmp <- tmp[!tmp==4]
full_balance <- c(tmp, 28) 

converged <- auto_names[auto_names %in% full_balance]

feasible_samples <- CA_top_thresholdsmult_samples[converged]

save(feasible_samples, file="feasible_samples.RData")
