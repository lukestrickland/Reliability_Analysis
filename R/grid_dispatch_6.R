rm(list=ls())
source("dmc/dmc.R")

load_model("LBA", "lba_B.R")

run.grid.dmc("CA_top_samples_A",model.dir ="LBA",
             model.file="lba_B.R",user="ljs392",
             n.add=60, wall.hours = 300,
             GB = 3, max.try=5)
# 
# run.grid.dmc("CA_fixed_thresholds_samples",model.dir ="LBA",
#              model.file="lba_B.R",user="ljs392",
#              n.add=60, wall.hours = 300,
#              GB = 3, max.try=5)
# 
# 
# rm(list=ls())
# source("dmc/dmc.R")
# 
# load_model("LBA", "lba_B_automation.R")
# run.grid.dmc("auto_top_samples",model.dir ="LBA",
#              model.file="lba_B_automation.R",user="ljs392",
#              n.add=60, wall.hours = 300,
#              GB = 3, max.try=5)
