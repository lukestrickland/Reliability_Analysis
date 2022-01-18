
source("dmc/dmc.R")

load_model("LBA", "lba_B.R")
run.grid.dmc("CA_sess_samples",model.dir ="LBA",
             model.file="lba_B.R",user="ljs392",
             n.add=60, wall.hours = 300,
             GB = 5, max.try=10)