#This Script sets up the final reported model set for the paper which are then
# dispatched on a grid system using pbs pro (see grid_dispatch.R)

##Create model-ready data frame and load dmc functions

#load dmc
source("dmc/dmc.R")

create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                           labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")


#Full classic LBA model of the experiment
load_model("LBA", "lba_B.R")

CA_top_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
    sd_v = c("M"), st0 = "1"),
  match.map = list(
    M = list(nn = "N", cc="C")
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(CA_top_model, "p.vector")

pnames[grep("A", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0


CA_top_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    #True rates, use grepl to avoid sd_v
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))]

    )
  ]
  
  
check.p.vector(CA_top_p.vector, CA_top_model)

CA_top_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_p.vector)),
  p1=CA_top_p.vector,                           
  p2=c(1,1,1,rep(1, 18), rep(2, 24)),
  lower=c(0.1, 0,0, rep(0, 18), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_top_p.vector)-2))
)

CA_top_dm <- data.model.dmc(cleandats,
                                   CA_top_model)

CA_top_samples <- h.samples.dmc(nmc = 180,
                                          CA_top_p.prior,
                                          CA_top_dm, thin=20)

save(CA_top_samples, file="CA_top_samples.RData")


##Carryover hypothesis - try trimming the first 100 trials of each block.

cleandats_lesscarryover <- cleandats %>% filter(trialnum>100)

CA__lesscarryover_dm <- data.model.dmc(cleandats_lesscarryover,
                            CA_top_model)

CA_lesscarryover_samples <- h.samples.dmc(nmc = 180,
                                CA_top_p.prior,
                                CA__lesscarryover_dm, thin=20)

save(CA_lesscarryover_samples, file="CA_lesscarryover_samples.RData")



#Threshold shift to approximate quick information boost



#Full classic LBA model of the experiment
load_model("LBA", "lba_B_autothres.R")


tmap <-
  empty.map(list(
         S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
            failtrial=c("nonf", "fail"),
             R = c("N", "C")
  ),
    levels=c(
             "MAN",
             "LnNS","LnNF",
             "LnCS", "LnCF",
              "LcNS","LcNF",
             "LcCS", "LcCF",
             
             "HnNS","HnNF",
             "HnCS", "HnCF",
             "HcNS","HcNF",
             "HcCS", "HcCF"
             
             ))

tmap[1:72] <- c(
  "HnNS","HcNS","LnNS","LcNS","MAN", "MAN",
  "HnNS","HcNS","LnNS","LcNS","MAN", "MAN",
  "HnNS","HcNS","LnNS","LcNS","MAN", "MAN",
  
  "HnNF","HcNF","LnNF","LcNF","MAN", "MAN",
  "HnNF","HcNF","LnNF","LcNF","MAN", "MAN",
  "HnNF","HcNF","LnNF","LcNF","MAN", "MAN",

  "HnCS","HcCS","LnCS","LcCS","MAN", "MAN",
  "HnCS","HcCS","LnCS","LcCS","MAN", "MAN",
  "HnCS","HcCS","LnCS","LcCS","MAN", "MAN",

  "HnCF","HcCF","LnCF","LcCF","MAN", "MAN",
  "HnCF","HcCF","LnCF","LcCF","MAN", "MAN",
  "HnCF","HcCF","LnCF","LcCF","MAN", "MAN"

)

#checks on match.map

tmap[grepl("M", names(tmap)) ]

tmap[
  grepl("nn", names(tmap)) & grepl("nonf", names(tmap)) & 
    grepl("L", names(tmap)) 
]

tmap[
  grepl("nn", names(tmap)) & grepl("nonf", names(tmap)) & 
    grepl("H", names(tmap)) 
]

tmap[
  grepl("cc", names(tmap)) & grepl("nonf", names(tmap)) &
    grepl("L", names(tmap)) 
]

tmap[
  grepl("cc", names(tmap)) & grepl("nonf", names(tmap)) &
    grepl("H", names(tmap)) 
]

tmap[
  grepl("nn", names(tmap)) & grepl("fail", names(tmap)) &
    grepl("L", names(tmap)) 
]

tmap[
  grepl("nn", names(tmap)) & grepl("fail", names(tmap)) &
    grepl("H", names(tmap)) 
]

tmap[
  grepl("cc", names(tmap)) & grepl("fail", names(tmap)) &
    grepl("L", names(tmap)) 
]

tmap[
  grepl("cc", names(tmap)) & grepl("fail", names(tmap)) &
    grepl("H", names(tmap)) 
]


CA_top_thresholdsmult_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
    sd_v = c("M"), st0 = "1", tb=c("TMAP")),
  match.map = list(
    M = list(nn = "N", cc="C"),
    TMAP=tmap
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1, tb.MAN=1
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(CA_top_thresholdsmult_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("tb", names(pnames))] <- 1

CA_top_thresholdsmult_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A$", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    names(pnames)[grep("tb", names(pnames))],
    #True rates, use grepl to avoid sd_v
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))]
    
  )
  ]


check.p.vector(CA_top_thresholdsmult_p.vector, CA_top_thresholdsmult_model)

CA_top_thresholdsmult_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_thresholdsmult_p.vector)),
  p1=CA_top_thresholdsmult_p.vector,                           
  p2=c(1,1,1,rep(1, 34), rep(2, 24)),
  lower=c(0.1, 0,0, rep(0, 34), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_top_thresholdsmult_p.vector)-2))
)

CA_top_thresholdsmult_dm <- data.model.dmc(cleandats,
                            CA_top_thresholdsmult_model)

CA_top_thresholdsmult_samples <- h.samples.dmc(nmc = 180,
                                CA_top_thresholdsmult_p.prior,
                                CA_top_thresholdsmult_dm, thin=20)

save(CA_top_thresholdsmult_samples, file="CA_top_thresholdsmult_samples.RData")



<<<<<<< HEAD




=======
# # Reviewer comment: Inhibition/excitation mechanisms by day (practice effects?)
# 
# CA_sess_model <- model.dmc(
#   p.map = list(
#     A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "sess", "failtrial", "M"),
#     sd_v = c("M"), st0 = "1"),
#   match.map = list(
#     M = list(nn = "N", cc="C")
#   ),
#   factors = list(
#     S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
#     failtrial=c("nonf", "fail")
#   ),
#   constants = c(st0 = 0, sd_v.false = 1
#   ),
#   responses = c("N", "C"),type = "norm"
# )
# 
# 
# CA_sess_p.vector  <- c(t0=0.3,A=3,
#                                 sd_v.true = 1,
#                
#   B.A.1.N=2, B.M.1.N=2,             
#    B.A.2.N=2, B.M.2.N=2, B.A.1.C=2,             
#   B.M.1.C=2, B.A.2.C=2,   B.M.2.C=2, 
#   
#   mean_v.nn.A.1.nonf.true=1,  mean_v.cc.A.1.nonf.true=1, 
#  mean_v.nn.M.1.nonf.true=1,  mean_v.cc.M.1.nonf.true=1,
#  mean_v.nn.A.1.fail.true=1, mean_v.cc.A.1.fail.true=1, 
#  mean_v.nn.M.1.fail.true=1,  mean_v.cc.M.1.fail.true=1, 
#  mean_v.nn.A.1.nonf.false=0, mean_v.cc.A.1.nonf.false=0,
#  mean_v.nn.M.1.nonf.false=0,mean_v.cc.M.1.nonf.false=0, 
#  mean_v.nn.A.1.fail.false=0, mean_v.cc.A.1.fail.false=0,
#  mean_v.nn.M.1.fail.false=0, mean_v.cc.M.1.fail.false=0,
#  
#  
#  mean_v.nn.A.2.nonf.true=1,  mean_v.cc.A.2.nonf.true=1, 
#  mean_v.nn.M.2.nonf.true=1,  mean_v.cc.M.2.nonf.true=1,
#  mean_v.nn.A.2.fail.true=1, mean_v.cc.A.2.fail.true=1, 
#  mean_v.nn.M.2.fail.true=1,  mean_v.cc.M.2.fail.true=1, 
#  mean_v.nn.A.2.nonf.false=0, mean_v.cc.A.2.nonf.false=0,
#  mean_v.nn.M.2.nonf.false=0,mean_v.cc.M.2.nonf.false=0, 
#  mean_v.nn.A.2.fail.false=0, mean_v.cc.A.2.fail.false=0,
#  mean_v.nn.M.2.fail.false=0, mean_v.cc.M.2.fail.false=0
#  
#  
#  )
# 
# check.p.vector(CA_sess_p.vector, CA_sess_model)
# 
# CA_sess_p.prior <- prior.p.dmc(
#   dists = rep("tnorm", length(CA_sess_p.vector)),
#   p1=CA_sess_p.vector,                           
#   p2=c(1,1,1,rep(1, 8), rep(2, 32)),
#   lower=c(0.1, 0,0, rep(0, 8), rep(NA, 32)),
#   upper=c(5,10, rep(Inf, length(CA_sess_p.vector)-2))
# )
# 
# CA_sess_dm <- data.model.dmc(cleandats,
#                                    CA_sess_model)
# 
# CA_sess_samples <- h.samples.dmc(nmc = 180,
#                                           CA_sess_p.prior,
#                                           CA_sess_dm, thin=20)
# 
# save(CA_sess_samples, file="CA_sess_samples.RData")
# 
# # Reviewer comment: Threshold shifts to approximate a quick information boost
# 
# # Thresholds can be affected by automation via a threshold multiplier
# # use of multiplier assures we can keep the constraint that b>A
# 
# load_model("LBA", "lba_B_autothres.R")
# 
# tmap <-
#   empty.map(list(
#          S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
#             failtrial=c("nonf", "fail"),
#              R = c("N", "C")
#   ), 
#     levels=c(
#              "MAN",
#              "AnNS","AnNF",
#              "AnCS", "AnCF",
#               "AcNS","AcNF",
#              "AcCS", "AcCF"
#              ))
# 
# tmap[1:32] <- c(
#   "AnNS","AcNS","MAN", "MAN",
#   "AnNS","AcNS","MAN", "MAN",
#   
#   "AnNF","AcNF","MAN", "MAN",
#   "AnNF","AcNF","MAN", "MAN",
#   
#   "AnCS","AcCS","MAN", "MAN",
#   "AnCS","AcCS","MAN", "MAN",
#   
#   "AnCF","AcCF","MAN", "MAN",
#   "AnCF","AcCF","MAN", "MAN"
#   
# )
# 
# 
# CA_top_thresholdsmult_model <- model.dmc(
#   p.map = list(
#     A = "1",B = c("sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
#     sd_v = c("M"), st0 = "1", tb="TMAP"),
#   match.map = list(
#     M = list(nn = "N", cc="C"),
#     TMAP=tmap
#   ),
#   factors = list(
#     S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
#     failtrial=c("nonf", "fail")
#   ),
#   constants = c(st0 = 0, sd_v.false = 1, tb.MAN=1
#   ),
#   responses = c("N", "C"),type = "norm"
# )
# 
# 
# CA_top_thresholdsmult_p.vector  <- c(t0=0.3,A=3,
#                                 sd_v.true = 1,
#                
#   B.1.N=2,B.2.N=2, B.1.C=2, B.2.C=2, 
#           
#  tb.AnNS=1,                tb.AnNF=1,                tb.AnCS=1,               
# tb.AnCF=1,                tb.AcNS=1,                tb.AcNF=1,               
#  tb.AcCS=1,                tb.AcCF=1,   
#   
#   mean_v.nn.A.nonf.true=1,  mean_v.cc.A.nonf.true=1, 
#  mean_v.nn.M.nonf.true=1,  mean_v.cc.M.nonf.true=1,
#  mean_v.nn.A.fail.true=1, mean_v.cc.A.fail.true=1, 
#  mean_v.nn.M.fail.true=1,  mean_v.cc.M.fail.true=1, 
#  mean_v.nn.A.nonf.false=0, mean_v.cc.A.nonf.false=0,
#  mean_v.nn.M.nonf.false=0,mean_v.cc.M.nonf.false=0, 
#  mean_v.nn.A.fail.false=0, mean_v.cc.A.fail.false=0,
#  mean_v.nn.M.fail.false=0, mean_v.cc.M.fail.false=0
#  )
# 
# check.p.vector(CA_top_thresholdsmult_p.vector, CA_top_thresholdsmult_model)
# 
# CA_top_thresholdsmult_p.prior <- prior.p.dmc(
#   dists = rep("tnorm", length(CA_top_thresholdsmult_p.vector)),
#   p1=CA_top_thresholdsmult_p.vector,                           
#   p2=c(1,1,1,rep(1, 12), rep(2, 16)),
#   lower=c(0.1, 0,0, rep(0, 12), rep(NA, 16)),
#   upper=c(5,10, rep(Inf, length(CA_top_thresholdsmult_p.vector)-2))
# )
# 
# CA_top_thresholdsmult_dm <- data.model.dmc(cleandats,
#                                    CA_top_thresholdsmult_model)
# 
# CA_top_thresholdsmult_samples <- h.samples.dmc(nmc = 180,
#                                           CA_top_thresholdsmult_p.prior,
#                                           CA_top_thresholdsmult_dm, thin=20)
# 
# save(CA_top_thresholdsmult_samples, file="CA_top_thresholdsmult_samples.RData")
# 
# 
# 
# 
# 
# 
# 
# 

rm(list=ls())

source("dmc/dmc.R")

create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                        labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")


#Full classic LBA model of the experiment
load_model("LBA", "lba_B.R")

CA_fixed_thresholds_model <- model.dmc(
  p.map = list(
    A = "1",B = c("sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
    sd_v = c("M"), st0 = "1"),
  match.map = list(
    M = list(nn = "N", cc="C")
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(CA_fixed_thresholds_model, "p.vector")

pnames[grep("A", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0


CA_fixed_thresholds_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    #True rates, use grepl to avoid sd_v
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))]
    
  )
  ]


check.p.vector(CA_fixed_thresholds_p.vector, CA_fixed_thresholds_model)

CA_fixed_thresholds_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_fixed_thresholds_p.vector)),
  p1=CA_fixed_thresholds_p.vector,                           
  p2=c(1,1,1,rep(1, 6), rep(2, 24)),
  lower=c(0.1, 0,0, rep(0, 6), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_fixed_thresholds_p.vector)-2))
)

CA_fixed_thresholds_dm <- data.model.dmc(cleandats,
                            CA_fixed_thresholds_model)

CA_fixed_thresholds_samples <- h.samples.dmc(nmc = 180,
                                CA_fixed_thresholds_p.prior,
                                CA_fixed_thresholds_dm, thin=20)

save(CA_fixed_thresholds_samples, file="CA_fixed_thresholds_samples.RData")







rm(list=ls())
>>>>>>> 7b3829e497183ee3df41face34b70f3ebd63b375

source("dmc/dmc.R")

create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                        labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")


#Full classic LBA model of the experiment
load_model("LBA", "lba_B.R")

<<<<<<< HEAD
CA_sess_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "sess", "failtrial", "M"),
=======
CA_fixed_thresholds_model <- model.dmc(
  p.map = list(
    A = "1",B = c("sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
>>>>>>> 7b3829e497183ee3df41face34b70f3ebd63b375
    sd_v = c("M"), st0 = "1"),
  match.map = list(
    M = list(nn = "N", cc="C")
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1
  ),
  responses = c("N", "C"),type = "norm"
)


<<<<<<< HEAD
pnames <- attr(CA_sess_model, "p.vector")
=======
pnames <- attr(CA_fixed_thresholds_model, "p.vector")
>>>>>>> 7b3829e497183ee3df41face34b70f3ebd63b375

pnames[grep("A", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0


<<<<<<< HEAD
CA_sess_p.vector  <- 
=======
CA_fixed_thresholds_p.vector  <- 
>>>>>>> 7b3829e497183ee3df41face34b70f3ebd63b375
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    #True rates, use grepl to avoid sd_v
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))]
    
  )
  ]


<<<<<<< HEAD
check.p.vector(CA_sess_p.vector, CA_sess_model)

CA_sess_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_sess_p.vector)),
  p1=CA_sess_p.vector,                           
  p2=c(1,1,1,rep(1, 18), rep(2, 72)),
  lower=c(0.1, 0,0, rep(0, 18), rep(NA, 72)),
  upper=c(5,10, rep(Inf, length(CA_sess_p.vector)-2))
)

CA_sess_dm <- data.model.dmc(cleandats,
                            CA_sess_model)

CA_sess_samples <- h.samples.dmc(nmc = 180,
                                CA_sess_p.prior,
                                CA_sess_dm, thin=20)

save(CA_sess_samples, file="CA_sess_samples.RData")
=======
check.p.vector(CA_fixed_thresholds_p.vector, CA_fixed_thresholds_model)

CA_fixed_thresholds_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_fixed_thresholds_p.vector)),
  p1=CA_fixed_thresholds_p.vector,                           
  p2=c(1,1,1,rep(1, 6), rep(2, 24)),
  lower=c(0.1, 0,0, rep(0, 6), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_fixed_thresholds_p.vector)-2))
)

CA_fixed_thresholds_dm <- data.model.dmc(cleandats,
                                         CA_fixed_thresholds_model)

CA_fixed_thresholds_samples <- h.samples.dmc(nmc = 180,
                                             CA_fixed_thresholds_p.prior,
                                             CA_fixed_thresholds_dm, thin=20)

save(CA_fixed_thresholds_samples, file="CA_fixed_thresholds_samples.RData")






rm(list=ls())


create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                        labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")

source("dmc/dmc.R")

load_model("LBA", "lba_B_automation_diff.R")

mapauto <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ), 
  levels=c(
    "man",
    "HanT","HanF",
    "HacT", "HacF",
    "LanT","LanF",
    "LacT", "LacF"))

mapauto[1:72] <- c(
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man"
  
)

#to make it less constrained (basically the equivalent of freely estimated accumulation rates) I would be able to separate
# auto evidence for failures and successes - psychologically implausible?

#A few checks on mapmeanv due to mind-bending twisting of levels

mapauto[
  grepl("M", names(mapauto)) ]

mapauto[
  grepl("nn", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("nn", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

#parameter to index differences across auto success/failure trials
# in inhibition and excitation

mapautodiff <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ), 
  c(
    "none",
    "hanT","hanF",
    "hacT", "hacF",
    "lanT","lanF",
    "lacT", "lacF"))

mapautodiff[1:72] <- c(
  "hanT","hacF","lanT","lacF","none", "none",
  "hanT","hacF","lanT","lacF","none", "none",
  "hanT","hacF","lanT","lacF","none", "none",
  
  "hacF","hanT","lacF","lanT","none", "none",
  "hacF","hanT","lacF","lanT","none", "none",
  "hacF","hanT","lacF","lanT","none", "none",
  
  "hanF","hacT","lanF","lacT","none", "none",
  "hanF","hacT","lanF","lacT","none", "none",
  "hanF","hacT","lanF","lacT","none", "none",
  
  "hacT","hanF","lacT","lanF","none", "none",
  "hacT","hanF","lacT","lanF","none", "none",
  "hacT","hanF","lacT","lanF","none", "none"
  
)


mapautodiff[grep("M", names(mapautodiff))] <- "none"
mapautodiff[grep("nonf", names(mapautodiff))] <- "none"




mapautodiff[
  grepl("M", names(mapautodiff)) ]

mapautodiff[
  grepl("nn", names(mapautodiff)) & grepl("nonf", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("cc", names(mapautodiff)) & grepl("nonf", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("nn", names(mapautodiff)) & grepl("fail", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("cc", names(mapautodiff)) & grepl("fail", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]



auto_top_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "M"),
    sd_v = c("M"), st0 = "1", a= c("MAPAUTO"), diff= c("MAPAUTODIFF")),
  match.map = list(
    M = list(nn = "N", cc="C"),
    MAPAUTO = mapauto,
    MAPAUTODIFF = mapautodiff
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1, a.man=0, diff.none=0
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(auto_top_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("a\\.", names(pnames))] <- 0
pnames[grep("diff", names(pnames))] <- 0


auto_top_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A$", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))],
    names(pnames)[grep("a\\.", names(pnames))],
    names(pnames)[grep("diff", names(pnames))]
    
  )
  ]


check.p.vector(auto_top_p.vector, auto_top_model)

auto_top_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(auto_top_p.vector)),
  p1=auto_top_p.vector,                           
  p2=c(1,1,1,rep(1, 18), rep(2, 20)),
  lower=c(0.1, 0,0, rep(0, 18), rep(NA, 20)),
  upper=c(5,10, rep(Inf, length(auto_top_p.vector)-2))
)

auto_top_dm <- data.model.dmc(cleandats,
                              auto_top_model)

auto_top_samples <- h.samples.dmc(nmc = 180,
                                  auto_top_p.prior,
                                  auto_top_dm, thin=20)

save(auto_top_samples, file="auto_top_samples.RData")




rm(list=ls())


create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                        labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")

source("dmc/dmc.R")

load_model("LBA", "lba_B_automation_diff.R")

mapauto <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ), 
  levels=c(
    "man",
    "HanT","HanF",
    "HacT", "HacF",
    "LanT","LanF",
    "LacT", "LacF"))

mapauto[1:72] <- c(
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man"
  
)

#to make it less constrained (basically the equivalent of freely estimated accumulation rates) I would be able to separate
# auto evidence for failures and successes - psychologically implausible?

#A few checks on mapmeanv due to mind-bending twisting of levels

mapauto[
  grepl("M", names(mapauto)) ]

mapauto[
  grepl("nn", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("nn", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

#parameter to index differences across auto success/failure trials
# in inhibition and excitation

mapautodiff <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ), 
  c(
    "none",
    "hanT","hanF",
    "hacT", "hacF",
    "lanT","lanF",
    "lacT", "lacF"))

mapautodiff[1:72] <- c(
  "hanT","hacF","lanT","lacF","none", "none",
  "hanT","hacF","lanT","lacF","none", "none",
  "hanT","hacF","lanT","lacF","none", "none",
  
  "hacF","hanT","lacF","lanT","none", "none",
  "hacF","hanT","lacF","lanT","none", "none",
  "hacF","hanT","lacF","lanT","none", "none",
  
  "hanF","hacT","lanF","lacT","none", "none",
  "hanF","hacT","lanF","lacT","none", "none",
  "hanF","hacT","lanF","lacT","none", "none",
  
  "hacT","hanF","lacT","lanF","none", "none",
  "hacT","hanF","lacT","lanF","none", "none",
  "hacT","hanF","lacT","lanF","none", "none"
  
)


mapautodiff[grep("M", names(mapautodiff))] <- "none"
mapautodiff[grep("nonf", names(mapautodiff))] <- "none"




mapautodiff[
  grepl("M", names(mapautodiff)) ]

mapautodiff[
  grepl("nn", names(mapautodiff)) & grepl("nonf", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("cc", names(mapautodiff)) & grepl("nonf", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("nn", names(mapautodiff)) & grepl("fail", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]

mapautodiff[
  grepl("cc", names(mapautodiff)) & grepl("fail", names(mapautodiff)) &
    !grepl("M", names(mapautodiff)) 
]



auto_noexL_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "M"),
    sd_v = c("M"), st0 = "1", a= c("MAPAUTO"), diff= c("MAPAUTODIFF")),
  match.map = list(
    M = list(nn = "N", cc="C"),
    MAPAUTO = mapauto,
    MAPAUTODIFF = mapautodiff
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1, a.man=0, diff.none=0,
                a.LacT = 0, a.LanT= 0, diff.lacT=0, diff.lanT =0
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(auto_noexL_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("a\\.", names(pnames))] <- 0
pnames[grep("diff", names(pnames))] <- 0


auto_noexL_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A$", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))],
    names(pnames)[grep("a\\.", names(pnames))],
    names(pnames)[grep("diff", names(pnames))]
    
  )
  ]


check.p.vector(auto_noexL_p.vector, auto_noexL_model)

auto_noexL_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(auto_noexL_p.vector)),
  p1=auto_noexL_p.vector,                           
  p2=c(1,1,1,rep(1, 18), rep(2, 16)),
  lower=c(0.1, 0,0, rep(0, 18), rep(NA, 16)),
  upper=c(5,10, rep(Inf, length(auto_noexL_p.vector)-2))
)

auto_noexL_dm <- data.model.dmc(cleandats,
                              auto_noexL_model)

auto_noexL_samples <- h.samples.dmc(nmc = 180,
                                  auto_noexL_p.prior,
                                  auto_noexL_dm, thin=20)

save(auto_noexL_samples, file="auto_noexL_samples.RData")





rm(list=ls())


create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO_H","AUTO_L", "MANUAL"),
                           labels=c("H", "L", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                        labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")

source("dmc/dmc.R")

load_model("LBA", "lba_B_automation.R")

mapauto <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ), 
  levels=c(
    "man",
    "HanT","HanF",
    "HacT", "HacF",
    "LanT","LanF",
    "LacT", "LacF"))

mapauto[1:72] <- c(
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  "HanT","HacF","LanT","LacF","man", "man",
  
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  "HacF","HanT","LacF","LanT","man", "man",
  
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  "HanF","HacT","LanF","LacT","man", "man",
  
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man",
  "HacT","HanF","LacT","LanF","man", "man"
  
)

#to make it less constrained (basically the equivalent of freely estimated accumulation rates) I would be able to separate
# auto evidence for failures and successes - psychologically implausible?

#A few checks on mapmeanv due to mind-bending twisting of levels

mapauto[
  grepl("M", names(mapauto)) ]

mapauto[
  grepl("nn", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("nonf", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("nn", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]

mapauto[
  grepl("cc", names(mapauto)) & grepl("fail", names(mapauto)) &
    !grepl("M", names(mapauto)) 
]




auto_sym_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "M"),
    sd_v = c("M"), st0 = "1", a= c("MAPAUTO")),
  match.map = list(
    M = list(nn = "N", cc="C"),
    MAPAUTO = mapauto
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1, a.man=0
  ),
  responses = c("N", "C"),type = "norm"
)


pnames <- attr(auto_sym_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 3
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("a\\.", names(pnames))] <- 0



auto_sym_p.vector  <- 
  pnames[c(
    names(pnames)[grep("t0", names(pnames))],
    names(pnames)[grep("A$", names(pnames))],
    names(pnames)[grep("sd_v", names(pnames))],
    names(pnames)[grep("B", names(pnames))],
    names(pnames)[grep("true", names(pnames))][
      !grepl("sd_v", names(pnames)[grep("true", names(pnames))])
    ],
    names(pnames)[grep("false", names(pnames))],
    names(pnames)[grep("a\\.", names(pnames))]
    
  )
  ]


check.p.vector(auto_sym_p.vector, auto_sym_model)

auto_sym_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(auto_sym_p.vector)),
  p1=auto_sym_p.vector,                           
  p2=c(1,1,1,rep(1, 18), rep(2, 12)),
  lower=c(0.1, 0,0, rep(0, 18), rep(NA, 12)),
  upper=c(5,10, rep(Inf, length(auto_sym_p.vector)-2))
)

auto_sym_dm <- data.model.dmc(cleandats,
                              auto_sym_model)

auto_sym_samples <- h.samples.dmc(nmc = 180,
                                  auto_sym_p.prior,
                                  auto_sym_dm, thin=20)

save(auto_sym_samples, file="auto_sym_samples.RData")


>>>>>>> 7b3829e497183ee3df41face34b70f3ebd63b375



