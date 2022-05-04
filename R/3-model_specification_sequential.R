
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


CA_top_thresholdsmult_simple_model <- model.dmc(
  p.map = list(
    A = "1",B = c("sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
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


pnames <- attr(CA_top_thresholdsmult_simple_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 2.5
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("tb", names(pnames))] <- 1

CA_top_thresholdsmult_simple_p.vector  <-
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


check.p.vector(CA_top_thresholdsmult_simple_p.vector,
               CA_top_thresholdsmult_simple_model)

CA_top_thresholdsmult_simple_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_thresholdsmult_simple_p.vector)),
  p1=CA_top_thresholdsmult_simple_p.vector,
  p2=c(1,0.25,1,rep(1, 22), rep(2, 24)),
  lower=c(0.1, 0.01,0, rep(0, 22), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_top_thresholdsmult_simple_p.vector)-2))
)

CA_top_thresholdsmult_simple_dm <- data.model.dmc(cleandats,
                                                  CA_top_thresholdsmult_simple_model)

CA_top_thresholdsmult_simple_samples <- h.samples.dmc(nmc = 180,
                                                      CA_top_thresholdsmult_simple_p.prior,
                                                      CA_top_thresholdsmult_simple_dm, thin=20)

CA_top_thresholdsmult_simple_samples_A_lb <-
  CA_top_thresholdsmult_simple_samples

save(CA_top_thresholdsmult_simple_samples_A_lb,
     file="CA_top_thresholdsmult_simple_samples_A_lb.RData")
