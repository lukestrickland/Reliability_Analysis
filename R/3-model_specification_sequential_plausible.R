
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


#nS = decision aid recommends N cF = decision aid recommends N
#nF = decision aid recommends C cS = decision aid recommends C

tmap <-
  empty.map(list(
    S = c("nn", "cc"), cond = c("H", "L", "M"), sess = c("1", "2", "3"),
    failtrial=c("nonf", "fail"),
    R = c("N", "C")
  ),
  levels=c(
    "MAN",
    "L_dn_N","L_dc_N",
    "L_dn_C", "L_dc_C",
    
    "H_dn_N","H_dc_N",
    "H_dn_C", "H_dc_C"
    
  ))

tmap[1:72] <- c(
  "H_dn_N","H_dc_N","L_dn_N","L_dc_N","MAN", "MAN",
  "H_dn_N","H_dc_N","L_dn_N","L_dc_N","MAN", "MAN",
  "H_dn_N","H_dc_N","L_dn_N","L_dc_N","MAN", "MAN",
  
  "H_dc_N","H_dn_N","L_dc_N","L_dn_N","MAN", "MAN",
  "H_dc_N","H_dn_N","L_dc_N","L_dn_N","MAN", "MAN",
  "H_dc_N","H_dn_N","L_dc_N","L_dn_N","MAN", "MAN",
  
  "H_dn_C","H_dc_C","L_dn_C","L_dc_C","MAN", "MAN",
  "H_dn_C","H_dc_C","L_dn_C","L_dc_C","MAN", "MAN",
  "H_dn_C","H_dc_C","L_dn_C","L_dc_C","MAN", "MAN",
  
  "H_dc_C","H_dn_C","L_dc_C","L_dn_C","MAN", "MAN",
  "H_dc_C","H_dn_C","L_dc_C","L_dn_C","MAN", "MAN",
  "H_dc_C","H_dn_C","L_dc_C","L_dn_C","MAN", "MAN"
  
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


CA_top_information_boost_model <- model.dmc(
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


pnames <- attr(CA_top_information_boost_model, "p.vector")

pnames[grep("A$", names(pnames))] <- 2.5
pnames[grep("B", names(pnames))] <- 2
pnames[grep("t0", names(pnames))] <- 0.3
pnames[grep("true", names(pnames))] <- 1
pnames[grep("false", names(pnames))] <- 0
pnames[grep("tb", names(pnames))] <- 1

CA_top_information_boost_p.vector  <-
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


check.p.vector(CA_top_information_boost_p.vector,
               CA_top_information_boost_model)

CA_top_information_boost_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_information_boost_p.vector)),
  p1=CA_top_information_boost_p.vector,
  p2=c(1,0.25,1,rep(1, 14), rep(2, 24)),
  lower=c(0.1, 0.01,0, rep(0, 14), rep(NA, 24)),
  upper=c(5,10, rep(Inf, length(CA_top_information_boost_p.vector)-2))
)

CA_top_information_boost_dm <- data.model.dmc(cleandats,
                                                  CA_top_information_boost_model)


#Condition assigments are correct
tmap[grepl("L", names(tmap))]
tmap[grepl("H", names(tmap))]

#Response assigments are correct
tmap[grepl("N", names(tmap))& !grepl("M", names(tmap))]
tmap[grepl("C", names(tmap))& !grepl("M", names(tmap))]

#Manual assigment is correct
tmap[grepl("M", names(tmap))]

#Tricky one - automation recommendation assigment

#automation recommended non-conflict should be dn
tmap[grepl("nonf", names(tmap))& 
       grepl("nn", names(tmap))& 
       !grepl("M", names(tmap))]

tmap[grepl("fail", names(tmap))& 
       grepl("cc", names(tmap))& 
       !grepl("M", names(tmap))]

#automation recommended conflict should be dc
tmap[grepl("nonf", names(tmap))& 
       grepl("cc", names(tmap))& 
       !grepl("M", names(tmap))]

tmap[grepl("fail", names(tmap))& 
       grepl("nn", names(tmap))& 
       !grepl("M", names(tmap))]


CA_top_information_boost_samples <- h.samples.dmc(nmc = 180,
                                                  CA_top_information_boost_p.prior,
                                                  CA_top_information_boost_dm, 
                                                  thin=20)


save(CA_top_information_boost_samples,
     file="CA_top_information_boost_samples.RData")
