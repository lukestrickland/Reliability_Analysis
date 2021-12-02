library(tidyverse)
library(abind)

source("dmc/dmc.R")
source("dmc/dmc_extras.R")

load("samples_data/CA_sess_samples.RData")

#Find LF blocks that follow HF blocks

pnum_samples <- CA_sess_samples

for (i in 1:length(pnum_samples)){
  pnum_samples[[i]]$data$pnum <- names(pnum_samples)[i]
}

data_list <- lapply(pnum_samples, function(x) x$data)

data_df <- do.call("rbind", data_list)

data_df$block <- as.numeric(factor(data_df$block, levels=c("one", "two", "three"),
                                      labels = c("1", "2", "3")))

carryover_sess_subjs <- data_df %>% group_by(pnum, sess) %>% mutate(
  Lblock=unique(block[cond=="L"]), Hblock=unique(block[cond=="H"]))

sess_subjs_carryover <-  carryover_sess_subjs %>% filter(Lblock>Hblock)
sess_subjs_nocarryover <-  carryover_sess_subjs %>% filter(Hblock>Lblock)

carryover <- distinct_at(sess_subjs_carryover,vars(sess,pnum))
nocarryover <- distinct_at(sess_subjs_nocarryover,vars(sess,pnum))

get_full_theta_list <- function(samples, carryover, relstring) {
  
  
  fullthetas_list <- list()  
  
  for (i in 1:length(samples)) {
    
    carryover_sessions <- carryover %>% filter(pnum==names(samples)[i])
    
    sessions <- as.character(
      carryover_sessions$sess
    )
    
    possible_names <- colnames(
      samples[[i]]$theta
    )
    
    carryover_names <- 
      possible_names[
        grep(paste("[", paste(sessions, collapse=""), "]", sep =""), possible_names)]
    
    carryover_greps <- unique(
      str_replace(
        carryover_names,
        paste("[", paste(sessions, collapse=""), "]", sep =""), 
        paste("[", paste(sessions, collapse=""), "]", sep =""))
    )
    
    #for each grep get the average parameter value across the dimensions for each theta
    #which I will the nlater abind together
    posts <- samples[[i]]$theta
    for(j in 1:length(carryover_greps)){
      avthetas <- apply(
        posts[, grep(carryover_greps[j], colnames(posts)),, drop=FALSE], c(1,3),
        mean
      )
      dim(avthetas) <- c(dim(avthetas)[1], 1, dim(avthetas)[-1])
      
      if(j==1) fullthetas <- avthetas else {
        fullthetas <- abind(fullthetas, avthetas, along=2)
      }
      
    }
    
    new_p_names <- str_replace(carryover_names,
                               paste("[", paste(sessions, collapse=""), "]", sep =""),
                               relstring)
    
    colnames(fullthetas) <- unique(new_p_names)
    
    fullthetas_list[[i]] <- fullthetas
    
    print(i)
    
  }
  fullthetas_list
}



fullthetas_carryover <- get_full_theta_list(CA_sess_samples, carryover, 
                                            relstring="carryover")

fullthetas_nocarryover <- get_full_theta_list(CA_sess_samples, nocarryover, 
                                            relstring="nocarryover")

new_samples <- CA_sess_samples

for(i in 1:length(new_samples)) {
  
  fullthetas_all_i <- abind(fullthetas_carryover[[i]], 
                            fullthetas_nocarryover[[i]], 
                            along=2)
  new_samples[[i]]$theta <- fullthetas_all_i 
  
}

test <- get.msds(new_samples)

msds <- test

Vs <- msds[grep("mean_v", rownames(msds)),]

Vs$Cond <- "Manual" 
Vs$Cond[grep("L", rownames(Vs))] <- "Low"
Vs$Cond[grep("H", rownames(Vs))] <- "High"
Vs$Auto <- "Automation Correct"
Vs$Auto[grep("fail", rownames(Vs))] <- "Automation Incorrect"
Vs$S <- "Conflict"
Vs$S[grep("nn", rownames(Vs))] <- "Non-conflict"
Vs$match <- "Match"
Vs$match[grep("false", rownames(Vs))] <- "Mismatch"
Vs$carryover <- "carryover"
Vs$carryover[grep("nocarryover", rownames(Vs))]<- "nocarryover"


names(Vs)[names(Vs)=="Cond"] <- "Condition"

Vs$exinh <- NA
Vs$exinh[Vs$Auto=="Automation Correct" & Vs$match=="Match"] <- "Excitation"
Vs$exinh[Vs$Auto=="Automation Incorrect" & Vs$match=="Match"] <- "Inhibition"
Vs$exinh[Vs$Auto=="Automation Correct" & Vs$match=="Mismatch"] <- "Inhibition"
Vs$exinh[Vs$Auto=="Automation Incorrect" & Vs$match=="Mismatch"] <- "Excitation"


ggplot(Vs, aes(factor(Auto),M)) + 
  geom_point(stat = "identity",aes(shape=carryover, col=Condition), size=2.5) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.3, col=Condition))+ 
  ylab("Accumulation Rate") + xlab("")+
  facet_grid(S ~ match,scales = "free", space = "free") +
  theme(text = element_text(size = 12))

