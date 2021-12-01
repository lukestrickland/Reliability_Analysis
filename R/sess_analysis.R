library(tidyverse)
library(abind)

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

new_samples <- CA_sess_samples

for (i in 1:length(CA_sess_samples)) {
  
  carryover_sessions <- carryover %>% filter(pnum==names(CA_sess_samples)[i])
  
  sessions <- as.character(
    carryover_sessions$sess
  )
  
  possible_names <- colnames(
    CA_sess_samples[[i]]$theta
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
  posts <- CA_sess_samples[[i]]$theta
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
                             "carryover")
  
  colnames(fullthetas) <- new_p_names
  
  new_samples[[i]]$theta <- fullthetas
  

}