###################################################################################################
################################# ATC-LAB DATA EXTRACTION ROUTINE #################################
###################################################################################################

### clear workspace ###
rm(list = ls())

#requires tidyverse
source('R/0-extract_functions.R')

xml_files <- dir("XML_logs")

for(i in 1:length(xml_files)){
    splitstrings <- strsplit(xml_files[i], split="_")[[1]]
    
    #Fix up for new automation condition names
    if(length(splitstrings)>5) 
      splitstrings[5] <- paste(
        splitstrings[5], splitstrings[6], sep="_")
    ppt <- gsub("ppt", "", splitstrings[1])
    sess <- gsub("sess", "", splitstrings[2])
    block <- gsub("block", "", splitstrings[3])
    cond <- gsub("\\.xml\\.log", "", splitstrings[5])
    current_df <- extract_one(ppt=ppt, sess=sess, block= block,  cond=cond)
    print(i)
    if(i==1) full_df <- current_df else full_df <- rbind(full_df, current_df)
}

dats <- full_df[, c("ppt", "sess", "block", "cond", "stimulus",
                    "failtrial",
                    "response", "rkey", "RT", "score", 
                    "cumulative_score")]

#Get responses from Rkey
dats$rkey[!dats$rkey %in% c("F", "J")] <- NA

#convert response key to responses
dats$R <- NA
key_set <- c(rep(1:2, 6), rep(2:1, 6))


for (i in unique(dats$ppt)) {
  counterbalance <- (as.numeric(i) %% 24) + 1
  cb_key <- key_set[counterbalance]
  if (cb_key==1) {
    dats$R[dats$rkey=="F" & dats$ppt==i] <- "C"
    dats$R[dats$rkey=="J"& dats$ppt==i] <- "N" }
  else {
    dats$R[dats$rkey=="F"& dats$ppt==i] <- "N"
    dats$R[dats$rkey=="J"& dats$ppt==i] <- "C"     
  }
  
}

dats$stimulus <- factor(dats$stimulus, levels = c("nonconflict", "conflict"),
                        labels=c("n", "c"))

save(dats, file="img/dats_raw.RData")

load("img/dats_raw.RData")

#DATA EXCLUSIONS - NON RESPONSES AND <200MS
100* length(dats$RT[is.na(dats$RT)])/length(dats$RT) 
cleandats <- dats[!is.na(dats$RT),]

oldlen <- length(cleandats$RT)
cleandats <-  cleandats[cleandats$RT>200,]

(oldlen -length(cleandats$RT))/length(cleandats$RT) *100

###CLEAN UP DATA FRAME

cleandats$cond <- factor(cleandats$cond)
cleandats$failtrial <- factor(cleandats$failtrial,
                              levels=c("FALSE", "TRUE"),
                              labels=c("nonf", "fail"))
cleandats$sess <- factor(cleandats$sess)
cleandats$block <- factor(cleandats$block, levels=c("1", "2", "3"),
                          labels= c("one", "two", "three"))

colnames(cleandats)[colnames(cleandats)=='ppt'] <- 's'
colnames(cleandats)[colnames(cleandats)=='stimulus'] <- 'S'
cleandats <- cleandats[,c("s", "sess", "block", "cond",
                          "failtrial", "S", "R", "RT")]

cleandats$RT <- cleandats$RT/1000

save(cleandats, file="img/cleandats.RData")



###A few checks that the data from the csvs with map/aircraft data
###align ok with data from xml logs


all(full_df$stimulus==full_df$conflict_status)
all(full_df$stimulus[full_df$DOMS<5]=="conflict")
all(full_df$stimulus[full_df$DOMS>5]=="nonconflict")

all(full_df$autorec[full_df$cond=="MANUAL"]=="########")

all(full_df$autorec[full_df$stimulus=="conflict" 
                    & full_df$failtrial & full_df$cond=="AUTO"]=="NON-CONF")

all(full_df$autorec[full_df$stimulus=="conflict" 
                    & !full_df$failtrial & full_df$cond=="AUTO"]=="CONFLICT")

all(full_df$autorec[full_df$stimulus=="nonconflict" 
                    & full_df$failtrial & full_df$cond=="AUTO"]=="CONFLICT")

all(full_df$autorec[full_df$stimulus=="nonconflict" 
                    & !full_df$failtrial & full_df$cond=="AUTO"]=="NON-CONF")


all(full_df$stimulus[full_df$response=="FALSE_ALARM_NON_CONFLICT"]=="conflict")
all(full_df$stimulus[full_df$response=="FALSE_ALARM_CONFLICT"]=="nonconflict")

all(full_df$stimulus[full_df$response=="HIT_CONFLICT_ONTIME"]=="conflict")
all(full_df$stimulus[full_df$response=="HIT_NON_CONFLICT"]=="nonconflict")

all(full_df$stimulus[full_df$response=="MISS_CONFLICT"]=="conflict")
all(full_df$stimulus[full_df$response=="MISS_NON_CONFLICT"]=="nonconflict")





dats2 <- full_df[, c("ppt", "sess", "block", "cond", "stimulus",
                    "failtrial",
                    "response", "rkey", "RT", "score", 
                    "cumulative_score", "response")]

#Get responses from Rkey
dats2$rkey[!dats2$rkey %in% c("F", "J")] <- NA

#convert response key to responses
dats2$R <- NA
key_set <- c(rep(1:2, 6), rep(2:1, 6))


for (i in unique(dats2$ppt)) {
  counterbalance <- (as.numeric(i) %% 24) + 1
  cb_key <- key_set[counterbalance]
  if (cb_key==1) {
    dats2$R[dats2$rkey=="F" & dats2$ppt==i] <- "C"
    dats2$R[dats2$rkey=="J"& dats2$ppt==i] <- "N" }
  else {
    dats2$R[dats2$rkey=="F"& dats2$ppt==i] <- "N"
    dats2$R[dats2$rkey=="J"& dats2$ppt==i] <- "C"     
  }
  
}

dats2$stimulus <- factor(dats2$stimulus, levels = c("nonconflict", "conflict"),
                        labels=c("n", "c"))


all(dats2$R[dats2$response=="HIT_CONFLICT_ONTIME"]=="C")
all(dats2$R[dats2$response=="FALSE_ALARM_CONFLICT"]=="C")
all(dats2$R[dats2$response=="FALSE_ALARM_NON_CONFLICT"]=="N")
all(dats2$R[dats2$response=="FALSE_ALARM_NON_CONFLICT"]=="N")

table(dats$sess, dats$ppt, dats$cond, dats$block)


##check trialnums valid

dats <-
  dats %>% group_by(ppt, sess, cond) %>% mutate(trial_counter=1:length(RT))

#check auto failure trials are same
# all(dats$trial_counter[dats$failtrial & dats$cond=="MANUAL"]==
#       dats$trial_counter[dats$failtrial & dats$cond=="AUTO"])
# 
# all(dats$stimulus[dats$failtrial & dats$cond=="MANUAL"]==
#       dats$stimulus[dats$failtrial & dats$cond=="AUTO"])
# 
# dats$DOMS <- full_df$DOMS
# dats$ac1speed <- full_df$ac1_speed
# 
# all(dats$DOMS[dats$failtrial & dats$cond=="MANUAL"]==
#       dats$DOMS[dats$failtrial & dats$cond=="AUTO"])
# 
# all(dats$ac1speed[dats$failtrial & dats$cond=="MANUAL"]==
#       dats$ac1speed[dats$failtrial & dats$cond=="AUTO"])
# 
# dats$stimulus[!dats$failtrial & dats$cond=="MANUAL"]==
#   dats$stimulus[!dats$failtrial & dats$cond=="AUTO"]
# 
# dats$ac1speed[!dats$failtrial & dats$cond=="MANUAL"]==
#   dats$ac1speed[!dats$failtrial & dats$cond=="AUTO"]
# 
# dats$DOMS[!dats$failtrial & dats$cond=="MANUAL"]==
#   dats$DOMS[!dats$failtrial & dats$cond=="AUTO"]

