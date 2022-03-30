

accs_2 <-
  cleandats %>% group_by(s, Stimulus, Condition) %>% 
  filter(!is.na(R)) %>% summarise(acc = mean(C)) %>%
  arrange(s)


mean_accs <- accs_2 %>% group_by(Stimulus, Condition) %>%
  summarise(meanacc=mean(acc))

searr= se2(accs, facs= c("Stimulus", "Condition"),
           dvnam="acc", sfac="s")
se_accs <- as.data.frame.table(searr)
colnames(se_accs)[colnames(se_accs)=="Freq"] <- "seacc"

plot.df <- full_join(mean_accs, se_accs)



ggplot(plot.df, aes(Condition, meanacc)) +geom_point(aes(col=Stimulus, shape=Stimulus), size=3)  + 
  geom_errorbar(aes(
    ymax = meanacc + seacc,
    ymin = meanacc - seacc,
    colour = Stimulus, width = 0.3
  )) +ylab ("Accuracy") +xlab("") +
  theme(text = element_text(size = 21))