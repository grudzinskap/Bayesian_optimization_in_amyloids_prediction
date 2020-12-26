#amyloidy

library(dplyr)
if(Sys.info()["nodename"] == "sobczyk-pc" )
  setwd("~/Dropbox/doktorat/moje_prace/amyloid_ngram2/")
amyloids <- read.csv("results/all_summaries.csv")

summary(amyloids)

#what kind of AUC can we expect?
amyloids %>% group_by(len_range) %>% 
  summarise(max(AUC_mean), mean(AUC_mean), min(AUC_mean)) 
#there is a huge gain when we take into account that
#each encoding has different competence regions

#if we would use just one encoding for all then differences in performance are very small
amyloids %>% group_by(enc) %>% 
  summarise(AUC_mean=mean(AUC_mean)) %>%
  summarise(min(AUC_mean), mean(AUC_mean), max(AUC_mean))

#even if we take into account non-homogenity of learning set
amyloids %>% group_by(enc, pos) %>% 
  summarise(AUC_mean=mean(AUC_mean)) -> AUC_means
summary(AUC_means$AUC_mean)

#if we learn on 15-s then we get highest sensitivity. Naturally we pay with specificity
amyloids %>% group_by(len_range, pos) %>%
  summarise(max(Sens_mean), mean(Sens_mean), min(Sens_mean),
            max(Spec_mean), mean(Spec_mean), min(Spec_mean)) 


ggplot(amyloids, aes(x = len_range, y = Spec_mean)) +
  geom_boxplot() +
  facet_wrap(~ pos, labeller = label_both) +
  ggtitle("Specificity")

png("sesp_plot.png", height = 1024, width = 1024)
ggplot(amyloids, aes(x = Sens_mean, y = Spec_mean)) +
  geom_point() +
  facet_grid(len_range~ pos, labeller = label_both) 
dev.off()

#wychodzi na to, że prawdziwa różnica jest na sensitivity. Różnice są naprawdę duże,
#co sugeruje, że warto jest robić osobne modele dla każdego typu długości sekwencji
