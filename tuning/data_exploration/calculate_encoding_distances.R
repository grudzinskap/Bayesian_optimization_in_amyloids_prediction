#compute encoding distances

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(cvTools)
require(ranger)
require(hmeasure)
require(e1071)
library(foreach)

#### loading data and identifying best encoding #####

amyloids <- read.csv("results/all_summaries.csv")
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

best_positions <- amyloids %>% 
  select(len_range, enc_adj, AUC_mean) %>%
  group_by(len_range) %>%
  mutate(position = rank(AUC_mean)) %>%
  ungroup %>%
  group_by(enc_adj) %>%
  summarise(cum_rank = sum(position)) %>%
  arrange(desc(cum_rank)) %>%
  slice(1L:5)

best_enc <- best_positions[["enc_adj"]][1]

##### calculating encoding distances ######

# very time consuming, load from file instead
distances <- foreach(enc = aa_groups) %dopar% calc_ed(aa_groups[[best_enc]], enc)

# distances = lapply(aa_groups, function(enc) calc_ed(aa_groups[[best_enc]], enc))
save(distances, file="results/distances2.Rdata")
# load(file = "results/distances.Rdata")