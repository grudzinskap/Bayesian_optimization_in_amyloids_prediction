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

if(Sys.info()["nodename"] %in% c("phobos", "michal-XPS14"))
  pathway <- "/home/michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "tobit" )
  pathway <- "~/Dropbox/amyloid_cv_results/"

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

best_enc <- best_positions[["enc_adj"]]

partition_enc_id <- filter(amyloids, enc_adj %in% best_enc) %>%
  select(partition, enc, enc_adj) %>%
  filter(!duplicated(.)) %>%
  slice(1) %>%
  select(partition, enc) %>%
  unlist


ngrams_best_enc <- lapply(1L:45, function(single_replicate_id) {
  load(paste0(pathway, "results/cv_results_full_", single_replicate_id, "_", 
              partition_enc_id[["partition"]],".Rdata"))
  
  sapply(cv_results[[partition_enc_id[["enc"]]]][[1]], function(single_fold)
    single_fold[[2]]) 
}) %>% 
  unlist %>% 
  table(dnn = "ngram") %>%
  data.frame %>%
  filter(Freq == 225) %>%
  select(ngram) %>%
  unlist %>%
  as.character


lapply(split(decode_ngrams(ngrams_best_enc), nchar(decode_ngrams(ngrams_best_enc)))[-1],
       function(single_ngram_length) {
         sapply(single_ngram_length, function(i)
           sapply(single_ngram_length, function(j)
             hamming.distance(strsplit(i, "")[[1]], strsplit(j, "")[[1]])
           )
         )
       })
