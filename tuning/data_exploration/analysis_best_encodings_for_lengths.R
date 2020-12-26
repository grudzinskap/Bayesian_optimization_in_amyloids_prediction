source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")


if(Sys.info()["nodename"] == "tobit" )
  setwd("~/Dropbox/doktorat/moje_prace/amyloid_ngram2/")


library(dplyr)
library(biogram) #calc_ed
require(seqinr) #for choose_properties
require(ggplot2)

# read data ----------------------------------------

amyloids <- read.csv("results/all_summaries.csv")
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)


#############this part is copied from function choose_properties ################
prop_MK <- read.csv2("AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1]) %>%
  as.numeric

prop_MK <- cbind(prop_MK, years = years) %>% filter(years >= 1980) %>%
  select(X, name) %>% rbind(., 
                            data.frame(X = 545:550, 
                                       name =  c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")))
#################################################################################

trait_tab <- read.csv(file = "./results/trait_tab.csv")

# start of analysis ----------------------------------------


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

distances <- sapply(best_enc, function(i)
  sapply(best_enc, function(j)
    calc_ed(aa_groups[[i]], aa_groups[[j]])
  )
)

# properties - everything -----------------------------

load("./results/enc_dupes.RData")
#remember, indices in enc_adj and are equal to indices in enc_dupes + 2 and trait_tab  + 2

best_dupes <- lapply(enc_dupes[["aa_duplicates"]][best_enc], function(single_enc)
  single_enc %>%
    sapply(function(i) substr(i, 3, nchar(i))) %>%
    strsplit("K") %>%
    lapply(as.numeric) %>%
    do.call(rbind, .) %>%
    data.frame %>%
    select(X1) %>%
    unlist(use.names = FALSE) - 2)

frequencies <- data.frame(trait_tab[unlist(best_dupes, use.names = FALSE), ] %>%
                            apply(1, na.omit) %>% unlist %>% table)
colnames(frequencies) <- c("X", "Freq")
frequencies$X <- as.numeric(as.character(frequencies$X))

trait_props_all <- prop_MK %>%
  inner_join(frequencies) %>%
  arrange(desc(Freq)) %>%
  select(name, Freq) %>%
  mutate(prop = Freq/length(unlist(best_dupes, use.names = FALSE))) %>%
  droplevels()

trait_props_all[["name"]] <- factor(trait_props_all[["name"]],
                                    levels = rev(as.character(trait_props_all[["name"]])))

levels(trait_props_all[["name"]]) <- sapply(strwrap(levels(trait_props_all[["name"]]), 30, simplify = FALSE),
       function(i) paste0(i, collapse = "\n"))

ggplot(trait_props_all, aes(x = name, y = prop)) +
  geom_bar(stat = "identity") +
  coord_flip()


